//! Multi-modular reconstruction of a reduced rational F4 basis.
use crate::f4::{F4Trace, PRIMES, groebner_basis_f4_direct, learn};
use crate::{Fp, GroebnerError, Monomial, Polynomial, par};
use num_rational::BigRational;
use polycore::Primes;
use polycore::crt::{self, MixedRadixAccumulator as Accumulator};
use std::collections::HashMap;

type Layout = Vec<Vec<Monomial>>;
// An image, and the trace learned computing it, if any.
type Learned = (Vec<Polynomial<Fp>>, Option<F4Trace>);
// A basis modulo the prime it was computed at.
type Residues = (u64, Vec<Polynomial<Fp>>);

/// Compute a rational basis by CRT and rational reconstruction, checked at a fresh prime.
/// With `certify`, prove equality with the input ideal over Q: check Buchberger's
/// criterion, reduce every input by the candidate, and verify exact cofactor
/// certificates expressing every candidate polynomial in the input generators.
/// Certification builds an auxiliary Buchberger basis with cofactor tracking and
/// can be substantially more expensive than the modular computation.
pub fn groebner_basis_f4_rational(
    polynomials: Vec<Polynomial<BigRational>>,
    certify: bool,
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    groebner_basis_f4_rational_with(polynomials, certify, RationalOptions::default())
}

/// Tuning for [`groebner_basis_f4_rational_with`].
#[derive(Clone, Copy, Debug, Default)]
pub struct RationalOptions {
    /// Primes replayed per round after the first image. By default rounds start at
    /// eight primes and grow by about 12.5%, which suits systems that need few primes
    /// or have large matrices. When the first image is quick, rounds instead fill the
    /// thread pool with one lane chunk per thread, amortizing their overhead.
    pub batch: Option<usize>,
}

/// [`groebner_basis_f4_rational`] with an explicit replay schedule.
pub fn groebner_basis_f4_rational_with(
    polynomials: Vec<Polynomial<BigRational>>,
    certify: bool,
    options: RationalOptions,
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    let input = crate::groebner::prepare_input(polynomials)?;
    let primitive: Vec<_> = input.iter().map(Polynomial::primitive).collect();
    // 23-bit primes allow fully deferred u64 accumulation for matrices with up to
    // 2^18 columns, trading more CRT images for much cheaper row elimination.
    let basis = reconstruct_inner(
        &primitive,
        Primes::below(1 << 23),
        options.batch,
        true,
        |_, _| {},
        &mut ReconstructionStats::default(),
    )?;
    if certify && !certify_basis(&input, &basis)? {
        return Err(GroebnerError::ReconstructionFailed);
    }
    Ok(basis)
}

/// Establish both ideal inclusions independently of modular reconstruction.
fn certify_basis(
    input: &[Polynomial<BigRational>],
    basis: &[Polynomial<BigRational>],
) -> Result<bool, GroebnerError> {
    if !crate::is_groebner_basis(basis)? || input.iter().any(|p| !p.reduce(basis).is_zero()) {
        return Ok(false);
    }
    // Input reduction proves <input> is contained in <basis>. In particular,
    // [1] passes those checks for every input, so the reverse inclusion is vital.
    let lifted = crate::LiftBasis::new(input)?;
    for polynomial in basis {
        let Some(cofactors) = lifted.lift(polynomial)? else {
            return Ok(false);
        };
        if !crate::verify_lift(input, &cofactors, polynomial) {
            return Ok(false);
        }
    }
    Ok(true)
}

// Counts of reconstruction events, which tests check.
#[derive(Default, Debug)]
struct ReconstructionStats {
    images: usize,
    trace_replacements: usize,
    validation_failures: usize,
    recovery_restarts: usize,
}

#[derive(Default)]
struct TraceState {
    current: Option<F4Trace>,
    failed_leading: Option<Vec<Monomial>>,
}
impl TraceState {
    fn observe(&mut self, leading: &[Monomial], learned: Option<F4Trace>) -> bool {
        let Some(learned) = learned else {
            self.failed_leading = None;
            return false;
        };
        if self.current.is_none() {
            self.current = Some(learned);
            return false;
        }
        // Two consecutive full fallbacks with the same final heads corroborate a
        // replacement. One unlucky replay prime must not evict a working trace.
        if self.failed_leading.as_deref() == Some(leading) {
            self.current = Some(learned);
            self.failed_leading = None;
            return true;
        }
        self.failed_leading = Some(leading.to_vec());
        false
    }
}

struct Group {
    layout: Layout,
    primary: Accumulator,
    recovery: Option<Accumulator>,
    recovery_limit: usize,
}
impl Group {
    fn new(size: usize) -> Self {
        Self {
            layout: vec![Vec::new(); size],
            primary: Accumulator::new(0),
            recovery: None,
            recovery_limit: 32,
        }
    }

    fn add(&mut self, p: u64, basis: &[Polynomial<Fp>]) -> Result<bool, GroebnerError> {
        // Extend support without losing earlier images: every newly encountered
        // coefficient was zero modulo all preceding primes in this group.
        let mut positions = Vec::new();
        let mut offset = 0;
        let mut changed = false;
        for (support, f) in self.layout.iter_mut().zip(basis) {
            if support.iter().eq(f.terms.iter().map(|t| &t.0)) {
                positions.extend((offset..offset + support.len()).map(Some));
                offset += support.len();
                continue;
            }
            let old: HashMap<_, _> = support
                .iter()
                .cloned()
                .enumerate()
                .map(|(i, m)| (m, offset + i))
                .collect();
            offset += support.len();
            for (m, _) in &f.terms {
                if !old.contains_key(m) {
                    support.push(m.clone());
                    changed = true;
                }
            }
            support.sort_by(|a, b| f.order.compare(b, a));
            positions.extend(support.iter().map(|m| old.get(m).copied()));
        }
        if changed {
            self.primary.remap(&positions);
            if let Some(recovery) = &mut self.recovery {
                recovery.remap(&positions);
            }
        }
        let mut values = Vec::with_capacity(positions.len());
        for (support, f) in self.layout.iter().zip(basis) {
            let mut terms = f.terms.iter().peekable();
            for m in support {
                if let Some((_, coefficient)) = terms.next_if(|t| &t.0 == m) {
                    values.push(coefficient.value());
                } else {
                    values.push(0);
                }
            }
        }
        let restart = self
            .recovery
            .as_ref()
            .is_some_and(|a| a.image_count() >= self.recovery_limit);
        if restart {
            self.recovery_limit = self.recovery_limit.saturating_mul(2);
            self.recovery = Some(Accumulator::new(values.len()));
        } else if self.recovery.is_none() && self.primary.image_count() >= 32 {
            self.recovery = Some(Accumulator::new(values.len()));
        }
        // Primes come from one descending sequence, so they are distinct words.
        let merge = |acc: &mut Accumulator| {
            acc.add(p, &values)
                .map_err(|_| GroebnerError::ReconstructionFailed)
        };
        merge(&mut self.primary)?;
        if let Some(recovery) = &mut self.recovery {
            merge(recovery)?;
        }
        Ok(restart)
    }

    fn width(&self) -> usize {
        self.layout.iter().map(Vec::len).sum()
    }

    fn decode(
        &self,
        values: Vec<BigRational>,
        input: &Polynomial<BigRational>,
    ) -> Vec<Polynomial<BigRational>> {
        let mut values = values.into_iter();
        self.layout
            .iter()
            .map(|support| {
                let terms = support.iter().cloned().zip(values.by_ref()).collect();
                Polynomial::new(terms, input.nvars, input.order.clone())
            })
            .collect()
    }
}

fn heads(basis: &[Polynomial<Fp>]) -> Vec<Monomial> {
    basis.iter().filter_map(|f| f.lm().cloned()).collect()
}

#[cfg(test)]
fn reconstruct(
    input: &[Polynomial<BigRational>],
    primes: impl IntoIterator<Item = u64>,
    inspect_image: impl FnMut(u64, &mut Vec<Polynomial<Fp>>),
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    reconstruct_inner(
        input,
        primes,
        None,
        false,
        inspect_image,
        &mut ReconstructionStats::default(),
    )
}

// A first round quicker than this learns a trace whose narrow replay rounds would
// take about a third as long, mostly fork-join and reconstruction overhead, so the
// replay rounds widen to one lane chunk per thread. Deciding before any replay avoids
// overshooting the primes needed with a late wide round.
const WIDE_LEARN: std::time::Duration = std::time::Duration::from_millis(150);

// After the learned image, replay whole lane chunks and grow by about 12.5%, rather
// than doubling the work near completion, unless the caller fixed the width or a quick
// first image widened it.
fn schedule(count: usize, batch: Option<usize>, wide: bool) -> usize {
    let narrow = count
        .div_ceil(8)
        .max(2 * PRIMES)
        .next_multiple_of(PRIMES)
        .min(32);
    match (count, batch) {
        (0, _) => 1,
        (_, Some(batch)) => batch.max(1),
        _ if wide => narrow.max(PRIMES * par::threads()),
        _ => narrow,
    }
}

fn reconstruct_inner(
    input: &[Polynomial<BigRational>],
    primes: impl IntoIterator<Item = u64>,
    batch: Option<usize>,
    widen: bool,
    inspect_image: impl FnMut(u64, &mut Vec<Polynomial<Fp>>),
    stats: &mut ReconstructionStats,
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    let mut run = Reconstruction {
        input,
        primes: primes.into_iter(),
        groups: HashMap::new(),
        trace: TraceState::default(),
        check: None,
        spare: None,
        inspect_image,
        stats,
    };
    let mut wide = false;
    loop {
        let count = run.count();
        let ps: Vec<_> = run
            .primes
            .by_ref()
            .take(schedule(count, batch, wide))
            .collect();
        if ps.is_empty() {
            return Err(GroebnerError::ReconstructionFailed);
        }
        let start = std::time::Instant::now();
        run.round(&ps, count == 0)?;
        if count == 0 {
            wide = widen && start.elapsed() < WIDE_LEARN;
        }
        let Some((key, candidates)) = run.candidates()? else {
            continue;
        };
        for candidate in candidates {
            if run.validate(&key, &candidate)? {
                return Ok(candidate);
            }
        }
    }
}

// Leading monomials of a group, and the candidates reconstructed from it.
type Candidates = (Vec<Monomial>, Vec<Vec<Polynomial<BigRational>>>);

struct Reconstruction<'a, P, I> {
    input: &'a [Polynomial<BigRational>],
    primes: P,
    groups: HashMap<Vec<Monomial>, Group>,
    trace: TraceState,
    // The full check needs no candidate until the comparison, so it runs at a
    // reserved prime alongside the first learned image, which leaves most of the
    // pool idle.
    check: Option<Residues>,
    // The newest replayed image is held out of the accumulators as the next
    // candidate's agreement check, which saves replaying a fresh prime.
    spare: Option<Residues>,
    inspect_image: I,
    stats: &'a mut ReconstructionStats,
}

impl<P, I> Reconstruction<'_, P, I>
where
    P: Iterator<Item = u64>,
    I: FnMut(u64, &mut Vec<Polynomial<Fp>>),
{
    fn count(&self) -> usize {
        self.groups
            .values()
            .map(|g| g.primary.image_count())
            .max()
            .unwrap_or(0)
    }

    // Accumulate the images at `ps`, holding the newest back as the spare. The first
    // round also runs the full check.
    fn round(&mut self, ps: &[u64], first: bool) -> Result<(), GroebnerError> {
        let (input, trace) = (self.input, self.trace.current.as_ref());
        let reserved = first
            .then(|| {
                self.primes
                    .by_ref()
                    .find_map(|p| Some((p, map_input(input, p)?)))
            })
            .flatten();
        let (images, reserved) = par::join(
            || par::flat_map_chunks(ps, PRIMES, |ps| chunk_images(input, trace, ps)),
            || reserved.map(full_check).transpose(),
        );
        if let Some(reserved) = reserved? {
            self.check = Some(reserved);
        }
        if let Some((p, basis)) = self.spare.take() {
            self.add(p, &basis)?;
        }
        let last = (ps.len() > 1).then(|| ps.len() - 1);
        for (i, (&p, image)) in ps.iter().zip(images).enumerate() {
            let Some((mut basis, learned)) = image? else {
                continue;
            };
            self.stats.images += 1;
            (self.inspect_image)(p, &mut basis);
            self.observe(&basis, learned);
            if Some(i) == last {
                self.spare = Some((p, basis));
                continue;
            }
            self.add(p, &basis)?;
        }
        Ok(())
    }

    fn add(&mut self, p: u64, basis: &[Polynomial<Fp>]) -> Result<(), GroebnerError> {
        let group = self
            .groups
            .entry(heads(basis))
            .or_insert_with(|| Group::new(basis.len()));
        if group.add(p, basis)? {
            self.stats.recovery_restarts += 1;
        }
        Ok(())
    }

    fn observe(&mut self, basis: &[Polynomial<Fp>], learned: Option<F4Trace>) {
        if self.trace.observe(&heads(basis), learned) {
            self.stats.trace_replacements += 1;
        }
    }

    fn candidates(&mut self) -> Result<Option<Candidates>, GroebnerError> {
        let mut found = self.attempt();
        // Without the spare's prime reconstruction may fall just short. Failed
        // attempts are cheap, so retry with it rather than replay another batch.
        if found.as_ref().is_none_or(|(_, c)| c.is_empty())
            && let Some((p, basis)) = self.spare.take()
        {
            self.add(p, &basis)?;
            found = self.attempt();
        }
        Ok(found)
    }

    fn attempt(&self) -> Option<Candidates> {
        let (key, group) = self
            .groups
            .iter()
            .max_by_key(|(_, g)| g.primary.image_count())?;
        let candidates: Vec<_> = std::iter::once(&group.primary)
            .chain(group.recovery.iter())
            .filter_map(Accumulator::reconstruct)
            .map(|values| group.decode(values, &self.input[0]))
            .collect();
        Some((key.clone(), candidates))
    }

    // One cheap agreement image, followed by one FULL independent F4 run. Request them
    // individually: a candidate never triggers a whole batch.
    fn validate(
        &mut self,
        key: &[Monomial],
        candidate: &[Polynomial<BigRational>],
    ) -> Result<bool, GroebnerError> {
        Ok(self.agrees(false, key, candidate)? && self.agrees(true, key, candidate)?)
    }

    // Whether the candidate matches the image at a new prime, replayed or, with `full`,
    // computed without the trace. The image joins the accumulators unless it accepts.
    fn agrees(
        &mut self,
        full: bool,
        key: &[Monomial],
        candidate: &[Polynomial<BigRational>],
    ) -> Result<bool, GroebnerError> {
        let input = self.input;
        let reduce = |p| reduce_candidate(candidate, p);
        let reserved = if full {
            self.check.take()
        } else {
            self.spare.take_if(|(_, basis)| heads(basis) == key)
        };
        let reserved = reserved
            .and_then(|(p, basis)| Some((p, map_input(input, p)?, reduce(p)?, Some(basis))));
        let Some((p, mapped, reduced, known)) = reserved.or_else(|| {
            self.primes
                .by_ref()
                .find_map(|p| Some((p, map_input(input, p)?, reduce(p)?, None)))
        }) else {
            return Err(GroebnerError::ReconstructionFailed);
        };
        let (basis, learned) = match known {
            Some(basis) => (basis, None),
            // The untraced run is faster and shares nothing with the trace.
            None if full => (groebner_basis_f4_direct(mapped.clone(), true)?, None),
            None => replay_or_learn(self.trace.current.as_ref(), mapped.clone())?,
        };
        self.stats.images += 1;
        let agrees = basis == reduced;
        if agrees && full {
            return Ok(true);
        }
        if agrees {
            self.observe(&basis, learned);
        } else {
            // Learn a trace from a failed full check only to diagnose it.
            let learned = match learned {
                None if full => Some(learn(mapped.clone())?.1),
                learned => learned,
            };
            self.reject(key, mapped, &basis, learned)?;
        }
        self.add(p, &basis)?;
        Ok(agrees)
    }

    // Diagnose a trace only against a full result: a bad rational guess alone is not
    // evidence that its trace is broken. Recovery can then discard poisoned tails,
    // while the primary accumulator continues to retain every matching image.
    fn reject(
        &mut self,
        key: &[Monomial],
        mapped: Vec<Polynomial<Fp>>,
        basis: &[Polynomial<Fp>],
        learned: Option<F4Trace>,
    ) -> Result<(), GroebnerError> {
        self.stats.validation_failures += 1;
        if let Some(learned) = learned
            && self
                .trace
                .current
                .as_ref()
                .and_then(|t| t.replay(mapped))
                .as_deref()
                != Some(basis)
        {
            self.trace.current = Some(learned);
            self.trace.failed_leading = None;
            self.stats.trace_replacements += 1;
        }
        let group = self
            .groups
            .get_mut(key)
            .ok_or(GroebnerError::ReconstructionFailed)?;
        let width = group.width();
        group
            .recovery
            .get_or_insert_with(|| Accumulator::new(width));
        Ok(())
    }
}

// The image of mapped inputs by replaying the trace, or by learning a new one.
fn replay_or_learn(
    trace: Option<&F4Trace>,
    mapped: Vec<Polynomial<Fp>>,
) -> Result<Learned, GroebnerError> {
    if let Some(basis) = trace.and_then(|t| t.replay(mapped.clone())) {
        return Ok((basis, None));
    }
    learn(mapped).map(|(basis, trace)| (basis, Some(trace)))
}

// A full chunk of primes replays in one pass; the rest go one at a time.
fn chunk_images(
    input: &[Polynomial<BigRational>],
    trace: Option<&F4Trace>,
    ps: &[u64],
) -> Vec<Result<Option<Learned>, GroebnerError>> {
    let lanes = <[u64; PRIMES]>::try_from(ps)
        .ok()
        .zip(trace)
        .and_then(|(primes, t)| {
            let inputs = ps
                .iter()
                .map(|&p| map_input(input, p))
                .collect::<Option<_>>()?;
            t.replay_lanes(inputs, primes)
        });
    match lanes {
        Some(bases) => bases.into_iter().map(|b| Ok(Some((b, None)))).collect(),
        None => ps
            .iter()
            .map(|&p| {
                map_input(input, p)
                    .map(|fs| replay_or_learn(trace, fs))
                    .transpose()
            })
            .collect(),
    }
}

fn full_check((p, mapped): (u64, Vec<Polynomial<Fp>>)) -> Result<Residues, GroebnerError> {
    groebner_basis_f4_direct(mapped, true).map(|basis| (p, basis))
}

fn reduce_candidate(candidate: &[Polynomial<BigRational>], p: u64) -> Option<Vec<Polynomial<Fp>>> {
    par::map(candidate, |f| f.try_map(|c| Fp::from_rational(c, p)))
        .into_iter()
        .collect()
}

fn map_input(input: &[Polynomial<BigRational>], p: u64) -> Option<Vec<Polynomial<Fp>>> {
    input
        .iter()
        .map(|f| {
            // Primitive integer inputs have no bad denominators, but their heads can vanish.
            if f.lc().is_some_and(|c| crt::reduce(c, p) == Some(0)) {
                return None;
            }
            f.try_map(|c| Fp::from_rational(c, p))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{MonomialOrder, PolynomialRing};
    use num_traits::One;
    fn input(text: &str) -> Vec<Polynomial<BigRational>> {
        PolynomialRing::new(["x", "y"], MonomialOrder::Lex)
            .expect("ring")
            .parse_many(text)
            .expect("polynomials")
    }
    #[test]
    fn certification_rejects_strictly_larger_ideals() {
        let polys = input("x^2; y");
        for candidate in [input("1"), input("x; y")] {
            // Both old checks pass, despite these candidates generating a larger ideal.
            assert!(crate::is_groebner_basis(&candidate).unwrap());
            assert!(polys.iter().all(|p| p.reduce(&candidate).is_zero()));
            assert!(!certify_basis(&polys, &candidate).unwrap());
        }
    }

    #[test]
    fn certification_rejects_smaller_ideals_and_non_groebner_generators() {
        assert!(!certify_basis(&input("x; y"), &input("x^2; y")).unwrap());
        let polys = input("x^2 - y; x*y - 1");
        assert!(!crate::is_groebner_basis(&polys).unwrap());
        assert!(!certify_basis(&polys, &polys).unwrap());
    }

    #[test]
    fn exact_certification_accepts_reconstructed_bases() {
        for text in ["2/3*x^2 - 2/3*y; -3/5*x*y + 3/5", "x^2; y", "x; x - 1"] {
            let polys = input(text);
            let basis = groebner_basis_f4_rational(polys.clone(), true).unwrap();
            for batch in [1, 5, 40] {
                let options = RationalOptions { batch: Some(batch) };
                let wide = groebner_basis_f4_rational_with(polys.clone(), true, options);
                assert_eq!(wide.unwrap(), basis);
            }
            assert_eq!(basis, groebner_basis_f4_direct(polys, true).unwrap());
        }
    }

    #[test]
    fn missing_tail_support_is_zero_filled() {
        let polys = input("x + 6*y + 1");
        let primes = [2, 3, 5, 7].into_iter().chain(Primes::below(100000));
        let basis = reconstruct(&polys, primes, |_, _| {}).expect("reconstruction");
        assert_eq!(
            basis,
            groebner_basis_f4_direct(polys, true).expect("direct")
        );
    }
    #[test]
    fn leading_coefficient_primes_are_skipped() {
        let polys = input("6*x + y; y^2 - 1");
        let primes = [2, 3].into_iter().chain(Primes::below(100000));
        assert_eq!(
            reconstruct(&polys, primes, |_, _| {}).expect("reconstruction"),
            groebner_basis_f4_direct(polys, true).expect("direct")
        );
    }
    #[test]
    fn unlucky_leading_set_is_outvoted() {
        let polys = input("x + y; x - y");
        let primes = [2].into_iter().chain(Primes::below(100000));
        assert_eq!(
            reconstruct(&polys, primes, |_, _| {}).expect("reconstruction"),
            groebner_basis_f4_direct(polys, true).expect("direct")
        );
    }
    #[test]
    fn poisoned_tail_recovers_without_discarding_primary() {
        let polys = input("x + y + 1");
        let mut images = 0;
        let result = reconstruct(&polys, Primes::below(1 << 31), |_, basis| {
            images += 1;
            if images == 1 {
                basis[0].terms[1].1 = basis[0].terms[1].1 + Fp::one();
            }
        })
        .expect("restarts");
        assert!(
            images < 32,
            "a failed candidate should start recovery promptly"
        );
        assert_eq!(
            result,
            groebner_basis_f4_direct(polys, true).expect("direct")
        );
    }
    #[test]
    fn large_coefficients_keep_all_accumulated_images() {
        let mut polys = input("x + y");
        polys[0].terms[1].1 = BigRational::from_integer(num_bigint::BigInt::from(1u32) << 900);
        let mut stats = ReconstructionStats::default();
        let result = reconstruct_inner(
            &polys,
            Primes::below(1 << 31),
            None,
            false,
            |_, _| {},
            &mut stats,
        )
        .expect("reconstruction");
        assert_eq!(result, polys);
        assert!(
            stats.images <= 75,
            "healthy images must not be recomputed: {stats:?}"
        );
    }

    #[test]
    fn quick_first_image_widens_rounds() {
        let mut polys = input("x + y");
        polys[0].terms[1].1 = BigRational::from_integer(num_bigint::BigInt::from(1u32) << 900);
        let mut stats = ReconstructionStats::default();
        let result = reconstruct_inner(
            &polys,
            Primes::below(1 << 31),
            None,
            true,
            |_, _| {},
            &mut stats,
        )
        .expect("reconstruction");
        assert_eq!(result, polys);
        assert!(
            stats.images > PRIMES * par::threads(),
            "a quick first image should widen the replay rounds: {stats:?}"
        );
    }

    #[test]
    fn poisoned_large_coefficients_recover_with_growing_windows() {
        let mut polys = input("x + y");
        polys[0].terms[1].1 = BigRational::from_integer(num_bigint::BigInt::from(1u32) << 900);
        let mut stats = ReconstructionStats::default();
        let mut images = 0;
        let result = reconstruct_inner(
            &polys,
            Primes::below(1 << 31),
            None,
            false,
            |_, basis| {
                images += 1;
                if images == 1 {
                    basis[0].terms[1].1 = basis[0].terms[1].1 + Fp::one();
                }
            },
            &mut stats,
        )
        .expect("recovery");
        assert_eq!(result, polys);
        assert!(stats.recovery_restarts > 0, "{stats:?}");
        assert!(stats.images < 180, "{stats:?}");
    }

    #[test]
    fn independent_check_replaces_a_self_consistent_bad_trace() {
        let polys = input("x + y; x - y");
        let mut stats = ReconstructionStats::default();
        let primes = [2, 3, 5, 7].into_iter().chain(Primes::below(100000));
        let result = reconstruct_inner(&polys, primes, None, false, |_, _| {}, &mut stats)
            .expect("recovery");
        assert_eq!(
            result,
            groebner_basis_f4_direct(polys, true).expect("direct")
        );
        assert!(stats.validation_failures > 0);
        assert!(stats.trace_replacements > 0);
    }

    #[test]
    fn repeated_fallbacks_replace_a_trace_but_one_failure_does_not() {
        let polys = input("x + y; x - y");
        let mut trace = TraceState::default();
        for (p, expected) in [(3, false), (2, false), (5, false), (7, true)] {
            let (basis, learned) = learn(map_input(&polys, p).expect("map")).expect("learn");
            assert_eq!(trace.observe(&heads(&basis), Some(learned)), expected);
        }
        assert_eq!(
            trace
                .current
                .unwrap()
                .replay(map_input(&polys, 11).expect("map")),
            Some(
                groebner_basis_f4_direct(map_input(&polys, 11).expect("map"), true).expect("full")
            )
        );
    }

    #[test]
    fn exhausted_primes_report_failure() {
        assert!(reconstruct(&input("x + y"), [], |_, _| {}).is_err());
    }
}
