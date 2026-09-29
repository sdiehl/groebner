//! Multi-modular reconstruction of a reduced rational F4 basis.
use crate::f4::{F4Trace, PRIMES, groebner_basis_f4_direct, learn};
use crate::{Fp, GroebnerError, Monomial, Polynomial, lehmer};
use num_rational::BigRational;
use polycore::{Primes, crt, modp};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::HashMap;

type Layout = Vec<Vec<Monomial>>;
type Image = Option<(Vec<Polynomial<Fp>>, Option<F4Trace>)>;
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
    /// or have large matrices. Small systems with tall coefficients finish sooner
    /// with a round wide enough to fill the thread pool, such as four per thread.
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
        |_, _| {},
        #[cfg(test)]
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

#[cfg(test)]
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

/// Coefficients in Garner's mixed radix form, one row of digits per prime:
/// `x = d[0] + p[0] * (d[1] + p[1] * (d[2] + ...))`. Adding a prime only takes
/// machine word arithmetic; big integers appear when reconstructing.
#[derive(Default)]
struct Accumulator {
    primes: Vec<u64>,
    digits: Vec<Vec<u32>>,
}

impl Accumulator {
    fn count(&self) -> usize {
        self.primes.len()
    }

    fn add(&mut self, p: u64, values: &[u64]) {
        assert!(p < 1 << 32, "mixed radix digits are u32");
        // `x mod p` is the dot product of the digits with these weights.
        let mut weights = Vec::with_capacity(self.primes.len());
        let mut m = 1;
        for &q in &self.primes {
            weights.push(m as u32);
            m = modp::mul(m, q % p, p);
        }
        let minv = modp::inv(m, p);
        let max = self.primes.iter().copied().fold(p, u64::max);
        let lazy = u128::from(max).pow(2) * (self.primes.len() as u128 + 1) < 1 << 64;
        let digits = &self.digits;
        let fill = |(c, out): (usize, &mut [u32])| {
            let span = c * CHUNK..c * CHUNK + out.len();
            let mut acc = vec![0u64; out.len()];
            for (row, &w) in digits.iter().zip(&weights) {
                let row = &row[span.clone()];
                if lazy {
                    for (a, &t) in acc.iter_mut().zip(row) {
                        *a += u64::from(t) * u64::from(w);
                    }
                } else {
                    for (a, &t) in acc.iter_mut().zip(row) {
                        *a = modp::add(*a, modp::mul(u64::from(t), u64::from(w), p), p);
                    }
                }
            }
            for ((o, &v), a) in out.iter_mut().zip(&values[span]).zip(acc) {
                *o = modp::mul(modp::sub(v, a % p, p), minv, p) as u32;
            }
        };
        const CHUNK: usize = 4096;
        let mut next = vec![0u32; values.len()];
        #[cfg(feature = "parallel")]
        next.par_chunks_mut(CHUNK).enumerate().for_each(fill);
        #[cfg(not(feature = "parallel"))]
        next.chunks_mut(CHUNK).enumerate().for_each(fill);
        self.digits.push(next);
        self.primes.push(p);
    }

    fn remap(&mut self, positions: &[Option<usize>]) {
        for row in &mut self.digits {
            *row = positions.iter().map(|i| i.map_or(0, |i| row[i])).collect();
        }
    }

    fn reconstruct(&self) -> Option<Vec<BigRational>> {
        if self.primes.is_empty() {
            return None;
        }
        // Consecutive primes whose product fits a word form one big radix digit.
        let mut radices: Vec<(u64, std::ops::Range<usize>)> = Vec::new();
        for (j, &p) in self.primes.iter().enumerate() {
            match radices.last_mut() {
                Some((r, js)) if r.checked_mul(p).is_some() => {
                    *r *= p;
                    js.end = j + 1;
                }
                _ => radices.push((p, j..j + 1)),
            }
        }
        let value = |i: usize| -> num_bigint::BigInt {
            let mut x = num_bigint::BigUint::default();
            for (r, js) in radices.iter().rev() {
                let d = js
                    .clone()
                    .rev()
                    .fold(0, |d, j| d * self.primes[j] + u64::from(self.digits[j][i]));
                x *= *r;
                x += d;
            }
            x.into()
        };
        let modulus: num_bigint::BigInt = radices
            .iter()
            .map(|(r, _)| *r)
            .product::<num_bigint::BigUint>()
            .into();
        let n = self.digits[0].len();
        // Probe spread-out coefficients before reconstructing a potentially huge basis.
        // Checking only the final coefficient often checks a trivial zero or one.
        // These probes are a heuristic: an unprobed large coefficient can still
        // make the full reconstruction fail again on the next batch.
        if modulus <= num_bigint::BigInt::from(1u32) {
            return None;
        }
        let bound = (&modulus / 2u32).sqrt();
        let probes = 64.min(n);
        for i in 0..probes {
            let index = i * (n - 1) / probes.saturating_sub(1).max(1);
            lehmer::wang(&value(index), &modulus, &bound)?;
        }
        // A failure anywhere fails the attempt, so the other chunks stop early.
        let failed = std::sync::atomic::AtomicBool::new(false);
        let run = |c: usize| -> Option<Vec<BigRational>> {
            // Coefficients of one polynomial mostly share denominators. Once `d` holds
            // theirs, `x * d` is a small integer and needs no half extended gcd. Both
            // parts are within the Wang bound, so this is the fraction Wang would find.
            let mut d = num_bigint::BigInt::from(1u32);
            (c * 256..n.min(c * 256 + 256))
                .map(|i| {
                    if failed.load(std::sync::atomic::Ordering::Relaxed) {
                        return None;
                    }
                    let x = value(i);
                    let y = crt::symmetric(&(&x * &d), &modulus);
                    if y.magnitude() <= bound.magnitude() {
                        return Some(lehmer::fraction(y, d.magnitude()));
                    }
                    let Some(c) = lehmer::wang(&x, &modulus, &bound) else {
                        failed.store(true, std::sync::atomic::Ordering::Relaxed);
                        return None;
                    };
                    let g = lehmer::gcd(d.magnitude(), c.denom().magnitude());
                    let lcm = &d / num_bigint::BigInt::from(g) * c.denom();
                    d = if lcm <= bound { lcm } else { c.denom().clone() };
                    Some(c)
                })
                .collect()
        };
        #[cfg(feature = "parallel")]
        return (0..n.div_ceil(256))
            .into_par_iter()
            .map(run)
            .collect::<Option<Vec<_>>>()
            .map(|v| v.concat());
        #[cfg(not(feature = "parallel"))]
        (0..n.div_ceil(256))
            .map(run)
            .collect::<Option<Vec<_>>>()
            .map(|v| v.concat())
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
            primary: Accumulator::default(),
            recovery: None,
            recovery_limit: 32,
        }
    }

    fn add(&mut self, p: u64, basis: &[Polynomial<Fp>]) -> bool {
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
            .is_some_and(|a| a.count() >= self.recovery_limit);
        if restart {
            self.recovery_limit = self.recovery_limit.saturating_mul(2);
            self.recovery = Some(Accumulator::default());
        } else if self.recovery.is_none() && self.primary.count() >= 32 {
            self.recovery = Some(Accumulator::default());
        }
        self.primary.add(p, &values);
        if let Some(recovery) = &mut self.recovery {
            recovery.add(p, &values);
        }
        restart
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

fn add_image(
    groups: &mut HashMap<Vec<Monomial>, Group>,
    p: u64,
    basis: &[Polynomial<Fp>],
    #[cfg(test)] stats: &mut ReconstructionStats,
) {
    let group = groups
        .entry(heads(basis))
        .or_insert_with(|| Group::new(basis.len()));
    if group.add(p, basis) {
        #[cfg(test)]
        {
            stats.recovery_restarts += 1;
        }
    }
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
        inspect_image,
        #[cfg(test)]
        &mut ReconstructionStats::default(),
    )
}

// After the learned image, replay whole lane chunks and grow by about 12.5%, rather
// than doubling the work near completion, unless the caller fixed the width.
fn schedule(count: usize, batch: Option<usize>) -> usize {
    match (count, batch) {
        (0, _) => 1,
        (_, Some(batch)) => batch.max(1),
        _ => count
            .div_ceil(8)
            .max(2 * PRIMES)
            .next_multiple_of(PRIMES)
            .min(32),
    }
}

fn reconstruct_inner(
    input: &[Polynomial<BigRational>],
    primes: impl IntoIterator<Item = u64>,
    batch: Option<usize>,
    mut inspect_image: impl FnMut(u64, &mut Vec<Polynomial<Fp>>),
    #[cfg(test)] stats: &mut ReconstructionStats,
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    let mut primes = primes.into_iter();
    let mut groups: HashMap<Vec<Monomial>, Group> = HashMap::new();
    let mut trace = TraceState::default();
    // The full check needs no candidate until the comparison, so it runs at a
    // reserved prime alongside the first learned image, which leaves most of the
    // pool idle.
    let mut check = None;
    // The newest replayed image is held out of the accumulators as the next
    // candidate's agreement check, which saves replaying a fresh prime.
    let mut spare: Option<Residues> = None;
    loop {
        let count = groups
            .values()
            .map(|g| g.primary.count())
            .max()
            .unwrap_or(0);
        let ps: Vec<_> = primes.by_ref().take(schedule(count, batch)).collect();
        if ps.is_empty() {
            return Err(GroebnerError::ReconstructionFailed);
        }

        let chunk = |ps: &[u64]| chunk_images(input, trace.current.as_ref(), ps);
        let reserved = (count == 0).then(|| {
            primes
                .by_ref()
                .find_map(|p| Some((p, map_input(input, p)?)))
        });
        #[cfg(feature = "parallel")]
        let (images, reserved) = rayon::join(
            || {
                ps.par_chunks(PRIMES)
                    .flat_map_iter(chunk)
                    .collect::<Vec<_>>()
            },
            || reserved.flatten().map(full_check).transpose(),
        );
        #[cfg(not(feature = "parallel"))]
        let (images, reserved) = (
            ps.chunks(PRIMES).flat_map(chunk).collect::<Vec<_>>(),
            reserved.flatten().map(full_check).transpose(),
        );
        if let Some(reserved) = reserved? {
            check = Some(reserved);
        }

        if let Some((p, basis)) = spare.take() {
            add_image(
                &mut groups,
                p,
                &basis,
                #[cfg(test)]
                stats,
            );
        }
        let last = (ps.len() > 1).then(|| ps.len() - 1);
        for (i, (p, image)) in ps.into_iter().zip(images).enumerate() {
            let Some((mut basis, learned)) = image? else {
                continue;
            };
            #[cfg(test)]
            {
                stats.images += 1;
            }
            inspect_image(p, &mut basis);
            if trace.observe(&heads(&basis), learned) {
                #[cfg(test)]
                {
                    stats.trace_replacements += 1;
                }
            }
            if Some(i) == last {
                spare = Some((p, basis));
                continue;
            }
            add_image(
                &mut groups,
                p,
                &basis,
                #[cfg(test)]
                stats,
            );
        }
        let attempt = |groups: &HashMap<Vec<Monomial>, Group>| {
            let (key, group) = groups.iter().max_by_key(|(_, g)| g.primary.count())?;
            let candidates: Vec<_> = std::iter::once(&group.primary)
                .chain(group.recovery.iter())
                .filter_map(Accumulator::reconstruct)
                .map(|values| group.decode(values, &input[0]))
                .collect();
            Some((key.clone(), candidates))
        };
        let mut found = attempt(&groups);
        // Without the spare's prime reconstruction may fall just short. Failed
        // attempts are cheap, so retry with it rather than replay another batch.
        if found.as_ref().is_none_or(|(_, c)| c.is_empty())
            && let Some((p, basis)) = spare.take()
        {
            add_image(
                &mut groups,
                p,
                &basis,
                #[cfg(test)]
                stats,
            );
            found = attempt(&groups);
        }
        let Some((key, candidates)) = found else {
            continue;
        };

        for candidate in candidates {
            // One cheap agreement image, followed by one FULL independent F4 run.
            // Request them individually: a candidate never triggers a whole batch.
            let reduce = |p| reduce_candidate(&candidate, p);
            for stage in 0..2 {
                let reserved = if stage == 0 {
                    spare.take_if(|(_, basis)| heads(basis) == key)
                } else {
                    check.take()
                };
                let reserved = reserved.and_then(|(p, basis)| {
                    Some((p, map_input(input, p)?, reduce(p)?, Some(basis)))
                });
                let Some((p, mapped, reduced, full)) = reserved.or_else(|| {
                    primes
                        .by_ref()
                        .find_map(|p| Some((p, map_input(input, p)?, reduce(p)?, None)))
                }) else {
                    return Err(GroebnerError::ReconstructionFailed);
                };

                let (basis, learned) = if let Some(basis) = full {
                    (basis, None)
                } else if stage == 0 {
                    if let Some(basis) = trace
                        .current
                        .as_ref()
                        .and_then(|t| t.replay(mapped.clone()))
                    {
                        (basis, None)
                    } else {
                        let (basis, learned) = learn(mapped.clone())?;
                        (basis, Some(learned))
                    }
                } else {
                    // The untraced run is faster and shares nothing with the trace.
                    (groebner_basis_f4_direct(mapped.clone(), true)?, None)
                };
                #[cfg(test)]
                {
                    stats.images += 1;
                }
                let agrees = basis == reduced;
                if !agrees {
                    // Learn a trace from a failed full check only to diagnose it.
                    let learned = match learned {
                        None if stage == 1 => Some(learn(mapped.clone())?.1),
                        learned => learned,
                    };
                    #[cfg(test)]
                    {
                        stats.validation_failures += 1;
                    }
                    // Diagnose a trace only against a full result. A bad rational
                    // guess alone is not evidence that its trace is broken.
                    if let Some(learned) = learned
                        && trace
                            .current
                            .as_ref()
                            .and_then(|t| t.replay(mapped))
                            .as_ref()
                            != Some(&basis)
                    {
                        trace.current = Some(learned);
                        trace.failed_leading = None;
                        #[cfg(test)]
                        {
                            stats.trace_replacements += 1;
                        }
                    }
                    // Recovery can discard poisoned tails, while the primary
                    // accumulator continues to retain every matching image.
                    let group = groups
                        .get_mut(&key)
                        .ok_or(GroebnerError::ReconstructionFailed)?;
                    if group.recovery.is_none() {
                        group.recovery = Some(Accumulator::default());
                    }
                } else if stage == 0 && trace.observe(&heads(&basis), learned) {
                    #[cfg(test)]
                    {
                        stats.trace_replacements += 1;
                    }
                }

                if agrees && stage == 1 {
                    return Ok(candidate);
                }
                add_image(
                    &mut groups,
                    p,
                    &basis,
                    #[cfg(test)]
                    stats,
                );
                if !agrees {
                    break;
                }
            }
        }
    }
}

// A full chunk of primes replays in one pass; the rest go one at a time.
fn chunk_images(
    input: &[Polynomial<BigRational>],
    trace: Option<&F4Trace>,
    ps: &[u64],
) -> Vec<Result<Image, GroebnerError>> {
    let image = |&p: &u64| -> Result<Image, GroebnerError> {
        map_input(input, p)
            .map(|fs| {
                if let Some(t) = trace
                    && let Some(basis) = t.replay(fs.clone())
                {
                    return Ok((basis, None));
                }
                learn(fs).map(|(basis, trace)| (basis, Some(trace)))
            })
            .transpose()
    };
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
        None => ps.iter().map(image).collect(),
    }
}

fn full_check((p, mapped): (u64, Vec<Polynomial<Fp>>)) -> Result<Residues, GroebnerError> {
    groebner_basis_f4_direct(mapped, true).map(|basis| (p, basis))
}

fn reduce_candidate(candidate: &[Polynomial<BigRational>], p: u64) -> Option<Vec<Polynomial<Fp>>> {
    let reduce = |f: &Polynomial<BigRational>| f.try_map(|c| Fp::from_rational(c, p));
    #[cfg(feature = "parallel")]
    return candidate.par_iter().map(reduce).collect();
    #[cfg(not(feature = "parallel"))]
    candidate.iter().map(reduce).collect()
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
        let result = reconstruct_inner(&polys, Primes::below(1 << 31), None, |_, _| {}, &mut stats)
            .expect("reconstruction");
        assert_eq!(result, polys);
        assert!(
            stats.images <= 75,
            "healthy images must not be recomputed: {stats:?}"
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
        let result =
            reconstruct_inner(&polys, primes, None, |_, _| {}, &mut stats).expect("recovery");
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
