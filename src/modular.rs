//! Multi-modular reconstruction of a reduced rational F4 basis.
use crate::f4::{F4Trace, learn};
use crate::{Fp, GroebnerError, Monomial, Polynomial};
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::Signed;
use polycore::{Primes, crt, modp};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::HashMap;

type Layout = Vec<Vec<Monomial>>;
type Image = Option<(Vec<Polynomial<Fp>>, Option<F4Trace>)>;
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
    let input = crate::groebner::prepare_input(polynomials)?;
    let primitive: Vec<_> = input.iter().map(Polynomial::primitive).collect();
    // 23-bit primes allow fully deferred u64 accumulation for matrices with up to
    // 2^18 columns, trading more CRT images for much cheaper row elimination.
    let basis = reconstruct(&primitive, Primes::below(1 << 23), |_, _| {})?;
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

struct Accumulator {
    residues: Vec<num_bigint::BigInt>,
    modulus: num_bigint::BigInt,
    count: usize,
}
impl Default for Accumulator {
    fn default() -> Self {
        Self {
            residues: Vec::new(),
            modulus: 1u32.into(),
            count: 0,
        }
    }
}
// `x mod p` for `x >= 0`, by Horner over the limbs without allocating.
fn residue(x: &num_bigint::BigInt, p: u64) -> u64 {
    let p = u128::from(p);
    x.magnitude()
        .iter_u64_digits()
        .rev()
        .fold(0, |r, d| ((u128::from(r) << 64 | u128::from(d)) % p) as u64)
}

impl Accumulator {
    fn add(&mut self, p: u64, values: &[u64]) {
        if self.count == 0 {
            self.residues = values.iter().copied().map(Into::into).collect();
            self.modulus = p.into();
        } else {
            let m = &self.modulus;
            let minv = modp::inv(residue(m, p), p);
            let step = |(x, &v): (&mut num_bigint::BigInt, &u64)| {
                let t = modp::mul(modp::sub(v, residue(x, p), p), minv, p);
                if t != 0 {
                    *x += m * t;
                }
            };
            #[cfg(feature = "parallel")]
            if self.residues.len() >= 1024 {
                self.residues.par_iter_mut().zip(values).for_each(step);
            } else {
                self.residues.iter_mut().zip(values).for_each(step);
            }
            #[cfg(not(feature = "parallel"))]
            self.residues.iter_mut().zip(values).for_each(step);
            self.modulus *= p;
        }
        self.count += 1;
    }

    fn remap(&mut self, positions: &[Option<usize>]) {
        if self.count != 0 {
            self.residues = positions
                .iter()
                .map(|i| i.map_or_else(|| 0u32.into(), |i| self.residues[i].clone()))
                .collect();
        }
    }

    fn reconstruct(&self) -> Option<Vec<BigRational>> {
        if self.count == 0 {
            return None;
        }
        // Probe spread-out coefficients before reconstructing a potentially huge basis.
        // Checking only the final coefficient often checks a trivial zero or one.
        // These probes are a heuristic: an unprobed large coefficient can still
        // make the full reconstruction fail again on the next batch.
        let context = crt::WangContext::new(&self.modulus)?;
        let probes = 8.min(self.residues.len());
        for i in 0..probes {
            let index = i * (self.residues.len() - 1) / probes.saturating_sub(1).max(1);
            context.reconstruct(&self.residues[index])?;
        }
        let bound = (&self.modulus / 2u32).sqrt();
        let run = |xs: &[num_bigint::BigInt]| -> Option<Vec<BigRational>> {
            // Coefficients of one polynomial mostly share denominators. Once `d` holds
            // theirs, `x * d` is a small integer and needs no half extended gcd. Both
            // parts are within the Wang bound, so this is the fraction Wang would find.
            let mut d = num_bigint::BigInt::from(1u32);
            xs.iter()
                .map(|x| {
                    let y = crt::symmetric(&(x * &d), &self.modulus);
                    if y.abs() <= bound {
                        return Some(BigRational::new(y, d.clone()));
                    }
                    let c = context.reconstruct(x)?;
                    let lcm = d.lcm(c.denom());
                    d = if lcm <= bound { lcm } else { c.denom().clone() };
                    Some(c)
                })
                .collect()
        };
        #[cfg(feature = "parallel")]
        return self
            .residues
            .par_chunks(256)
            .map(run)
            .collect::<Option<Vec<_>>>()
            .map(|v| v.concat());
        #[cfg(not(feature = "parallel"))]
        run(&self.residues)
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
            .is_some_and(|a| a.count >= self.recovery_limit);
        if restart {
            self.recovery_limit = self.recovery_limit.saturating_mul(2);
            self.recovery = Some(Accumulator::default());
        } else if self.recovery.is_none() && self.primary.count >= 32 {
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

fn reconstruct(
    input: &[Polynomial<BigRational>],
    primes: impl IntoIterator<Item = u64>,
    inspect_image: impl FnMut(u64, &mut Vec<Polynomial<Fp>>),
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    reconstruct_inner(
        input,
        primes,
        inspect_image,
        #[cfg(test)]
        &mut ReconstructionStats::default(),
    )
}

fn reconstruct_inner(
    input: &[Polynomial<BigRational>],
    primes: impl IntoIterator<Item = u64>,
    mut inspect_image: impl FnMut(u64, &mut Vec<Polynomial<Fp>>),
    #[cfg(test)] stats: &mut ReconstructionStats,
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    let mut primes = primes.into_iter();
    let mut groups: HashMap<Vec<Monomial>, Group> = HashMap::new();
    let mut trace = TraceState::default();
    loop {
        let count = groups.values().map(|g| g.primary.count).max().unwrap_or(0);
        // Grow by about 12.5%, rather than doubling the work near completion. Once a
        // pool's worth of primes did not suffice, take at least one prime per worker:
        // a few replays alone leave the pool idle.
        #[cfg(feature = "parallel")]
        let width = rayon::current_num_threads();
        #[cfg(not(feature = "parallel"))]
        let width = 2;
        let batch = if count == 0 {
            1
        } else {
            let step = count.div_ceil(8);
            if count >= width {
                step.max(width)
            } else {
                step
            }
            .clamp(2, 32)
        };
        let ps: Vec<_> = primes.by_ref().take(batch).collect();
        if ps.is_empty() {
            return Err(GroebnerError::ReconstructionFailed);
        }

        let image = |&p: &u64| -> Result<Image, GroebnerError> {
            map_input(input, p)
                .map(|fs| {
                    if let Some(t) = &trace.current
                        && let Some(basis) = t.replay(fs.clone())
                    {
                        return Ok((basis, None));
                    }
                    learn(fs).map(|(basis, trace)| (basis, Some(trace)))
                })
                .transpose()
        };
        #[cfg(feature = "parallel")]
        let images: Vec<_> = ps.par_iter().map(image).collect();
        #[cfg(not(feature = "parallel"))]
        let images: Vec<_> = ps.iter().map(image).collect();

        for (p, image) in ps.into_iter().zip(images) {
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
            add_image(
                &mut groups,
                p,
                &basis,
                #[cfg(test)]
                stats,
            );
        }
        let Some((key, group)) = groups.iter().max_by_key(|(_, g)| g.primary.count) else {
            continue;
        };
        let key = key.clone();

        let candidates: Vec<_> = std::iter::once(&group.primary)
            .chain(group.recovery.iter())
            .filter_map(Accumulator::reconstruct)
            .map(|values| group.decode(values, &input[0]))
            .collect();

        for candidate in candidates {
            // One cheap agreement image, followed by one FULL independent F4 run.
            // Request them individually: a candidate never triggers a whole batch.
            for stage in 0..2 {
                let Some((p, mapped, reduced)) = primes.by_ref().find_map(|p| {
                    let mapped = map_input(input, p)?;
                    let reduced: Option<Vec<_>> = candidate
                        .iter()
                        .map(|f| f.try_map(|c| Fp::from_rational(c, p)))
                        .collect();
                    Some((p, mapped, reduced?))
                }) else {
                    return Err(GroebnerError::ReconstructionFailed);
                };

                let (basis, learned) = if stage == 0 {
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
                    let (basis, learned) = learn(mapped.clone())?;
                    (basis, Some(learned))
                };
                #[cfg(test)]
                {
                    stats.images += 1;
                }
                let agrees = basis == reduced;
                if !agrees {
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
    use crate::f4::groebner_basis_f4_direct;
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
        let result = reconstruct_inner(&polys, Primes::below(1 << 31), |_, _| {}, &mut stats)
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
        let result = reconstruct_inner(&polys, primes, |_, _| {}, &mut stats).expect("recovery");
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
