//! Multi-modular reconstruction of a reduced rational F4 basis.
use crate::f4::{F4Trace, groebner_basis_f4_direct, learn};
use crate::{Fp, GroebnerError, Monomial, Polynomial};
use num_rational::BigRational;
use polycore::{Primes, crt};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

type Layout = Vec<Vec<Monomial>>;
type Image = Option<(Vec<Polynomial<Fp>>, Option<F4Trace>)>;
#[derive(Default)]
struct Schema {
    layout: Layout,
    revision: usize,
}
// Layout revisions isolate CRT accumulators when a previously vanishing tail appears.
// Within each leading-monomial group all images use the union of observed supports.
#[derive(Clone, PartialEq, Eq, Hash)]
struct Key {
    leading: Vec<Monomial>,
    revision: usize,
}

/// Compute a rational basis by CRT and rational reconstruction, checked at a fresh prime.
/// With `certify`, also check Buchberger's criterion and reduction of every input over Q.
pub fn groebner_basis_f4_rational(
    polynomials: Vec<Polynomial<BigRational>>,
    certify: bool,
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    let input = crate::groebner::prepare_input(polynomials)?;
    let primitive: Vec<_> = input.iter().map(Polynomial::primitive).collect();
    // 23-bit primes allow fully deferred u64 accumulation for matrices with up to
    // 2^18 columns, trading more CRT images for much cheaper row elimination.
    let basis = reconstruct(&primitive, Primes::below(1 << 23), |_, _| {})?;
    if certify
        && (!crate::is_groebner_basis(&basis)? || input.iter().any(|p| !p.reduce(&basis).is_zero()))
    {
        return Err(GroebnerError::ReconstructionFailed);
    }
    Ok(basis)
}

fn reconstruct(
    input: &[Polynomial<BigRational>],
    primes: impl IntoIterator<Item = u64>,
    mut inspect_image: impl FnMut(u64, &mut Vec<Polynomial<Fp>>),
) -> Result<Vec<Polynomial<BigRational>>, GroebnerError> {
    let mut primes = primes.into_iter().peekable();
    let mut layouts: HashMap<Vec<Monomial>, Schema> = HashMap::new();
    // Bound each CRT epoch: even a bad tail with the right leading monomials cannot
    // poison an accumulator forever. Increasing the budget also admits large coefficients.
    // A power-of-two limit leaves one verification prime after the complete
    // 1 + 2 + ... + 64 reconstruction batches, rather than discarding their candidate.
    let mut budget = 32usize;
    let mut trace: Option<F4Trace> = None;
    loop {
        if primes.peek().is_none() {
            return Err(GroebnerError::ReconstructionFailed);
        }
        let mut failure = None;
        let result = crt::reconstruct_voted(
            |ps| {
                let image = |&p: &u64| -> Result<Image, GroebnerError> {
                    let mapped = map_input(input, p);
                    mapped
                        .map(|fs| {
                            if let Some(t) = &trace
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
                let images: Vec<_> = images
                    .into_iter()
                    .zip(ps)
                    .map(|(image, &p)| {
                        image.map(|image| {
                            image.map(|(mut basis, learned)| {
                                inspect_image(p, &mut basis);
                                if trace.is_none() {
                                    trace = learned;
                                }
                                basis
                            })
                        })
                    })
                    .collect();
                // Expand all schemas before encoding any images in this batch.
                for basis in images
                    .iter()
                    .filter_map(|r| r.as_ref().ok().and_then(Option::as_ref))
                {
                    let leading: Vec<_> = basis.iter().filter_map(|f| f.lm().cloned()).collect();
                    let schema = layouts.entry(leading).or_insert_with(|| Schema {
                        layout: vec![Vec::new(); basis.len()],
                        revision: 0,
                    });
                    let mut changed = false;
                    for (support, f) in schema.layout.iter_mut().zip(basis) {
                        if support.iter().eq(f.terms.iter().map(|t| &t.0)) {
                            continue;
                        }
                        let mut known: HashSet<_> = support.iter().cloned().collect();
                        for (m, _) in &f.terms {
                            if known.insert(m.clone()) {
                                support.push(m.clone());
                                changed = true;
                            }
                        }
                        support.sort_by(|a, b| f.order.compare(b, a));
                    }
                    if changed {
                        schema.revision += 1;
                    }
                }
                images
                    .into_iter()
                    .map(|image| {
                        let basis = match image {
                            Ok(basis) => basis?,
                            Err(e) => {
                                failure = Some(e);
                                return None;
                            }
                        };
                        let leading: Vec<_> =
                            basis.iter().filter_map(|f| f.lm().cloned()).collect();
                        let schema = layouts.get(&leading)?;
                        let mut values = Vec::new();
                        for (support, f) in schema.layout.iter().zip(&basis) {
                            let mut terms = f.terms.iter().peekable();
                            for m in support {
                                if terms.peek().is_some_and(|t| &t.0 == m) {
                                    values.push(terms.next()?.1.value());
                                } else {
                                    values.push(0);
                                }
                            }
                        }
                        Some((
                            Key {
                                leading,
                                revision: schema.revision,
                            },
                            values,
                        ))
                    })
                    .collect()
            },
            primes.by_ref().take(budget),
        );
        if let Some(e) = failure {
            return Err(e);
        }
        if let Some((key, values)) = result {
            let mut values = values.into_iter();
            let schema = layouts
                .remove(&key.leading)
                .ok_or(GroebnerError::ReconstructionFailed)?;
            if schema.revision != key.revision {
                return Err(GroebnerError::ReconstructionFailed);
            }
            let candidate: Vec<_> = schema
                .layout
                .into_iter()
                .map(|support| {
                    let terms = support.into_iter().zip(values.by_ref()).collect();
                    Polynomial::new(terms, input[0].nvars, input[0].order.clone())
                })
                .collect();
            // A fresh full run validates the candidate independently of the learned trace.
            for p in primes.by_ref() {
                let Some(mapped) = map_input(input, p) else {
                    continue;
                };
                let Some(reduced): Option<Vec<_>> = candidate
                    .iter()
                    .map(|f| f.try_map(|c| Fp::from_rational(c, p)))
                    .collect()
                else {
                    continue;
                };
                if groebner_basis_f4_direct(mapped, true)? == reduced {
                    return Ok(candidate);
                }
                break;
            }
        }
        trace = None;
        budget = budget.saturating_mul(2);
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
    use crate::{MonomialOrder, PolynomialRing};
    use num_traits::One;
    fn input(text: &str) -> Vec<Polynomial<BigRational>> {
        PolynomialRing::new(["x", "y"], MonomialOrder::Lex)
            .expect("ring")
            .parse_many(text)
            .expect("polynomials")
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
    fn poisoned_tail_restarts_without_changing_leading_set() {
        let polys = input("x + y + 1");
        let mut images = 0;
        let result = reconstruct(&polys, Primes::below(1 << 31), |_, basis| {
            images += 1;
            if images == 1 {
                basis[0].terms[1].1 = basis[0].terms[1].1 + Fp::one();
            }
        })
        .expect("restarts");
        assert!(images > 32, "poisoned accumulator must be discarded");
        assert_eq!(
            result,
            groebner_basis_f4_direct(polys, true).expect("direct")
        );
    }
    #[test]
    fn last_batch_candidate_gets_a_verification_prime_before_retry() {
        let mut polys = input("x + y");
        polys[0].terms[1].1 = BigRational::from_integer(num_bigint::BigInt::from(1u32) << 450);
        let mut images = 0;
        let result = reconstruct(&polys, Primes::below(1 << 31), |_, _| images += 1)
            .expect("reconstruction");
        assert_eq!(
            images, 32,
            "the candidate from 31 images should be checked before restarting"
        );
        assert_eq!(result, polys);
    }

    #[test]
    fn exhausted_primes_report_failure() {
        assert!(reconstruct(&input("x + y"), [], |_, _| {}).is_err());
    }
}
