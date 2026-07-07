//! Buchberger-style Groebner basis algorithms.
//!
//! An implementation of Buchberger's algorithm and related utilities for computing Groebner bases
//! of polynomial ideals. It also provides functions for checking if a set of polynomials forms a Groebner basis.
//!
//! # Algorithms
//! - Buchberger's algorithm (with monic and minimal basis options)
//! - Basis minimization
//! - S-polynomial computation and reduction
//! - Canonicalization of polynomials
//!
//! # Example
//! ```
//! use groebner::{groebner_basis, groebner_basis_incremental, MonomialOrder, PolynomialRing};
//! use num_rational::BigRational;
//!
//! let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)?;
//! let f1 = ring.parse("x^2 - y")?;
//! let f2 = ring.parse("x*y - 1")?;
//! let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true)?;
//! assert!(!basis.is_empty());
//!
//! let f3 = ring.parse("y^2 - x")?;
//! let updated = groebner_basis_incremental(basis, vec![f3], ring.order(), true)?;
//! assert!(!updated.is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::field::Field;
use crate::grebauer_moller;
use crate::monomial::Monomial;
use crate::polynomial::Polynomial;
use crate::sugar::{select_next_by_sugar, SugaredPolynomial};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fmt;

#[derive(Debug)]
pub enum GroebnerError {
    NoLeadingMonomial(usize),
    EmptyInput,
    Polynomial(crate::polynomial::PolynomialError),
    InvalidPrimeField { expected: u32, actual: u32 },
}

impl From<crate::polynomial::PolynomialError> for GroebnerError {
    fn from(e: crate::polynomial::PolynomialError) -> Self {
        GroebnerError::Polynomial(e)
    }
}

impl fmt::Display for GroebnerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GroebnerError::NoLeadingMonomial(idx) => {
                write!(f, "Polynomial at index {idx} has no leading monomial")
            }
            GroebnerError::EmptyInput => write!(f, "Input polynomial list is empty"),
            GroebnerError::Polynomial(e) => write!(f, "Polynomial error: {e}"),
            GroebnerError::InvalidPrimeField { expected, actual } => write!(
                f,
                "runtime prime {actual} does not match coefficient field modulus {expected}"
            ),
        }
    }
}

impl std::error::Error for GroebnerError {}

#[derive(Clone, Debug)]
pub struct CriticalPair {
    pub i: usize,
    pub j: usize,
    pub lcm: Monomial,
    pub degree: u32,
}

impl CriticalPair {
    fn new<F: Field>(
        i: usize,
        j: usize,
        poly_i: &Polynomial<F>,
        poly_j: &Polynomial<F>,
    ) -> Result<Self, GroebnerError> {
        let lm_i = poly_i
            .leading_monomial()
            .ok_or(GroebnerError::NoLeadingMonomial(i))?;
        let lm_j = poly_j
            .leading_monomial()
            .ok_or(GroebnerError::NoLeadingMonomial(j))?;
        let lcm = lm_i.lcm(lm_j);
        let degree = lcm.degree();
        Ok(Self { i, j, lcm, degree })
    }
}

impl PartialEq for CriticalPair {
    fn eq(&self, other: &Self) -> bool {
        self.degree == other.degree && self.lcm == other.lcm
    }
}
impl Eq for CriticalPair {}
impl PartialOrd for CriticalPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for CriticalPair {
    fn cmp(&self, other: &Self) -> Ordering {
        self.degree.cmp(&other.degree).reverse()
    }
}

/// Enum for S-polynomial selection strategy
pub enum SelectionStrategy {
    Degree,        // Default: by degree
    Sugar,         // Use sugar strategy
    GebauerMoller, // Use Gebauer–Möller criteria
}

/// Compute Groebner basis
pub fn groebner_basis<F: Field>(
    polynomials: Vec<Polynomial<F>>,
    order: crate::monomial::MonomialOrder,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    groebner_basis_with_strategy(polynomials, order, canonicalize, &SelectionStrategy::Degree)
}

/// Compute a Groebner basis using parallel batches of equal-degree critical pairs.
///
/// This uses Buchberger's algorithm with degree-based pair selection, but each selected
/// minimal-degree batch is reduced against a fixed snapshot of the current basis in parallel.
/// The final basis is minimized and canonicalized in the same shape as [`groebner_basis`].
///
/// This API is available with the default `parallel` feature and requires coefficient fields to
/// be [`Send`] and [`Sync`] so reductions can run on Rayon worker threads.
#[cfg(feature = "parallel")]
pub fn groebner_basis_parallel<F: Field + Send + Sync>(
    polynomials: Vec<Polynomial<F>>,
    order: crate::monomial::MonomialOrder,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    if polynomials.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }
    let mut basis: Vec<Polynomial<F>> = polynomials
        .into_iter()
        .filter(|p| !p.is_zero())
        .map(|p| p.make_monic())
        .collect();
    if basis.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }

    let mut pairs = all_critical_pairs(&basis)?;
    while !pairs.is_empty() {
        let selected = select_min_degree_critical_pairs(&mut pairs);
        let snapshot = basis.clone();
        let mut new_polynomials = selected
            .par_iter()
            .map(|pair| reduce_pair_against_basis(pair, &snapshot))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<_>>();

        if new_polynomials.is_empty() {
            continue;
        }

        new_polynomials.sort_by(|a, b| match (a.leading_monomial(), b.leading_monomial()) {
            (Some(ma), Some(mb)) => mb.compare(ma, order),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => std::cmp::Ordering::Equal,
        });
        new_polynomials.dedup_by(|a, b| a.leading_monomial() == b.leading_monomial());

        for polynomial in new_polynomials {
            let reduced = polynomial.reduce(&basis).map_err(GroebnerError::from)?;
            if reduced.is_zero() {
                continue;
            }
            let monic_reduced = reduced.make_monic();
            let new_index = basis.len();
            for (i, existing) in basis.iter().enumerate() {
                pairs.push(CriticalPair::new(i, new_index, existing, &monic_reduced)?);
            }
            basis.push(monic_reduced);
        }
    }

    finish_basis(basis, canonicalize)
}

/// Update a Groebner basis with additional generators.
///
/// The existing basis is reused as reduction context for the new generators before the combined
/// system is completed with Buchberger. This keeps the API conservative and correctness-oriented:
/// callers can incrementally add generators without having to rebuild the input list themselves.
pub fn groebner_basis_incremental<F: Field>(
    existing_basis: Vec<Polynomial<F>>,
    new_polynomials: Vec<Polynomial<F>>,
    order: crate::monomial::MonomialOrder,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    if existing_basis.is_empty() && new_polynomials.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }

    let mut combined: Vec<_> = existing_basis
        .into_iter()
        .filter(|polynomial| !polynomial.is_zero())
        .collect();

    for polynomial in new_polynomials
        .into_iter()
        .filter(|polynomial| !polynomial.is_zero())
    {
        let reduced = if combined.is_empty() {
            polynomial
        } else {
            polynomial.reduce(&combined).map_err(GroebnerError::from)?
        };
        if !reduced.is_zero() {
            combined.push(reduced);
        }
    }

    groebner_basis_with_strategy(combined, order, canonicalize, &SelectionStrategy::Degree)
}

/// Compute Groebner basis with a selectable S-polynomial selection strategy
#[allow(clippy::needless_range_loop)]
pub fn groebner_basis_with_strategy<F: Field>(
    polynomials: Vec<Polynomial<F>>,
    _order: crate::monomial::MonomialOrder,
    canonicalize: bool,
    strategy: &SelectionStrategy,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    if polynomials.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }
    let _nvars = polynomials[0].nvars;
    let mut basis: Vec<Polynomial<F>> = polynomials
        .into_iter()
        .filter(|p| !p.is_zero())
        .map(|p| p.make_monic())
        .collect();
    if basis.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }
    let mut pairs = BinaryHeap::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            if let Ok(pair) = CriticalPair::new(i, j, &basis[i], &basis[j]) {
                pairs.push(pair);
            } else {
                return Err(GroebnerError::NoLeadingMonomial(i));
            }
        }
    }
    // Optional sugar strategy
    let mut sugar_queue: Vec<SugaredPolynomial<F>> = Vec::new();
    if let SelectionStrategy::Sugar = *strategy {
        // Initialize sugar queue with all S-polynomials
        for i in 0..basis.len() {
            for j in i + 1..basis.len() {
                if let Ok(s_poly) = basis[i].s_polynomial(&basis[j]) {
                    let sugared = SugaredPolynomial::new(s_poly.clone());
                    sugar_queue.push(sugared);
                }
            }
        }
    }
    // Optional Gebauer–Möller criteria
    let mut gm_pairs: Vec<CriticalPair> = Vec::new();
    if let SelectionStrategy::GebauerMoller = *strategy {
        // Start with all pairs
        for i in 0..basis.len() {
            for j in i + 1..basis.len() {
                if let Ok(pair) = CriticalPair::new(i, j, &basis[i], &basis[j]) {
                    gm_pairs.push(pair);
                }
            }
        }
        gm_pairs = grebauer_moller::filter_gm_pairs(&basis, gm_pairs);
    }
    while match strategy {
        SelectionStrategy::Degree => !pairs.is_empty(),
        SelectionStrategy::Sugar => !sugar_queue.is_empty(),
        SelectionStrategy::GebauerMoller => !gm_pairs.is_empty(),
    } {
        let (_poly_i, _poly_j, s_poly) = match strategy {
            SelectionStrategy::Degree => {
                let Some(pair) = pairs.pop() else {
                    break;
                };
                if pair.i >= basis.len() || pair.j >= basis.len() {
                    continue;
                }
                let poly_i = &basis[pair.i];
                let poly_j = &basis[pair.j];
                let lm_i = poly_i
                    .leading_monomial()
                    .ok_or(GroebnerError::NoLeadingMonomial(pair.i))?;
                let lm_j = poly_j
                    .leading_monomial()
                    .ok_or(GroebnerError::NoLeadingMonomial(pair.j))?;
                let product = lm_i.multiply(lm_j);
                if pair.lcm == product {
                    continue;
                }
                let s_poly = poly_i.s_polynomial(poly_j).map_err(GroebnerError::from)?;
                (poly_i.clone(), poly_j.clone(), s_poly)
            }
            SelectionStrategy::Sugar => {
                let Some(sugared) = select_next_by_sugar(&mut sugar_queue) else {
                    break;
                };
                // Find which basis elements produced this S-polynomial (approximate: just use all pairs)
                // This is a simplification; for a more precise implementation, track indices.
                let mut found = false;
                let mut poly_i = None;
                let mut poly_j = None;
                for i in 0..basis.len() {
                    for j in i + 1..basis.len() {
                        if let Ok(s) = basis[i].s_polynomial(&basis[j]) {
                            if s == sugared.poly {
                                poly_i = Some(basis[i].clone());
                                poly_j = Some(basis[j].clone());
                                found = true;
                                break;
                            }
                        }
                    }
                    if found {
                        break;
                    }
                }
                if !found {
                    continue;
                }
                if poly_i.is_none() || poly_j.is_none() {
                    continue;
                } else if let (Some(poly_i), Some(poly_j)) = (poly_i, poly_j) {
                    // Return the S-polynomial and the two polynomials that produced it
                    (poly_i, poly_j, sugared.poly)
                } else {
                    continue; // Should never happen
                }
            }
            SelectionStrategy::GebauerMoller => {
                let Some(pair) = gm_pairs.pop() else {
                    break;
                };
                if pair.i >= basis.len() || pair.j >= basis.len() {
                    continue;
                }
                let poly_i = &basis[pair.i];
                let poly_j = &basis[pair.j];
                let lm_i = poly_i
                    .leading_monomial()
                    .ok_or(GroebnerError::NoLeadingMonomial(pair.i))?;
                let lm_j = poly_j
                    .leading_monomial()
                    .ok_or(GroebnerError::NoLeadingMonomial(pair.j))?;
                let product = lm_i.multiply(lm_j);
                if pair.lcm == product {
                    continue;
                }
                let s_poly = poly_i.s_polynomial(poly_j).map_err(GroebnerError::from)?;
                (poly_i.clone(), poly_j.clone(), s_poly)
            }
        };
        let reduced = s_poly.reduce(&basis).map_err(GroebnerError::from)?;
        if !reduced.is_zero() {
            let monic_reduced = reduced.make_monic();
            let new_index = basis.len();
            for (i, existing) in basis.iter().enumerate() {
                if let Ok(new_pair) = CriticalPair::new(i, new_index, existing, &monic_reduced) {
                    if let SelectionStrategy::Degree = *strategy {
                        pairs.push(new_pair);
                    } else {
                        // For sugar, push new S-polys to sugar_queue
                        if let Ok(s_poly) = existing.s_polynomial(&monic_reduced) {
                            let sugared = SugaredPolynomial::new(s_poly.clone());
                            sugar_queue.push(sugared);
                        }
                    }
                } else {
                    return Err(GroebnerError::NoLeadingMonomial(i));
                }
            }
            basis.push(monic_reduced);
        }
    }
    finish_basis(basis, canonicalize)
}

fn minimize_basis<F: Field>(basis: &mut Vec<Polynomial<F>>) {
    let mut to_remove = Vec::new();
    for i in 0..basis.len() {
        for j in 0..basis.len() {
            if i != j {
                if let (Some(lm_i), Some(lm_j)) =
                    (basis[i].leading_monomial(), basis[j].leading_monomial())
                {
                    if lm_j.divides(lm_i) {
                        to_remove.push(i);
                        break;
                    }
                }
            }
        }
    }
    to_remove.sort_unstable();
    to_remove.reverse();
    for &i in &to_remove {
        basis.remove(i);
    }
}

fn finish_basis<F: Field>(
    mut basis: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    minimize_basis(&mut basis);
    if canonicalize {
        for poly in &mut basis {
            *poly = poly.make_monic();
        }
        basis.sort_by(|a, b| {
            let la = a.leading_monomial();
            let lb = b.leading_monomial();
            match (la, lb) {
                (Some(ma), Some(mb)) => mb.compare(ma, a.order),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            }
        });
        let mut i = 0;
        while i < basis.len() {
            let mut j = i + 1;
            while j < basis.len() {
                if basis[i].terms.len() == basis[j].terms.len()
                    && basis[i]
                        .terms
                        .iter()
                        .zip(&basis[j].terms)
                        .all(|(t1, t2)| t1.monomial == t2.monomial)
                {
                    let ratio = basis[i].terms[0]
                        .coefficient
                        .divide(&basis[j].terms[0].coefficient);
                    if basis[i]
                        .terms
                        .iter()
                        .zip(&basis[j].terms)
                        .all(|(t1, t2)| t1.coefficient.divide(&t2.coefficient) == ratio)
                    {
                        basis.remove(j);
                        continue;
                    }
                }
                j += 1;
            }
            i += 1;
        }
    }
    Ok(basis)
}

#[cfg(feature = "parallel")]
fn all_critical_pairs<F: Field>(
    basis: &[Polynomial<F>],
) -> Result<Vec<CriticalPair>, GroebnerError> {
    let mut pairs = Vec::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            pairs.push(CriticalPair::new(i, j, &basis[i], &basis[j])?);
        }
    }
    Ok(pairs)
}

#[cfg(feature = "parallel")]
fn select_min_degree_critical_pairs(pairs: &mut Vec<CriticalPair>) -> Vec<CriticalPair> {
    let Some(min_degree) = pairs.iter().map(|pair| pair.degree).min() else {
        return Vec::new();
    };

    let mut selected = Vec::new();
    let mut remaining = Vec::new();
    for pair in pairs.drain(..) {
        if pair.degree == min_degree {
            selected.push(pair);
        } else {
            remaining.push(pair);
        }
    }
    *pairs = remaining;
    selected
}

#[cfg(feature = "parallel")]
fn reduce_pair_against_basis<F: Field>(
    pair: &CriticalPair,
    basis: &[Polynomial<F>],
) -> Result<Option<Polynomial<F>>, GroebnerError> {
    if pair.i >= basis.len() || pair.j >= basis.len() {
        return Ok(None);
    }
    let poly_i = &basis[pair.i];
    let poly_j = &basis[pair.j];
    let lm_i = poly_i
        .leading_monomial()
        .ok_or(GroebnerError::NoLeadingMonomial(pair.i))?;
    let lm_j = poly_j
        .leading_monomial()
        .ok_or(GroebnerError::NoLeadingMonomial(pair.j))?;
    if pair.lcm == lm_i.multiply(lm_j) {
        return Ok(None);
    }
    let s_poly = poly_i.s_polynomial(poly_j).map_err(GroebnerError::from)?;
    let reduced = s_poly.reduce(basis).map_err(GroebnerError::from)?;
    Ok((!reduced.is_zero()).then(|| reduced.make_monic()))
}

pub fn is_groebner_basis<F: Field>(basis: &[Polynomial<F>]) -> Result<bool, GroebnerError> {
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            let s_poly = basis[i]
                .s_polynomial(&basis[j])
                .map_err(GroebnerError::from)?;
            let reduced = s_poly.reduce(basis).map_err(GroebnerError::from)?;
            if !reduced.is_zero() {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Check the Buchberger criterion in parallel.
///
/// This API is available with the default `parallel` feature and evaluates independent
/// S-polynomial reductions on Rayon worker threads.
#[cfg(feature = "parallel")]
pub fn is_groebner_basis_parallel<F: Field + Send + Sync>(
    basis: &[Polynomial<F>],
) -> Result<bool, GroebnerError> {
    let pairs = (0..basis.len())
        .flat_map(|i| (i + 1..basis.len()).map(move |j| (i, j)))
        .collect::<Vec<_>>();
    let checks = pairs
        .par_iter()
        .map(|&(i, j)| {
            let s_poly = basis[i]
                .s_polynomial(&basis[j])
                .map_err(GroebnerError::from)?;
            let reduced = s_poly.reduce(basis).map_err(GroebnerError::from)?;
            Ok(reduced.is_zero())
        })
        .collect::<Result<Vec<_>, GroebnerError>>()?;
    Ok(checks.into_iter().all(|is_zero| is_zero))
}
