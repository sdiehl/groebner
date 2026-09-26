//! Buchberger's algorithm and basis utilities.
//!
//! ```
//! use groebner::{groebner_basis, groebner_basis_incremental, MonomialOrder, PolynomialRing};
//! use num_rational::BigRational;
//!
//! let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)?;
//! let f1 = ring.parse("x^2 - y")?;
//! let f2 = ring.parse("x*y - 1")?;
//! let basis = groebner_basis(vec![f1, f2], true)?;
//! assert!(!basis.is_empty());
//!
//! let f3 = ring.parse("y^2 - x")?;
//! let updated = groebner_basis_incremental(basis, vec![f3], true)?;
//! assert!(!updated.is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::field::Field;
use crate::grebauer_moller;
use crate::monomial::Monomial;
use crate::polynomial::Polynomial;
use crate::sugar::{SugaredPolynomial, select_next_by_sugar};
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
    OrderMismatch,
    NotZeroDimensional,
    NoModulus,
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
            GroebnerError::OrderMismatch => {
                write!(
                    f,
                    "polynomials do not share a monomial order and variable count"
                )
            }
            GroebnerError::NotZeroDimensional => {
                write!(f, "ideal is not zero-dimensional")
            }
            GroebnerError::NoModulus => write!(f, "no modulus found in coefficients"),
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

/// S-polynomial selection strategy for [`groebner_basis_with_strategy`].
pub enum SelectionStrategy {
    Degree,
    Sugar,
    GebauerMoller,
}

/// Compute a Groebner basis with Buchberger's algorithm under the polynomials' own order.
///
/// With `canonicalize` the result is the unique reduced Groebner basis, sorted by leading monomial.
pub fn groebner_basis<F: Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    groebner_basis_with_strategy(polynomials, canonicalize, &SelectionStrategy::Degree)
}

/// Buchberger's algorithm with each minimal-degree batch of pairs reduced in parallel.
#[cfg(feature = "parallel")]
pub fn groebner_basis_parallel<F: Field + Send + Sync>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    let mut basis = prepare_input(polynomials)?;
    let order = basis[0].order.clone();
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

        new_polynomials.sort_by(|a, b| compare_leading(b, a, &order));
        new_polynomials.dedup_by(|a, b| a.leading_monomial() == b.leading_monomial());

        for polynomial in new_polynomials {
            let reduced = polynomial.reduce(&basis)?;
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

/// Extend an existing Groebner basis with new generators.
pub fn groebner_basis_incremental<F: Field>(
    existing_basis: Vec<Polynomial<F>>,
    new_polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
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
            polynomial.reduce(&combined)?
        };
        if !reduced.is_zero() {
            combined.push(reduced);
        }
    }
    groebner_basis_with_strategy(combined, canonicalize, &SelectionStrategy::Degree)
}

/// Buchberger's algorithm with a selectable pair selection strategy.
#[allow(clippy::needless_range_loop)]
pub fn groebner_basis_with_strategy<F: Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
    strategy: &SelectionStrategy,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    let mut basis = prepare_input(polynomials)?;
    let mut pairs = BinaryHeap::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            pairs.push(CriticalPair::new(i, j, &basis[i], &basis[j])?);
        }
    }
    let mut sugar_queue: Vec<SugaredPolynomial<F>> = Vec::new();
    if let SelectionStrategy::Sugar = *strategy {
        for i in 0..basis.len() {
            for j in i + 1..basis.len() {
                if let Ok(s_poly) = basis[i].s_polynomial(&basis[j]) {
                    sugar_queue.push(SugaredPolynomial::new(s_poly));
                }
            }
        }
    }
    let mut gm_pairs: Vec<CriticalPair> = Vec::new();
    if let SelectionStrategy::GebauerMoller = *strategy {
        for i in 0..basis.len() {
            for j in i + 1..basis.len() {
                gm_pairs.push(CriticalPair::new(i, j, &basis[i], &basis[j])?);
            }
        }
        gm_pairs = grebauer_moller::filter_gm_pairs(&basis, gm_pairs);
    }
    while match strategy {
        SelectionStrategy::Degree => !pairs.is_empty(),
        SelectionStrategy::Sugar => !sugar_queue.is_empty(),
        SelectionStrategy::GebauerMoller => !gm_pairs.is_empty(),
    } {
        let s_poly = match strategy {
            SelectionStrategy::Degree => {
                let Some(pair) = pairs.pop() else {
                    break;
                };
                match pair_s_polynomial(&pair, &basis)? {
                    Some(s_poly) => s_poly,
                    None => continue,
                }
            }
            SelectionStrategy::Sugar => {
                let Some(sugared) = select_next_by_sugar(&mut sugar_queue) else {
                    break;
                };
                sugared.poly
            }
            SelectionStrategy::GebauerMoller => {
                let Some(pair) = gm_pairs.pop() else {
                    break;
                };
                match pair_s_polynomial(&pair, &basis)? {
                    Some(s_poly) => s_poly,
                    None => continue,
                }
            }
        };
        let reduced = s_poly.reduce(&basis)?;
        if reduced.is_zero() {
            continue;
        }
        let monic_reduced = reduced.make_monic();
        let new_index = basis.len();
        for (i, existing) in basis.iter().enumerate() {
            let new_pair = CriticalPair::new(i, new_index, existing, &monic_reduced)?;
            match strategy {
                SelectionStrategy::Degree => pairs.push(new_pair),
                SelectionStrategy::GebauerMoller => gm_pairs.push(new_pair),
                SelectionStrategy::Sugar => {
                    if let Ok(s_poly) = existing.s_polynomial(&monic_reduced) {
                        sugar_queue.push(SugaredPolynomial::new(s_poly));
                    }
                }
            }
        }
        basis.push(monic_reduced);
    }
    finish_basis(basis, canonicalize)
}

pub(crate) fn prepare_input<F: Field>(
    polynomials: Vec<Polynomial<F>>,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    let basis: Vec<Polynomial<F>> = polynomials
        .into_iter()
        .filter(|p| !p.is_zero())
        .map(|p| p.make_monic())
        .collect();
    if basis.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }
    let (nvars, order) = (basis[0].nvars, &basis[0].order);
    if basis.iter().any(|p| p.nvars != nvars || &p.order != order) {
        return Err(GroebnerError::OrderMismatch);
    }
    Ok(basis)
}

fn pair_s_polynomial<F: Field>(
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
    if lm_i.is_coprime(lm_j) {
        return Ok(None);
    }
    Ok(Some(poly_i.s_polynomial(poly_j)?))
}

pub(crate) fn compare_leading<F: Field>(
    a: &Polynomial<F>,
    b: &Polynomial<F>,
    order: &crate::monomial::MonomialOrder,
) -> Ordering {
    match (a.leading_monomial(), b.leading_monomial()) {
        (Some(ma), Some(mb)) => ma.compare(mb, order),
        (Some(_), None) => Ordering::Greater,
        (None, Some(_)) => Ordering::Less,
        (None, None) => Ordering::Equal,
    }
}

fn minimize_basis<F: Field>(basis: &mut Vec<Polynomial<F>>) {
    let leads: Vec<_> = basis
        .iter()
        .filter_map(|p| p.leading_monomial().cloned())
        .collect();
    let mut keep = vec![true; basis.len()];
    for i in 0..leads.len() {
        for j in 0..leads.len() {
            if i != j && keep[j] && leads[j].divides(&leads[i]) && (leads[j] != leads[i] || j < i) {
                keep[i] = false;
                break;
            }
        }
    }
    let mut index = 0;
    basis.retain(|_| {
        index += 1;
        keep[index - 1]
    });
}

/// Minimize and, when `canonicalize` is set, interreduce and sort a Groebner basis.
pub(crate) fn finish_basis<F: Field>(
    mut basis: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    basis.retain(|p| !p.is_zero());
    minimize_basis(&mut basis);
    if !canonicalize {
        return Ok(basis);
    }
    let order = basis
        .first()
        .map(|p| p.order.clone())
        .ok_or(GroebnerError::EmptyInput)?;
    let mut reduced = Vec::with_capacity(basis.len());
    for i in 0..basis.len() {
        let others: Vec<_> = basis
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, p)| p.clone())
            .collect();
        let r = basis[i].reduce(&others)?.make_monic();
        if !r.is_zero() {
            reduced.push(r);
        }
    }
    reduced.sort_by(|a, b| compare_leading(b, a, &order));
    Ok(reduced)
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
    let (selected, remaining) = pairs.drain(..).partition(|pair| pair.degree == min_degree);
    *pairs = remaining;
    selected
}

#[cfg(feature = "parallel")]
fn reduce_pair_against_basis<F: Field>(
    pair: &CriticalPair,
    basis: &[Polynomial<F>],
) -> Result<Option<Polynomial<F>>, GroebnerError> {
    let Some(s_poly) = pair_s_polynomial(pair, basis)? else {
        return Ok(None);
    };
    let reduced = s_poly.reduce(basis)?;
    Ok((!reduced.is_zero()).then(|| reduced.make_monic()))
}

/// Check Buchberger's criterion: every S-polynomial reduces to zero.
pub fn is_groebner_basis<F: Field>(basis: &[Polynomial<F>]) -> Result<bool, GroebnerError> {
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            let reduced = basis[i].s_polynomial(&basis[j])?.reduce(basis)?;
            if !reduced.is_zero() {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Check Buchberger's criterion with S-polynomial reductions run in parallel.
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
            let reduced = basis[i].s_polynomial(&basis[j])?.reduce(basis)?;
            Ok(reduced.is_zero())
        })
        .collect::<Result<Vec<_>, GroebnerError>>()?;
    Ok(checks.into_iter().all(|is_zero| is_zero))
}
