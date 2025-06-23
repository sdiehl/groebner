//! Groebner basis algorithms.
//!
//! This module implements Buchberger's algorithm and related utilities for computing Groebner bases
//! of polynomial ideals. It also provides functions for checking if a set of polynomials forms a Groebner basis.
//!
//! # Algorithms
//! - Buchberger's algorithm (with monic and minimal basis options)
//! - Basis minimization
//! - S-polynomial computation and reduction
//!
//! # Example
//! ```
//! use groebner::{groebner_basis, Polynomial, Rational, MonomialOrder, Term, Monomial};
//! let f1 = Polynomial::new(
//!     vec![
//!         Term::new(Rational::new(1, 1), Monomial::new(vec![2, 0])),
//!         Term::new(Rational::new(-1, 1), Monomial::new(vec![0, 1]))
//!     ],
//!     2,
//!     MonomialOrder::Lex
//! );
//! let f2 = Polynomial::new(
//!     vec![
//!         Term::new(Rational::new(1, 1), Monomial::new(vec![1, 1])),
//!         Term::new(Rational::new(-1, 1), Monomial::new(vec![0, 0]))
//!     ],
//!     2,
//!     MonomialOrder::Lex
//! );
//! let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true);
//! assert!(!basis.is_empty());
//! ```

use crate::field::Field;
use crate::monomial::Monomial;
use crate::polynomial::Polynomial;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Debug, Clone)]
struct CriticalPair {
    i: usize,
    j: usize,
    lcm: Monomial,
    degree: u32,
}

impl CriticalPair {
    fn new(
        i: usize,
        j: usize,
        poly_i: &Polynomial<impl Field>,
        poly_j: &Polynomial<impl Field>,
    ) -> Self {
        let lm_i = poly_i.leading_monomial().unwrap();
        let lm_j = poly_j.leading_monomial().unwrap();
        let lcm = lm_i.lcm(lm_j);
        let degree = lcm.degree();
        Self { i, j, lcm, degree }
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

/// Compute Groebner basis using Buchberger's algorithm
///
/// # Example
/// ```
/// use groebner::{groebner_basis, Polynomial, Rational, MonomialOrder};
/// // x^2 - y, xy - 1
/// let f1 = Polynomial::new(
///     vec![
///         groebner::Term::new(Rational::new(1, 1), groebner::Monomial::new(vec![2, 0])),
///         groebner::Term::new(Rational::new(-1, 1), groebner::Monomial::new(vec![0, 1]))
///     ],
///     2,
///     MonomialOrder::Lex
/// );
/// let f2 = Polynomial::new(
///     vec![
///         groebner::Term::new(Rational::new(1, 1), groebner::Monomial::new(vec![1, 1])),
///         groebner::Term::new(Rational::new(-1, 1), groebner::Monomial::new(vec![0, 0]))
///     ],
///     2,
///     MonomialOrder::Lex
/// );
/// let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true);
/// assert!(!basis.is_empty());
/// ```
#[allow(clippy::needless_range_loop)]
pub fn groebner_basis<F: Field>(
    polynomials: Vec<Polynomial<F>>,
    _order: crate::monomial::MonomialOrder,
    canonicalize: bool,
) -> Vec<Polynomial<F>> {
    if polynomials.is_empty() {
        return Vec::new();
    }
    let _nvars = polynomials[0].nvars;
    let mut basis: Vec<Polynomial<F>> = polynomials
        .into_iter()
        .filter(|p| !p.is_zero())
        .map(|p| p.make_monic())
        .collect();
    if basis.is_empty() {
        return Vec::new();
    }
    let mut pairs = BinaryHeap::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            pairs.push(CriticalPair::new(i, j, &basis[i], &basis[j]));
        }
    }
    while let Some(pair) = pairs.pop() {
        if pair.i >= basis.len() || pair.j >= basis.len() {
            continue;
        }
        let poly_i = &basis[pair.i];
        let poly_j = &basis[pair.j];
        let lm_i = poly_i.leading_monomial().unwrap();
        let lm_j = poly_j.leading_monomial().unwrap();
        let product = lm_i.multiply(lm_j);
        if pair.lcm == product {
            continue;
        }
        let s_poly = poly_i.s_polynomial(poly_j);
        let reduced = s_poly.reduce(&basis);
        if !reduced.is_zero() {
            let monic_reduced = reduced.make_monic();
            let new_index = basis.len();
            for (i, existing) in basis.iter().enumerate() {
                pairs.push(CriticalPair::new(i, new_index, existing, &monic_reduced));
            }
            basis.push(monic_reduced);
        }
    }
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
    basis
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

pub fn is_groebner_basis<F: Field>(basis: &[Polynomial<F>]) -> bool {
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            let s_poly = basis[i].s_polynomial(&basis[j]);
            let reduced = s_poly.reduce(basis);
            if !reduced.is_zero() {
                return false;
            }
        }
    }
    true
}
