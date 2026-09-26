//! Ideal membership certificates via Buchberger's algorithm with cofactor tracking.
//!
//! ```
//! use groebner::{LiftBasis, MonomialOrder, PolynomialRing, verify_lift};
//! use num_rational::BigRational;
//!
//! let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::GRevLex)?;
//! let generators = ring.parse_many("x^2 - y; y^2 - x")?;
//! let f = ring.parse("x^4 - x")?;
//! let h = LiftBasis::new(&generators)?.lift(&f)?.ok_or("not a member")?;
//! assert!(verify_lift(&generators, &h, &f));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::field::Field;
use crate::groebner::{CriticalPair, GroebnerError, compare_leading, minimal_mask};
use crate::monomial::Monomial;
use crate::polynomial::{Polynomial, PolynomialError};
use std::collections::BinaryHeap;

type Row<F> = Vec<Polynomial<F>>;

/// Reduced Groebner basis of generators `F` with `basis[j] = sum_k cofactors[j][k] * F[k]`.
#[derive(Debug, Clone, PartialEq)]
pub struct LiftBasis<F> {
    pub basis: Vec<Polynomial<F>>,
    pub cofactors: Vec<Row<F>>,
}

impl<F: Field> LiftBasis<F> {
    /// Extended Buchberger over `generators`, which may include zeros and need not be monic.
    pub fn new(generators: &[Polynomial<F>]) -> Result<Self, GroebnerError> {
        let first = generators
            .iter()
            .find(|p| !p.is_zero())
            .ok_or(GroebnerError::EmptyInput)?;
        let (nvars, order) = (first.nvars, first.order.clone());
        if generators
            .iter()
            .any(|p| p.nvars != nvars || p.order != order)
        {
            return Err(GroebnerError::OrderMismatch);
        }
        let zero = Polynomial::zero(nvars, order.clone());
        let mut basis = Vec::new();
        let mut rows = Vec::new();
        for (k, f) in generators.iter().enumerate().filter(|(_, f)| !f.is_zero()) {
            let mut row = vec![zero.clone(); generators.len()];
            row[k] = Polynomial::constant(F::one(), nvars, order.clone());
            let (p, row) = monic(f.clone(), row);
            basis.push(p);
            rows.push(row);
        }
        let mut pairs = BinaryHeap::new();
        for j in 0..basis.len() {
            for i in 0..j {
                pairs.push(CriticalPair::new(i, j, &basis[i], &basis[j])?);
            }
        }
        while let Some(pair) = pairs.pop() {
            let (gi, gj) = (&basis[pair.i], &basis[pair.j]);
            let (Some(lm_i), Some(lm_j)) = (gi.leading_monomial(), gj.leading_monomial()) else {
                continue;
            };
            if lm_i.is_coprime(lm_j) {
                continue;
            }
            let mi = cofactor(&pair.lcm, lm_i)?;
            let mj = cofactor(&pair.lcm, lm_j)?;
            let s = shift_subtract(gi, &mi, gj, &mj);
            let row: Row<F> = rows[pair.i]
                .iter()
                .zip(&rows[pair.j])
                .map(|(a, b)| shift_subtract(a, &mi, b, &mj))
                .collect();
            let (r, row) = reduce_tracked(&s, row, &basis, &rows)?;
            if r.is_zero() {
                continue;
            }
            for (i, g) in basis.iter().enumerate() {
                pairs.push(CriticalPair::new(i, basis.len(), g, &r)?);
            }
            basis.push(r);
            rows.push(row);
        }

        let keep = minimal_mask(&basis);
        let (basis, rows): (Vec<_>, Vec<_>) = basis
            .into_iter()
            .zip(rows)
            .zip(keep)
            .filter_map(|(entry, k)| k.then_some(entry))
            .unzip();
        let mut divisors = basis.clone();
        let mut reduced = Vec::with_capacity(basis.len());
        for (i, (p, row)) in basis.iter().zip(&rows).enumerate() {
            let saved = std::mem::replace(&mut divisors[i], zero.clone());
            let (r, row) = reduce_tracked(p, row.clone(), &divisors, &rows)?;
            divisors[i] = saved;
            if !r.is_zero() {
                reduced.push((r, row));
            }
        }
        reduced.sort_by(|a, b| compare_leading(&b.0, &a.0, &order));
        let (basis, cofactors) = reduced.into_iter().unzip();
        Ok(Self { basis, cofactors })
    }

    /// Cofactors `h` with `f = sum h[k] * F[k]` over the generators, or `None` if `f` is not a member.
    pub fn lift(&self, f: &Polynomial<F>) -> Result<Option<Row<F>>, GroebnerError> {
        let (Some(g), Some(row)) = (self.basis.first(), self.cofactors.first()) else {
            return Err(GroebnerError::EmptyInput);
        };
        let (quotients, remainder) = f.reorder(g.order.clone()).divide(&self.basis)?;
        if !remainder.is_zero() {
            return Ok(None);
        }
        let zeros = vec![Polynomial::zero(g.nvars, g.order.clone()); row.len()];
        Ok(Some(subtract_combination(
            zeros,
            &quotients,
            &self.cofactors,
            &F::one().negate(),
        )))
    }
}

/// Check `f = sum cofactors[k] * generators[k]` using only ring addition and multiplication.
pub fn verify_lift<F: Field>(
    generators: &[Polynomial<F>],
    cofactors: &[Polynomial<F>],
    f: &Polynomial<F>,
) -> bool {
    generators.len() == cofactors.len()
        && generators
            .iter()
            .chain(cofactors)
            .all(|p| p.nvars == f.nvars && p.order == f.order)
        && generators
            .iter()
            .zip(cofactors)
            .fold(f.clone(), |acc, (g, h)| acc.subtract(&h.multiply(g)))
            .is_zero()
}

fn cofactor(lcm: &Monomial, lead: &Monomial) -> Result<Monomial, PolynomialError> {
    lcm.divide(lead).ok_or(PolynomialError::DivisionFailed)
}

fn shift_subtract<F: Field>(
    a: &Polynomial<F>,
    ma: &Monomial,
    b: &Polynomial<F>,
    mb: &Monomial,
) -> Polynomial<F> {
    a.multiply_monomial(ma).subtract_multiple(&F::one(), mb, b)
}

/// `acc - sign * sum_j quotients[j] * rows[j]`, componentwise.
fn subtract_combination<F: Field>(
    mut acc: Row<F>,
    quotients: &[Polynomial<F>],
    rows: &[Row<F>],
    sign: &F,
) -> Row<F> {
    for (q, row) in quotients.iter().zip(rows) {
        for t in &q.terms {
            let c = t.coefficient.multiply(sign);
            for (a, r) in acc.iter_mut().zip(row) {
                if !r.is_zero() {
                    *a = a.subtract_multiple(&c, &t.monomial, r);
                }
            }
        }
    }
    acc
}

fn reduce_tracked<F: Field>(
    p: &Polynomial<F>,
    row: Row<F>,
    basis: &[Polynomial<F>],
    rows: &[Row<F>],
) -> Result<(Polynomial<F>, Row<F>), PolynomialError> {
    let (quotients, r) = p.divide(basis)?;
    Ok(monic(
        r,
        subtract_combination(row, &quotients, rows, &F::one()),
    ))
}

fn monic<F: Field>(p: Polynomial<F>, row: Row<F>) -> (Polynomial<F>, Row<F>) {
    match p.leading_coefficient().and_then(Field::inverse) {
        Some(inv) => (
            p.multiply_scalar(&inv),
            row.iter().map(|h| h.multiply_scalar(&inv)).collect(),
        ),
        None => (p, row),
    }
}
