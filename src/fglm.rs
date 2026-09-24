//! FGLM change of ordering for zero-dimensional ideals.
//!
//! ```
//! use groebner::{fglm, groebner_basis_f4, MonomialOrder, PolynomialRing, PrimeField};
//!
//! let ring = PolynomialRing::<PrimeField<32003>>::new(["x", "y"], MonomialOrder::GRevLex)?;
//! let polys = ring.parse_many("x^2 + y^2 - 1; x - y^3")?;
//! let grevlex = groebner_basis_f4(polys, true)?;
//! let lex = fglm(&grevlex, &MonomialOrder::Lex)?;
//! assert!(lex.iter().any(|p| p.terms.iter().all(|t| t.monomial.exponents()[0] == 0)));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::field::Field;
use crate::groebner::{compare_leading, GroebnerError};
use crate::monomial::{Monomial, MonomialOrder};
use crate::polynomial::{Polynomial, Term};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet};

/// True when the leading monomials of `basis` contain a pure power of every variable.
pub fn is_zero_dimensional<F: Field>(basis: &[Polynomial<F>]) -> bool {
    let Some(nvars) = basis.first().map(|p| p.nvars) else {
        return false;
    };
    let leads: Vec<&Monomial> = basis
        .iter()
        .filter_map(Polynomial::leading_monomial)
        .collect();
    if leads.iter().any(|m| m.is_one()) {
        return true;
    }
    (0..nvars).all(|i| {
        leads.iter().any(|m| {
            m.exponents()
                .iter()
                .enumerate()
                .all(|(j, e)| j == i || *e == 0)
        })
    })
}

/// Monomials not divisible by any leading monomial of `basis`, or `None` if there are infinitely many.
pub fn standard_monomials<F: Field>(basis: &[Polynomial<F>]) -> Option<Vec<Monomial>> {
    if !is_zero_dimensional(basis) {
        return None;
    }
    let nvars = basis.first()?.nvars;
    let leads: Vec<&Monomial> = basis
        .iter()
        .filter_map(Polynomial::leading_monomial)
        .collect();
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    let mut stack = vec![Monomial::one(nvars)];
    while let Some(m) = stack.pop() {
        if !seen.insert(m.clone()) || leads.iter().any(|l| l.divides(&m)) {
            continue;
        }
        for i in 0..nvars {
            stack.push(m.multiply(&Monomial::variable(i, nvars)));
        }
        out.push(m);
    }
    Some(out)
}

struct Ordered {
    monomial: Monomial,
    order: MonomialOrder,
}

impl PartialEq for Ordered {
    fn eq(&self, other: &Self) -> bool {
        self.monomial == other.monomial
    }
}
impl Eq for Ordered {}
impl PartialOrd for Ordered {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Ordered {
    fn cmp(&self, other: &Self) -> Ordering {
        other.monomial.compare(&self.monomial, &self.order)
    }
}

/// Convert a reduced Groebner basis of a zero-dimensional ideal to the reduced basis under `target`.
pub fn fglm<F: Field>(
    basis: &[Polynomial<F>],
    target: &MonomialOrder,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    if basis.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }
    if !is_zero_dimensional(basis) {
        return Err(GroebnerError::NotZeroDimensional);
    }
    let nvars = basis[0].nvars;
    let mut index: HashMap<Monomial, usize> = HashMap::new();
    let mut rows: Vec<(Vec<F>, Vec<F>)> = Vec::new();
    let mut staircase: Vec<(Monomial, Polynomial<F>)> = Vec::new();
    let mut new_leads: Vec<Monomial> = Vec::new();
    let mut result: Vec<Polynomial<F>> = Vec::new();
    let mut processed: HashSet<Monomial> = HashSet::new();
    let mut queue = BinaryHeap::new();
    queue.push(Ordered {
        monomial: Monomial::one(nvars),
        order: target.clone(),
    });

    while let Some(Ordered { monomial: m, .. }) = queue.pop() {
        if !processed.insert(m.clone()) || new_leads.iter().any(|l| l.divides(&m)) {
            continue;
        }
        let nf = normal_form_of(&m, &staircase, basis)?;
        let mut vector = vec![F::zero(); index.len()];
        for t in &nf.terms {
            let next = index.len();
            let k = *index.entry(t.monomial.clone()).or_insert(next);
            if k >= vector.len() {
                vector.resize(k + 1, F::zero());
            }
            vector[k] = t.coefficient.clone();
        }
        let mut combination = vec![F::zero(); staircase.len() + 1];
        combination[staircase.len()] = F::one();
        for (row, row_combination) in &rows {
            let Some(pivot) = row.iter().position(|v| !v.is_zero()) else {
                continue;
            };
            if pivot >= vector.len() || vector[pivot].is_zero() {
                continue;
            }
            let c = vector[pivot].clone();
            for (k, v) in row.iter().enumerate() {
                vector[k] = vector[k].subtract(&c.multiply(v));
            }
            for (k, v) in row_combination.iter().enumerate() {
                combination[k] = combination[k].subtract(&c.multiply(v));
            }
        }
        if vector.iter().all(Field::is_zero) {
            let mut terms = vec![Term::new(F::one(), m.clone())];
            for (k, (b, _)) in staircase.iter().enumerate() {
                if !combination[k].is_zero() {
                    terms.push(Term::new(combination[k].clone(), b.clone()));
                }
            }
            result.push(Polynomial::new(terms, nvars, target.clone()));
            new_leads.push(m);
        } else {
            let Some(pivot) = vector.iter().position(|v| !v.is_zero()) else {
                continue;
            };
            let inv = vector[pivot]
                .inverse()
                .ok_or(crate::polynomial::PolynomialError::DivisionByZero)?;
            let vector: Vec<F> = vector.iter().map(|v| v.multiply(&inv)).collect();
            let combination: Vec<F> = combination.iter().map(|v| v.multiply(&inv)).collect();
            rows.push((vector, combination));
            for i in 0..nvars {
                queue.push(Ordered {
                    monomial: m.multiply(&Monomial::variable(i, nvars)),
                    order: target.clone(),
                });
            }
            staircase.push((m, nf));
        }
    }
    result.sort_by(|a, b| compare_leading(b, a, target));
    Ok(result)
}

fn normal_form_of<F: Field>(
    m: &Monomial,
    staircase: &[(Monomial, Polynomial<F>)],
    basis: &[Polynomial<F>],
) -> Result<Polynomial<F>, GroebnerError> {
    let order = basis[0].order.clone();
    for (b, nf) in staircase.iter().rev() {
        if let Some(q) = m.divide(b) {
            if q.degree() == 1 {
                return Ok(nf.multiply_monomial(&q).reduce(basis)?);
            }
        }
    }
    Ok(Polynomial::monomial(m.clone(), order).reduce(basis)?)
}
