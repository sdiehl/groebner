//! Multivariate polynomials over a [`Field`].
//!
//! Terms are stored in descending monomial order. Build polynomials with
//! [`crate::PolynomialRing::parse`] rather than by hand.
//!
//! ```
//! use groebner::{MonomialOrder, PolynomialRing};
//! use num_rational::BigRational;
//! let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)?;
//! let p = ring.parse("x^5 - x + 1")?;
//! assert_eq!(p.terms.len(), 3);
//! assert_eq!(p.nvars, 2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::field::Field;
use crate::monomial::{Monomial, MonomialOrder};
use std::cmp::Ordering;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Term<F> {
    pub coefficient: F,
    pub monomial: Monomial,
}

impl<F> Term<F> {
    pub fn new(coefficient: F, monomial: Monomial) -> Self {
        Self {
            coefficient,
            monomial,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<F> {
    pub terms: Vec<Term<F>>,
    pub nvars: usize,
    pub order: MonomialOrder,
}

#[derive(Debug)]
pub enum PolynomialError {
    NoLeadingMonomial,
    NoLeadingCoefficient,
    DivisionFailed,
    DivisionByZero,
}

impl fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolynomialError::NoLeadingMonomial => write!(f, "Polynomial has no leading monomial"),
            PolynomialError::NoLeadingCoefficient => {
                write!(f, "Polynomial has no leading coefficient")
            }
            PolynomialError::DivisionFailed => write!(f, "Division of monomials failed"),
            PolynomialError::DivisionByZero => write!(f, "Division by zero"),
        }
    }
}

impl std::error::Error for PolynomialError {}

impl<F: Field> Polynomial<F> {
    pub fn new(mut terms: Vec<Term<F>>, nvars: usize, order: MonomialOrder) -> Self {
        terms.retain(|t| !t.coefficient.is_zero());
        terms.sort_by(|a, b| b.monomial.compare(&a.monomial, &order));
        let mut combined: Vec<Term<F>> = Vec::with_capacity(terms.len());
        for term in terms {
            push_combined(&mut combined, term);
        }
        Self {
            terms: combined,
            nvars,
            order,
        }
    }

    pub fn zero(nvars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: Vec::new(),
            nvars,
            order,
        }
    }

    pub fn constant(coeff: F, nvars: usize, order: MonomialOrder) -> Self {
        if coeff.is_zero() {
            Self::zero(nvars, order)
        } else {
            Self::new(vec![Term::new(coeff, Monomial::one(nvars))], nvars, order)
        }
    }

    pub fn monomial(monomial: Monomial, order: MonomialOrder) -> Self {
        let nvars = monomial.nvars();
        Self::new(vec![Term::new(F::one(), monomial)], nvars, order)
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
    pub fn is_constant(&self) -> bool {
        self.terms.len() == 1 && self.terms[0].monomial.is_one()
    }
    pub fn leading_term(&self) -> Option<&Term<F>> {
        self.terms.first()
    }
    pub fn leading_monomial(&self) -> Option<&Monomial> {
        self.terms.first().map(|t| &t.monomial)
    }
    pub fn leading_coefficient(&self) -> Option<&F> {
        self.terms.first().map(|t| &t.coefficient)
    }
    pub fn total_degree(&self) -> u32 {
        self.terms
            .iter()
            .map(|t| t.monomial.degree())
            .max()
            .unwrap_or(0)
    }

    /// The same polynomial with its terms sorted under `order`.
    pub fn reorder(&self, order: MonomialOrder) -> Self {
        Self::new(self.terms.clone(), self.nvars, order)
    }

    pub fn map_coefficients<G: Field>(&self, f: impl Fn(&F) -> G) -> Polynomial<G> {
        let terms = self
            .terms
            .iter()
            .map(|t| Term::new(f(&t.coefficient), t.monomial.clone()))
            .collect();
        Polynomial::new(terms, self.nvars, self.order.clone())
    }

    pub fn make_monic(&self) -> Self {
        match self.leading_coefficient().and_then(Field::inverse) {
            Some(inv) => self.multiply_scalar(&inv),
            None => self.clone(),
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        self.merge(other, None, None)
    }

    pub fn subtract(&self, other: &Self) -> Self {
        self.merge(other, Some(&F::one().negate()), None)
    }

    /// `self - scale * monomial * other` in a single merge.
    pub fn subtract_multiple(&self, scale: &F, monomial: &Monomial, other: &Self) -> Self {
        self.merge(other, Some(&scale.negate()), Some(monomial))
    }

    fn merge(&self, other: &Self, scale: Option<&F>, shift: Option<&Monomial>) -> Self {
        debug_assert_eq!(self.nvars, other.nvars);
        let right = |t: &Term<F>| {
            let coefficient = match scale {
                Some(s) => t.coefficient.multiply(s),
                None => t.coefficient.clone(),
            };
            let monomial = match shift {
                Some(m) => t.monomial.multiply(m),
                None => t.monomial.clone(),
            };
            Term::new(coefficient, monomial)
        };
        let mut merged = Vec::with_capacity(self.terms.len() + other.terms.len());
        let (mut i, mut j) = (0, 0);
        while i < self.terms.len() && j < other.terms.len() {
            let r = right(&other.terms[j]);
            match self.terms[i].monomial.compare(&r.monomial, &self.order) {
                Ordering::Greater => {
                    merged.push(self.terms[i].clone());
                    i += 1;
                }
                Ordering::Less => {
                    merged.push(r);
                    j += 1;
                }
                Ordering::Equal => {
                    let c = self.terms[i].coefficient.add(&r.coefficient);
                    if !c.is_zero() {
                        merged.push(Term::new(c, r.monomial));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        merged.extend_from_slice(&self.terms[i..]);
        merged.extend(other.terms[j..].iter().map(right));
        Self {
            terms: merged,
            nvars: self.nvars,
            order: self.order.clone(),
        }
    }

    pub fn multiply_scalar(&self, scalar: &F) -> Self {
        if scalar.is_zero() {
            return Self::zero(self.nvars, self.order.clone());
        }
        let terms = self
            .terms
            .iter()
            .map(|t| Term::new(t.coefficient.multiply(scalar), t.monomial.clone()))
            .collect();
        Self {
            terms,
            nvars: self.nvars,
            order: self.order.clone(),
        }
    }

    pub fn multiply_monomial(&self, monomial: &Monomial) -> Self {
        let terms = self
            .terms
            .iter()
            .map(|t| Term::new(t.coefficient.clone(), t.monomial.multiply(monomial)))
            .collect();
        Self {
            terms,
            nvars: self.nvars,
            order: self.order.clone(),
        }
    }

    pub fn multiply(&self, other: &Self) -> Self {
        let mut acc = Self::zero(self.nvars, self.order.clone());
        for t in &other.terms {
            acc = acc.merge(self, Some(&t.coefficient), Some(&t.monomial));
        }
        acc
    }

    pub fn s_polynomial(&self, other: &Self) -> Result<Self, PolynomialError> {
        if self.is_zero() || other.is_zero() {
            return Ok(Self::zero(self.nvars, self.order.clone()));
        }
        let lt1 = self
            .leading_term()
            .ok_or(PolynomialError::NoLeadingMonomial)?;
        let lt2 = other
            .leading_term()
            .ok_or(PolynomialError::NoLeadingMonomial)?;
        let lcm = lt1.monomial.lcm(&lt2.monomial);
        let m1 = lcm
            .divide(&lt1.monomial)
            .ok_or(PolynomialError::DivisionFailed)?;
        let m2 = lcm
            .divide(&lt2.monomial)
            .ok_or(PolynomialError::DivisionFailed)?;
        let c1 = lt1
            .coefficient
            .inverse()
            .ok_or(PolynomialError::DivisionByZero)?;
        let c2 = lt2
            .coefficient
            .inverse()
            .ok_or(PolynomialError::DivisionByZero)?;
        Ok(self
            .multiply_monomial(&m1)
            .multiply_scalar(&c1)
            .subtract_multiple(&c2, &m2, other))
    }

    /// Full normal form of `self` modulo `basis`: no remaining term is divisible by any leading monomial.
    pub fn reduce(&self, basis: &[Self]) -> Result<Self, PolynomialError> {
        let divisors: Vec<_> = basis
            .iter()
            .filter_map(|g| g.leading_term().map(|lt| (lt, g)))
            .collect();
        let mut work = self.clone();
        let mut head = 0;
        let mut rest = Vec::new();
        while head < work.terms.len() {
            let lt = &work.terms[head];
            let reducer = divisors
                .iter()
                .filter(|(glt, _)| glt.monomial.divides(&lt.monomial))
                .min_by_key(|(_, g)| g.terms.len());
            if let Some((glt, g)) = reducer {
                let m = lt
                    .monomial
                    .divide(&glt.monomial)
                    .ok_or(PolynomialError::DivisionFailed)?;
                let c = lt
                    .coefficient
                    .divide(&glt.coefficient)
                    .ok_or(PolynomialError::DivisionByZero)?;
                let tail = Self {
                    terms: work.terms[head..].to_vec(),
                    nvars: work.nvars,
                    order: work.order.clone(),
                };
                work = tail.subtract_multiple(&c, &m, g);
                head = 0;
            } else {
                rest.push(lt.clone());
                head += 1;
            }
        }
        Ok(Self {
            terms: rest,
            nvars: self.nvars,
            order: self.order.clone(),
        })
    }
}

fn push_combined<F: Field>(terms: &mut Vec<Term<F>>, term: Term<F>) {
    if let Some(last) = terms.last_mut()
        && last.monomial == term.monomial
    {
        last.coefficient = last.coefficient.add(&term.coefficient);
        if last.coefficient.is_zero() {
            terms.pop();
        }
        return;
    }
    terms.push(term);
}

impl<F: Field> fmt::Display for Polynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        for (i, term) in self.terms.iter().enumerate() {
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "{}", term.coefficient)?;
            if !term.monomial.is_one() {
                write!(f, "*{}", term.monomial)?;
            }
        }
        Ok(())
    }
}
