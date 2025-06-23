//! Multivariate polynomial types and operations.
//!
//! This module defines the `Polynomial` struct and related types for representing and manipulating
//! multivariate polynomials over arbitrary fields. It supports arithmetic, reduction, and S-polynomial
//! computation, which are fundamental for Groebner basis algorithms.
//!
//! # Example
//! ```
//! use num_rational::BigRational;
//! use groebner::{Polynomial, Term, Monomial, MonomialOrder};
//! let p1 = Polynomial::new(
//!     vec![Term::new(BigRational::new(2.into(), 1.into()), Monomial::new(vec![1, 0]))],
//!     2,
//!     MonomialOrder::Lex
//! );
//! let p2 = Polynomial::new(
//!     vec![Term::new(BigRational::new(3.into(), 1.into()), Monomial::new(vec![0, 1]))],
//!     2,
//!     MonomialOrder::Lex
//! );
//! let sum = p1.add(&p2);
//! assert_eq!(sum.terms.len(), 2);
//! ```

use crate::field::Field;
use crate::monomial::{Monomial, MonomialOrder};
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
}

impl fmt::Display for PolynomialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolynomialError::NoLeadingMonomial => write!(f, "Polynomial has no leading monomial"),
            PolynomialError::NoLeadingCoefficient => {
                write!(f, "Polynomial has no leading coefficient")
            }
            PolynomialError::DivisionFailed => write!(f, "Division of monomials failed"),
        }
    }
}

impl std::error::Error for PolynomialError {}

impl<F: Field> Polynomial<F> {
    pub fn new(mut terms: Vec<Term<F>>, nvars: usize, order: MonomialOrder) -> Self {
        terms.retain(|t| !t.coefficient.is_zero());
        terms.sort_by(|a, b| b.monomial.compare(&a.monomial, order));
        let mut combined: Vec<Term<F>> = Vec::new();
        for term in terms {
            if let Some(last) = combined.last_mut() {
                if last.monomial == term.monomial {
                    last.coefficient = last.coefficient.add(&term.coefficient);
                    if last.coefficient.is_zero() {
                        combined.pop();
                    }
                    continue;
                }
            }
            combined.push(term);
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
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
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
    pub fn make_monic(&self) -> Self {
        if let Some(lc) = self.leading_coefficient() {
            if !lc.is_zero() {
                let inv_lc = lc.inverse();
                return self.multiply_scalar(&inv_lc);
            }
        }
        self.clone()
    }
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);
        let mut terms = self.terms.clone();
        terms.extend(other.terms.clone());
        Self::new(terms, self.nvars, self.order)
    }
    pub fn subtract(&self, other: &Self) -> Self {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);
        let mut terms = self.terms.clone();
        for term in &other.terms {
            terms.push(Term::new(term.coefficient.negate(), term.monomial.clone()));
        }
        Self::new(terms, self.nvars, self.order)
    }
    pub fn multiply_scalar(&self, scalar: &F) -> Self {
        if scalar.is_zero() {
            return Self::zero(self.nvars, self.order);
        }
        let terms = self
            .terms
            .iter()
            .map(|t| Term::new(t.coefficient.multiply(scalar), t.monomial.clone()))
            .collect();
        Self::new(terms, self.nvars, self.order)
    }
    pub fn multiply_monomial(&self, monomial: &Monomial) -> Self {
        let terms = self
            .terms
            .iter()
            .map(|t| Term::new(t.coefficient.clone(), t.monomial.multiply(monomial)))
            .collect();
        Self::new(terms, self.nvars, self.order)
    }
    pub fn multiply(&self, other: &Self) -> Self {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);
        if self.is_zero() || other.is_zero() {
            return Self::zero(self.nvars, self.order);
        }
        let mut terms = Vec::new();
        for t1 in &self.terms {
            for t2 in &other.terms {
                terms.push(Term::new(
                    t1.coefficient.multiply(&t2.coefficient),
                    t1.monomial.multiply(&t2.monomial),
                ));
            }
        }
        Self::new(terms, self.nvars, self.order)
    }
    pub fn s_polynomial(&self, other: &Self) -> Result<Self, PolynomialError> {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);
        if self.is_zero() || other.is_zero() {
            return Ok(Self::zero(self.nvars, self.order));
        }
        let lm1 = self
            .leading_monomial()
            .ok_or(PolynomialError::NoLeadingMonomial)?;
        let lm2 = other
            .leading_monomial()
            .ok_or(PolynomialError::NoLeadingMonomial)?;
        let lc1 = self
            .leading_coefficient()
            .ok_or(PolynomialError::NoLeadingCoefficient)?;
        let lc2 = other
            .leading_coefficient()
            .ok_or(PolynomialError::NoLeadingCoefficient)?;
        let lcm = lm1.lcm(lm2);
        let m1 = lcm.divide(lm1).ok_or(PolynomialError::DivisionFailed)?;
        let m2 = lcm.divide(lm2).ok_or(PolynomialError::DivisionFailed)?;
        let term1 = self.multiply_monomial(&m1).multiply_scalar(&lc2.inverse());
        let term2 = other.multiply_monomial(&m2).multiply_scalar(&lc1.inverse());
        Ok(term1.subtract(&term2))
    }
    pub fn reduce(&self, basis: &[Self]) -> Result<Self, PolynomialError> {
        let mut remainder = self.clone();
        while !remainder.is_zero() {
            let mut reduced = false;
            if let Some(leading_mono) = remainder.leading_monomial() {
                for divisor in basis {
                    if let Some(div_leading) = divisor.leading_monomial() {
                        if div_leading.divides(leading_mono) {
                            let quotient_mono = leading_mono
                                .divide(div_leading)
                                .ok_or(PolynomialError::DivisionFailed)?;
                            let lc_remainder = remainder
                                .leading_coefficient()
                                .ok_or(PolynomialError::NoLeadingCoefficient)?;
                            let lc_divisor = divisor
                                .leading_coefficient()
                                .ok_or(PolynomialError::NoLeadingCoefficient)?;
                            let quotient_coeff = lc_remainder.multiply(&lc_divisor.inverse());
                            let subtrahend = divisor
                                .multiply_monomial(&quotient_mono)
                                .multiply_scalar(&quotient_coeff);
                            remainder = remainder.subtract(&subtrahend);
                            reduced = true;
                            break;
                        }
                    }
                }
            }
            if !reduced {
                break;
            }
        }
        Ok(remainder)
    }
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
            if !term.monomial.exponents.iter().all(|&e| e == 0) {
                write!(f, "*{}", term.monomial)?;
            }
        }
        Ok(())
    }
}
