//! Groebner-specific operations on polycore polynomials.
use crate::monomial::{MonomialExt, divisibility_mask};
use crate::{Field, Monomial};
pub use polycore::{Poly as Polynomial, Term};
use std::fmt;

/// Construct a term with the coefficient first, for migration from 0.3.
pub fn term<F>(coefficient: F, monomial: Monomial) -> Term<F> {
    (monomial, coefficient)
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

/// Groebner algorithms and compatibility helpers for [`Polynomial`].
pub trait PolynomialExt<F: Field> {
    fn leading_term(&self) -> Option<&Term<F>>;
    fn leading_monomial(&self) -> Option<&Monomial>;
    fn leading_coefficient(&self) -> Option<&F>;
    fn make_monic(&self) -> Polynomial<F>;
    fn map_coefficients<G: Field>(&self, f: impl Fn(&F) -> G) -> Polynomial<G>;
    fn add(&self, other: &Polynomial<F>) -> Polynomial<F>;
    fn subtract(&self, other: &Polynomial<F>) -> Polynomial<F>;
    fn multiply(&self, other: &Polynomial<F>) -> Polynomial<F>;
    fn multiply_scalar(&self, c: &F) -> Polynomial<F>;
    fn multiply_monomial(&self, m: &Monomial) -> Polynomial<F>;
    fn subtract_multiple(&self, c: &F, m: &Monomial, other: &Polynomial<F>) -> Polynomial<F>;
    fn s_polynomial(&self, other: &Polynomial<F>) -> Result<Polynomial<F>, PolynomialError>;
    fn normal_form(&self, basis: &[Polynomial<F>]) -> Result<Polynomial<F>, PolynomialError>;
    fn divide_with_remainder(
        &self,
        basis: &[Polynomial<F>],
    ) -> Result<(Vec<Polynomial<F>>, Polynomial<F>), PolynomialError>;
}
impl<F: Field> PolynomialExt<F> for Polynomial<F> {
    fn leading_term(&self) -> Option<&Term<F>> {
        self.lt()
    }
    fn leading_monomial(&self) -> Option<&Monomial> {
        self.lm()
    }
    fn leading_coefficient(&self) -> Option<&F> {
        self.lc()
    }
    fn make_monic(&self) -> Self {
        self.monic()
    }
    fn map_coefficients<G: Field>(&self, f: impl Fn(&F) -> G) -> Polynomial<G> {
        self.map(f)
    }
    fn add(&self, other: &Self) -> Self {
        self + other
    }
    fn subtract(&self, other: &Self) -> Self {
        self - other
    }
    fn multiply(&self, other: &Self) -> Self {
        self * other
    }
    fn multiply_scalar(&self, c: &F) -> Self {
        self.scale(c)
    }
    fn multiply_monomial(&self, m: &Monomial) -> Self {
        self.mul_term(&F::one(), m)
    }
    fn subtract_multiple(&self, c: &F, m: &Monomial, other: &Self) -> Self {
        self.sub_mul(c, m, other)
    }
    fn s_polynomial(&self, other: &Polynomial<F>) -> Result<Polynomial<F>, PolynomialError> {
        if self.is_zero() || other.is_zero() {
            return Ok(Polynomial::zero(self.nvars, self.order.clone()));
        }
        let lt1 = self
            .leading_term()
            .ok_or(PolynomialError::NoLeadingMonomial)?;
        let lt2 = other
            .leading_term()
            .ok_or(PolynomialError::NoLeadingMonomial)?;
        let lcm = lt1.0.lcm(&lt2.0);
        let m1 = lcm.divide(&lt1.0).ok_or(PolynomialError::DivisionFailed)?;
        let m2 = lcm.divide(&lt2.0).ok_or(PolynomialError::DivisionFailed)?;
        let c1 = lt1.1.inverse().ok_or(PolynomialError::DivisionByZero)?;
        let c2 = lt2.1.inverse().ok_or(PolynomialError::DivisionByZero)?;
        Ok(self
            .multiply_monomial(&m1)
            .multiply_scalar(&c1)
            .subtract_multiple(&c2, &m2, other))
    }

    /// Full normal form of `self` modulo `basis`: no remaining term is divisible by any leading monomial.
    fn normal_form(&self, basis: &[Polynomial<F>]) -> Result<Polynomial<F>, PolynomialError> {
        divide_by(self, basis.iter().enumerate(), |_, _, _| {})
    }

    /// Multivariate division: quotients `q` and remainder `r` with `self = sum q[i] * divisors[i] + r`.
    fn divide_with_remainder(
        &self,
        divisors: &[Polynomial<F>],
    ) -> Result<(Vec<Polynomial<F>>, Polynomial<F>), PolynomialError> {
        let mut quotients = vec![Vec::new(); divisors.len()];
        let remainder = divide_by(self, divisors.iter().enumerate(), |i, c, m| {
            quotients[i].push(term(c.clone(), m.clone()));
        })?;
        let quotients = quotients
            .into_iter()
            .map(|terms| Polynomial {
                terms,
                nvars: self.nvars,
                order: self.order.clone(),
            })
            .collect();
        Ok((quotients, remainder))
    }
}

/// Division loop reporting each quotient term `c * m` against divisor `i`, in descending order.
pub(crate) fn divide_by<'a, F: Field + 'a>(
    polynomial: &Polynomial<F>,
    divisors: impl Iterator<Item = (usize, &'a Polynomial<F>)>,
    mut record: impl FnMut(usize, &F, &Monomial),
) -> Result<Polynomial<F>, PolynomialError> {
    let leads: Vec<_> = divisors
        .filter_map(|(i, g)| {
            g.leading_term()
                .map(|lt| (i, lt, g, divisibility_mask(&lt.0)))
        })
        .collect();
    let mut work = polynomial.clone();
    let mut head = 0;
    let mut rest = Vec::new();
    while head < work.terms.len() {
        let lt = &work.terms[head];
        let mask = divisibility_mask(&lt.0);
        let reducer = leads
            .iter()
            .filter(|(_, glt, _, divisor_mask)| divisor_mask & !mask == 0 && glt.0.divides(&lt.0))
            .min_by_key(|(_, _, g, _)| g.terms.len());
        if let Some((i, glt, g, _)) = reducer {
            let m = lt.0.divide(&glt.0).ok_or(PolynomialError::DivisionFailed)?;
            let c = glt
                .1
                .inverse()
                .map(|inverse| lt.1.clone() * inverse)
                .ok_or(PolynomialError::DivisionByZero)?;
            record(*i, &c, &m);
            work.terms.drain(..head);
            work = work.subtract_multiple(&c, &m, g);
            head = 0;
        } else {
            rest.push(lt.clone());
            head += 1;
        }
    }
    Ok(Polynomial {
        terms: rest,
        nvars: polynomial.nvars,
        order: polynomial.order.clone(),
    })
}
