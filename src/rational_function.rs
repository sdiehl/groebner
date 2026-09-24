//! The field Q(a) of rational functions in one parameter.
//!
//! ```
//! use groebner::{groebner_basis_f4, MonomialOrder, PolynomialRing, RationalFunction};
//!
//! let ring = PolynomialRing::<RationalFunction>::with_parameter(["x", "y"], MonomialOrder::Lex, "a")?;
//! let basis = groebner_basis_f4(ring.parse_many("x^2 - a; x*y - 1")?, true)?;
//! let printed: Vec<String> = basis.iter().map(|p| ring.format(p)).collect::<Result<_, _>>()?;
//! assert_eq!(printed, ["x - a*y", "y^2 - 1/a"]);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::f4::F4Field;
use crate::field::Field;
use crate::polynomial::{Polynomial, Term};
use num_rational::BigRational;
use num_traits::Signed;
use std::fmt;

type Poly = Vec<BigRational>;

/// A quotient `p(a) / q(a)` of univariate rational polynomials in lowest terms with `q` monic.
///
/// Coefficient vectors are in ascending degree. [`fmt::Display`] names the parameter `a`;
/// [`crate::PolynomialRing::format`] uses the ring's parameter name instead.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RationalFunction {
    numerator: Poly,
    denominator: Poly,
}

impl RationalFunction {
    /// `numerator / denominator` reduced to lowest terms, or `None` if the denominator is zero.
    pub fn new(numerator: Vec<BigRational>, denominator: Vec<BigRational>) -> Option<Self> {
        let numerator = trimmed(numerator);
        let denominator = trimmed(denominator);
        if denominator.is_empty() {
            return None;
        }
        if numerator.is_empty() {
            return Some(Self::zero());
        }
        let g = gcd(&numerator, &denominator);
        let numerator = divide_exact(&numerator, &g);
        let denominator = divide_exact(&denominator, &g);
        let lead = denominator.last()?.clone();
        Some(Self {
            numerator: scale(&numerator, &lead.recip()),
            denominator: scale(&denominator, &lead.recip()),
        })
    }

    pub fn constant(value: BigRational) -> Self {
        Self {
            numerator: trimmed(vec![value]),
            denominator: vec![BigRational::one()],
        }
    }

    /// The parameter `a` itself.
    pub fn parameter() -> Self {
        Self::parameter_power(1)
    }

    pub fn parameter_power(exponent: u32) -> Self {
        let mut numerator = vec![BigRational::zero(); exponent as usize];
        numerator.push(BigRational::one());
        Self {
            numerator,
            denominator: vec![BigRational::one()],
        }
    }

    pub fn numerator(&self) -> &[BigRational] {
        &self.numerator
    }

    pub fn denominator(&self) -> &[BigRational] {
        &self.denominator
    }

    /// Value at `a = value`, or `None` where the denominator vanishes.
    pub fn evaluate(&self, value: &BigRational) -> Option<BigRational> {
        let denominator = horner(&self.denominator, value);
        (!denominator.is_zero()).then(|| horner(&self.numerator, value) / denominator)
    }

    /// Plain-text form with the parameter printed as `name`, parenthesized when used as a factor.
    pub(crate) fn render(&self, name: &str, as_factor: bool) -> String {
        if self.numerator.is_empty() {
            return "0".to_string();
        }
        let negative = self.numerator.last().is_some_and(Signed::is_negative);
        let numerator = if negative {
            negated(&self.numerator)
        } else {
            self.numerator.clone()
        };
        let sign = if negative { "-" } else { "" };
        let text = poly_text(&numerator, name);
        if self.denominator.len() == 1 {
            return if terms(&numerator) > 1 && as_factor {
                format!("{sign}({text})")
            } else if negative {
                poly_text(&self.numerator, name)
            } else {
                text
            };
        }
        let numerator_text = if terms(&numerator) > 1 || !has_integer_coefficients(&numerator) {
            format!("({text})")
        } else {
            text
        };
        let denominator_text = poly_text(&self.denominator, name);
        let denominator_text = if terms(&self.denominator) > 1 {
            format!("({denominator_text})")
        } else {
            denominator_text
        };
        format!("{sign}{numerator_text}/{denominator_text}")
    }

    pub(crate) fn render_latex(&self, name: &str) -> String {
        if self.numerator.is_empty() {
            return "0".to_string();
        }
        let negative = self.numerator.last().is_some_and(Signed::is_negative);
        let numerator = if negative {
            negated(&self.numerator)
        } else {
            self.numerator.clone()
        };
        let sign = if negative { "-" } else { "" };
        let text = poly_latex(&numerator, name);
        if self.denominator.len() > 1 {
            let denominator = poly_latex(&self.denominator, name);
            format!("{sign}\\frac{{{text}}}{{{denominator}}}")
        } else if terms(&numerator) > 1 {
            format!("{sign}\\left({text}\\right)")
        } else {
            format!("{sign}{text}")
        }
    }
}

/// Substitute `a = value` into every coefficient, or `None` if some denominator vanishes there.
pub fn specialize(
    polynomial: &Polynomial<RationalFunction>,
    value: &BigRational,
) -> Option<Polynomial<BigRational>> {
    let terms = polynomial
        .terms
        .iter()
        .map(|t| {
            Some(Term::new(
                t.coefficient.evaluate(value)?,
                t.monomial.clone(),
            ))
        })
        .collect::<Option<Vec<_>>>()?;
    Some(Polynomial::new(
        terms,
        polynomial.nvars,
        polynomial.order.clone(),
    ))
}

impl Field for RationalFunction {
    fn zero() -> Self {
        Self {
            numerator: Vec::new(),
            denominator: vec![BigRational::one()],
        }
    }
    fn one() -> Self {
        Self::constant(BigRational::one())
    }
    fn is_zero(&self) -> bool {
        self.numerator.is_empty()
    }
    fn is_one(&self) -> bool {
        self.denominator.len() == 1 && self.numerator.len() == 1 && self.numerator[0].is_one()
    }
    fn add(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        let (numerator, denominator) = if self.denominator == other.denominator {
            (
                sum(&self.numerator, &other.numerator),
                self.denominator.clone(),
            )
        } else {
            (
                sum(
                    &product(&self.numerator, &other.denominator),
                    &product(&other.numerator, &self.denominator),
                ),
                product(&self.denominator, &other.denominator),
            )
        };
        Self::new(numerator, denominator).unwrap_or_else(Self::zero)
    }
    fn subtract(&self, other: &Self) -> Self {
        self.add(&other.negate())
    }
    fn multiply(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        Self::new(
            product(&self.numerator, &other.numerator),
            product(&self.denominator, &other.denominator),
        )
        .unwrap_or_else(Self::zero)
    }
    fn negate(&self) -> Self {
        Self {
            numerator: negated(&self.numerator),
            denominator: self.denominator.clone(),
        }
    }
    fn inverse(&self) -> Option<Self> {
        let lead = self.numerator.last()?.recip();
        Some(Self {
            numerator: scale(&self.denominator, &lead),
            denominator: scale(&self.numerator, &lead),
        })
    }
}

impl F4Field for RationalFunction {}

impl fmt::Display for RationalFunction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.render("a", false))
    }
}

fn trimmed(mut p: Poly) -> Poly {
    while p.last().is_some_and(Field::is_zero) {
        p.pop();
    }
    p
}

fn terms(p: &[BigRational]) -> usize {
    p.iter().filter(|c| !c.is_zero()).count()
}

fn has_integer_coefficients(p: &[BigRational]) -> bool {
    p.iter().all(num_rational::Ratio::is_integer)
}

fn scale(p: &[BigRational], c: &BigRational) -> Poly {
    p.iter().map(|x| x * c).collect()
}

fn negated(p: &[BigRational]) -> Poly {
    p.iter().map(|x| -x).collect()
}

fn sum(a: &[BigRational], b: &[BigRational]) -> Poly {
    let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut out = long.to_vec();
    for (x, y) in out.iter_mut().zip(short) {
        *x += y;
    }
    trimmed(out)
}

fn product(a: &[BigRational], b: &[BigRational]) -> Poly {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![BigRational::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    trimmed(out)
}

/// Quotient and remainder of `a` by a nonzero `b`.
fn divide(a: &[BigRational], b: &[BigRational]) -> (Poly, Poly) {
    let mut remainder = a.to_vec();
    let Some(lead) = b.last().map(BigRational::recip) else {
        return (Vec::new(), remainder);
    };
    if a.len() < b.len() {
        return (Vec::new(), remainder);
    }
    let mut quotient = vec![BigRational::zero(); a.len() - b.len() + 1];
    for shift in (0..quotient.len()).rev() {
        let c = &remainder[shift + b.len() - 1] * &lead;
        if c.is_zero() {
            continue;
        }
        for (k, y) in b.iter().enumerate() {
            remainder[shift + k] -= &c * y;
        }
        quotient[shift] = c;
    }
    remainder.truncate(b.len() - 1);
    (trimmed(quotient), trimmed(remainder))
}

fn divide_exact(a: &[BigRational], b: &[BigRational]) -> Poly {
    divide(a, b).0
}

fn gcd(a: &[BigRational], b: &[BigRational]) -> Poly {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    while !b.is_empty() {
        let r = divide(&a, &b).1;
        a = b;
        b = r;
    }
    match a.last().map(BigRational::recip) {
        Some(lead) => scale(&a, &lead),
        None => vec![BigRational::one()],
    }
}

fn horner(p: &[BigRational], value: &BigRational) -> BigRational {
    p.iter()
        .rev()
        .fold(BigRational::zero(), |acc, c| acc * value + c)
}

fn poly_text(p: &[BigRational], name: &str) -> String {
    let mut out = String::new();
    for (k, c) in p.iter().enumerate().rev().filter(|(_, c)| !c.is_zero()) {
        let magnitude = c.abs();
        let power = match k {
            0 => String::new(),
            1 => name.to_string(),
            _ => format!("{name}^{k}"),
        };
        let body = match (k, magnitude.is_one()) {
            (0, _) => magnitude.to_string(),
            (_, true) => power,
            _ => format!("{magnitude}*{power}"),
        };
        push_signed(&mut out, c.is_negative(), &body);
    }
    out
}

fn poly_latex(p: &[BigRational], name: &str) -> String {
    let mut out = String::new();
    for (k, c) in p.iter().enumerate().rev().filter(|(_, c)| !c.is_zero()) {
        let magnitude = c.abs();
        let coefficient = if magnitude.is_integer() {
            magnitude.to_string()
        } else {
            format!("\\frac{{{}}}{{{}}}", magnitude.numer(), magnitude.denom())
        };
        let power = match k {
            0 => String::new(),
            1 => name.to_string(),
            _ => format!("{name}^{{{k}}}"),
        };
        let body = match (k, magnitude.is_one()) {
            (0, _) => coefficient,
            (_, true) => power,
            _ => format!("{coefficient} {power}"),
        };
        push_signed(&mut out, c.is_negative(), &body);
    }
    out
}

fn push_signed(out: &mut String, negative: bool, body: &str) {
    match (out.is_empty(), negative) {
        (true, true) => out.push('-'),
        (true, false) => {}
        (false, true) => out.push_str(" - "),
        (false, false) => out.push_str(" + "),
    }
    out.push_str(body);
}
