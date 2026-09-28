//! The field Q(a) of rational functions in one parameter.
//!
//! A basis over Q(a) holds for a generic value of the parameter, and [`specialize`]
//! substitutes a number.
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
use crate::polynomial::Polynomial;
use crate::polynomial::term;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use polycore::{RatFunc, Uni};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A quotient `p(a) / q(a)` of univariate rational polynomials in lowest terms with `q` monic,
/// a thin wrapper over [`polycore::RatFunc`].
///
/// Coefficient vectors are in ascending degree. [`fmt::Display`] names the parameter `a`;
/// [`crate::PolynomialRing::format`] uses the ring's parameter name instead.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RationalFunction(RatFunc<BigRational>);

impl RationalFunction {
    /// `numerator / denominator` reduced to lowest terms, or `None` if the denominator is zero.
    pub fn new(numerator: Vec<BigRational>, denominator: Vec<BigRational>) -> Option<Self> {
        let denominator = Uni::new(denominator);
        (!denominator.is_zero()).then(|| Self(RatFunc::new(&Uni::new(numerator), &denominator)))
    }

    pub fn constant(value: BigRational) -> Self {
        Self(Uni::constant(value).into())
    }

    /// The parameter `a` itself.
    pub fn parameter() -> Self {
        Self(RatFunc::var())
    }

    pub fn parameter_power(exponent: u32) -> Self {
        let mut numerator = vec![BigRational::zero(); exponent as usize];
        numerator.push(BigRational::one());
        Self(Uni::new(numerator).into())
    }

    pub fn numerator(&self) -> &[BigRational] {
        &self.0.num().0
    }

    pub fn denominator(&self) -> &[BigRational] {
        &self.0.den().0
    }

    /// Value at `a = value`, or `None` where the denominator vanishes.
    pub fn evaluate(&self, value: &BigRational) -> Option<BigRational> {
        self.0.eval(value)
    }

    /// Plain-text form with the parameter printed as `name`, parenthesized when used as a factor.
    pub(crate) fn render(&self, name: &str, as_factor: bool) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        let negative = self.numerator().last().is_some_and(Signed::is_negative);
        let numerator = if negative {
            negated(self.numerator())
        } else {
            self.numerator().to_vec()
        };
        let sign = if negative { "-" } else { "" };
        let text = poly_text(&numerator, name);
        if self.denominator().len() == 1 {
            return if terms(&numerator) > 1 && as_factor {
                format!("{sign}({text})")
            } else if negative {
                poly_text(self.numerator(), name)
            } else {
                text
            };
        }
        let numerator_text = if terms(&numerator) > 1 || !has_integer_coefficients(&numerator) {
            format!("({text})")
        } else {
            text
        };
        let denominator_text = poly_text(self.denominator(), name);
        let denominator_text = if terms(self.denominator()) > 1 {
            format!("({denominator_text})")
        } else {
            denominator_text
        };
        format!("{sign}{numerator_text}/{denominator_text}")
    }

    pub(crate) fn render_latex(&self, name: &str) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        let negative = self.numerator().last().is_some_and(Signed::is_negative);
        let numerator = if negative {
            negated(self.numerator())
        } else {
            self.numerator().to_vec()
        };
        let sign = if negative { "-" } else { "" };
        let text = poly_latex(&numerator, name);
        if self.denominator().len() > 1 {
            let denominator = poly_latex(self.denominator(), name);
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
        .map(|t| Some(term(t.1.evaluate(value)?, t.0.clone())))
        .collect::<Option<Vec<_>>>()?;
    Some(Polynomial::new(
        terms,
        polynomial.nvars,
        polynomial.order.clone(),
    ))
}

impl Zero for RationalFunction {
    fn zero() -> Self {
        Self(RatFunc::zero())
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}
impl One for RationalFunction {
    fn one() -> Self {
        Self(RatFunc::one())
    }
}
impl Neg for RationalFunction {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0)
    }
}
impl Add for RationalFunction {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self(self.0 + other.0)
    }
}
impl Sub for RationalFunction {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}
impl Mul for RationalFunction {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self(self.0 * other.0)
    }
}
impl Div for RationalFunction {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        Self(self.0 / other.0)
    }
}

impl F4Field for RationalFunction {}

impl fmt::Display for RationalFunction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.render("a", false))
    }
}

fn terms(p: &[BigRational]) -> usize {
    p.iter().filter(|c| !c.is_zero()).count()
}

fn has_integer_coefficients(p: &[BigRational]) -> bool {
    p.iter().all(num_rational::Ratio::is_integer)
}

fn negated(p: &[BigRational]) -> Vec<BigRational> {
    p.iter().map(|x| -x).collect()
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
