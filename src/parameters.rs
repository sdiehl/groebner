//! Coefficients in `Q(a, b, ...)`, several parameters, as [`zippel_gcd::Frac`].
//!
//! Coefficients are cancelled by `zippel-gcd`, and the ring is built with
//! [`crate::PolynomialRing::with_parameters`].
//!
//! ```
//! use groebner::{Frac, MonomialOrder, PolynomialRing, groebner_basis_f4};
//!
//! let ring = PolynomialRing::<Frac>::with_parameters(["x", "y"], MonomialOrder::Lex, ["a", "b"])?;
//! let basis = groebner_basis_f4(ring.parse_many("x^2 - a; x*y - b")?, true)?;
//! let printed: Vec<String> = basis.iter().map(|p| ring.format(p)).collect::<Result<_, _>>()?;
//! assert_eq!(printed, ["x - a/b*y", "y^2 - b^2/a"]);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use num_rational::BigRational;
use num_traits::{One, Signed};
use polycore::{Order, Ring};
pub use zippel_gcd::Frac;

use crate::f4::F4Field;
use crate::ring::{ParseCoefficient, ParsePolynomialError};

type Poly = polycore::Poly<BigRational>;

impl F4Field for Frac {}

impl ParseCoefficient for Frac {
    fn parameter(index: usize, count: usize) -> Option<Self> {
        Some(Poly::var(index, count, Order::GRevLex).into())
    }

    fn format_coefficient(&self, parameter: &str) -> String {
        self.format_coefficient_in(&[parameter.to_owned()])
    }

    fn format_coefficient_in(&self, parameters: &[String]) -> String {
        let ring = names(self, parameters);
        let (sign, num) = unsigned(self.num());
        let text = ring.show(&num);
        if self.den().is_constant() {
            return if num.terms.len() > 1 {
                format!("{sign}({text})")
            } else {
                format!("{sign}{text}")
            };
        }
        let group = |p: &Poly, t: String| {
            let simple = p.terms.len() == 1 && p.terms[0].1.is_integer();
            if simple { t } else { format!("({t})") }
        };
        let den = ring.show(self.den());
        format!("{sign}{}/{}", group(&num, text), group(self.den(), den))
    }

    fn format_coefficient_latex_in(&self, parameters: &[String]) -> String {
        let ring = names(self, parameters);
        let (sign, num) = unsigned(self.num());
        let text = ring.latex(&num);
        if self.den().is_constant() {
            return if num.terms.len() > 1 {
                format!("{sign}\\left({text}\\right)")
            } else {
                format!("{sign}{text}")
            };
        }
        format!("{sign}\\frac{{{text}}}{{{}}}", ring.latex(self.den()))
    }

    fn parse_coefficient(
        numerator: &str,
        denominator: Option<&str>,
        modulus: Option<u64>,
    ) -> Result<Self, ParsePolynomialError> {
        BigRational::parse_coefficient(numerator, denominator, modulus).map(Self::from)
    }
}

/// A ring naming the parameters of `f`, padded if `f` has more variables than names.
fn names(f: &Frac, parameters: &[String]) -> Ring {
    let n = f.num().nvars.max(f.den().nvars);
    let mut names = parameters.to_vec();
    names.extend((names.len()..n).map(|i| format!("a{i}")));
    names.truncate(n.max(parameters.len()));
    Ring::new(names, f.num().order.clone())
}

/// The sign of the leading coefficient and `p` made positive there.
fn unsigned(p: &Poly) -> (&'static str, Poly) {
    if p.lc().is_some_and(Signed::is_negative) {
        ("-", p.scale(&-BigRational::one()))
    } else {
        ("", p.clone())
    }
}
