#![cfg(feature = "parameters")]

use groebner::{Frac, MonomialOrder, Polynomial, PolynomialRing, groebner_basis_f4};
use num_rational::BigRational;

fn q(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

/// The basis over `Q(a, b)` specializes to the basis of the specialized system over Q.
#[test]
fn specializes() -> Result<(), Box<dyn std::error::Error>> {
    let system = "x^2 + a*y - b; x*y - a*b*z; y^2 - z + b; x + y + z - a";
    let at = [q(2), q(-3)];
    let ring = PolynomialRing::<Frac>::with_parameters(
        ["x", "y", "z"],
        MonomialOrder::GRevLex,
        ["a", "b"],
    )?;
    let generic = groebner_basis_f4(ring.parse_many(system)?, true)?;
    let special: Vec<Polynomial<BigRational>> = generic
        .iter()
        .map(|p| p.try_map(|c| c.eval(&at)).ok_or("pole"))
        .collect::<Result<_, _>>()?;
    let rational = PolynomialRing::<BigRational>::new(["x", "y", "z"], MonomialOrder::GRevLex)?;
    let direct = system.replace('a', "(2)").replace('b', "(-3)");
    assert_eq!(
        special,
        groebner_basis_f4(rational.parse_many(&direct)?, true)?
    );
    Ok(())
}

#[test]
fn names() -> Result<(), Box<dyn std::error::Error>> {
    let ring = PolynomialRing::<Frac>::with_parameters(["x"], MonomialOrder::Lex, ["s", "t"])?;
    let p = ring.parse("(s^2 - t^2)/(s - t)*x - 2/t")?;
    assert_eq!(ring.format(&p)?, "(s + t)*x - 2/t");
    assert_eq!(
        ring.format_latex(&p)?,
        "\\left(s + t\\right) x - \\frac{2}{t}"
    );
    assert!(
        PolynomialRing::<Frac>::with_parameters(["x"], MonomialOrder::Lex, ["s", "s"]).is_err()
    );
    Ok(())
}
