use groebner::{
    groebner_basis, groebner_basis_f4, is_groebner_basis, specialize, Field, Ideal, MonomialOrder,
    ParsePolynomialError, Polynomial, PolynomialRing, RationalFunction,
};
use num_rational::BigRational;

type Result = std::result::Result<(), Box<dyn std::error::Error>>;

fn q(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

fn rf(numerator: &[i64], denominator: &[i64]) -> RationalFunction {
    RationalFunction::new(
        numerator.iter().map(|&n| q(n)).collect(),
        denominator.iter().map(|&n| q(n)).collect(),
    )
    .unwrap_or_else(RationalFunction::zero)
}

fn ring(order: MonomialOrder) -> PolynomialRing<RationalFunction> {
    PolynomialRing::with_parameter(["x", "y"], order, "a").unwrap_or_else(|e| unreachable!("{e}"))
}

fn show(
    ring: &PolynomialRing<RationalFunction>,
    basis: &[Polynomial<RationalFunction>],
) -> Vec<String> {
    basis
        .iter()
        .map(|p| ring.format(p).unwrap_or_default())
        .collect()
}

#[test]
fn arithmetic_is_in_lowest_terms() {
    let a = RationalFunction::parameter();
    let one = RationalFunction::one();
    let quotient = rf(&[-1, 0, 1], &[-1, 1]);
    assert_eq!(quotient, a.add(&one));
    assert_eq!(
        rf(&[2, 2], &[4, 4]),
        RationalFunction::constant(q(1) / q(2))
    );
    assert!(RationalFunction::new(vec![q(1)], vec![]).is_none());

    let f = rf(&[1], &[1, 1]);
    assert!(f
        .multiply(&f.inverse().unwrap_or_else(RationalFunction::zero))
        .is_one());
    assert!(f.subtract(&f).is_zero());
    assert_eq!(f.evaluate(&q(1)), Some(q(1) / q(2)));
    assert_eq!(f.evaluate(&q(-1)), None);
}

#[test]
fn display() {
    assert_eq!(rf(&[1, 1], &[1]).to_string(), "a + 1");
    assert_eq!(rf(&[1, -1], &[1]).to_string(), "-a + 1");
    assert_eq!(rf(&[0, 3], &[-2, 0, 1]).to_string(), "3*a/(a^2 - 2)");
    assert_eq!(rf(&[-1, -1], &[0, 1]).to_string(), "-(a + 1)/a");
    assert_eq!(rf(&[1], &[0, 2]).to_string(), "(1/2)/a");
}

#[test]
fn parse_and_format() -> Result {
    let ring = ring(MonomialOrder::Lex);
    let p = ring.parse("a^2*x - 1/2*a*y + 3")?;
    assert_eq!(ring.format(&p)?, "a^2*x - 1/2*a*y + 3");
    assert_eq!(ring.format(&ring.parse("a*x + x - y")?)?, "(a + 1)*x - y");
    assert_eq!(ring.format(&ring.parse("-a*x - x")?)?, "-(a + 1)*x");
    assert_eq!(
        ring.format_latex(&ring.parse("a^2*x - x + 1/3*a*y")?)?,
        "\\left(a^{2} - 1\\right) x + \\frac{1}{3} a y"
    );
    Ok(())
}

#[test]
fn parameter_errors() -> Result {
    let clash =
        PolynomialRing::<RationalFunction>::with_parameter(["x", "a"], MonomialOrder::Lex, "a");
    assert_eq!(
        clash.err(),
        Some(ParsePolynomialError::DuplicateVariable("a".into()))
    );

    let rational = PolynomialRing::<BigRational>::with_parameter(["x"], MonomialOrder::Lex, "a")?;
    assert_eq!(
        rational.parse("a*x").err(),
        Some(ParsePolynomialError::UnknownVariable("a".into()))
    );
    Ok(())
}

#[test]
fn basis_over_q_of_a() -> Result {
    let ring = ring(MonomialOrder::Lex);
    let expected = ["x - a*y", "y^2 - 1/a"];

    let f4 = groebner_basis_f4(ring.parse_many("x^2 - a; x*y - 1")?, true)?;
    assert_eq!(show(&ring, &f4), expected);
    assert!(is_groebner_basis(&f4)?);

    let buchberger = groebner_basis(ring.parse_many("x^2 - a; x*y - 1")?, true)?;
    assert_eq!(show(&ring, &buchberger), expected);
    Ok(())
}

#[test]
fn generic_basis_specializes() -> Result {
    let generic = ring(MonomialOrder::GRevLex);
    let input = "x^2 + y^2 - a; a*x*y - 1; x^3 - a*y + x";
    let basis = groebner_basis_f4(generic.parse_many(input)?, true)?;
    assert!(is_groebner_basis(&basis)?);

    let three = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::GRevLex)?;
    let direct = groebner_basis_f4(three.parse_many(&input.replace('a', "3"))?, true)?;
    let specialized: Vec<_> = basis.iter().filter_map(|p| specialize(p, &q(3))).collect();
    assert_eq!(specialized, direct);
    Ok(())
}

#[test]
fn ideal_over_q_of_a() -> Result {
    let ring = ring(MonomialOrder::GRevLex);
    let ideal = Ideal::new(ring.parse_many("x^2 - a; y^2 - a - 1")?)?;
    assert_eq!(ideal.vector_space_dimension(), Some(4));
    assert!(ideal.contains(&ring.parse("x^2*y^2 - a^2 - a")?)?);
    assert!(!ideal.contains(&ring.parse("x - a")?)?);
    Ok(())
}
