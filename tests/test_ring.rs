use groebner::{MonomialOrder, PolynomialRing};
use num_rational::BigRational;

#[test]
fn parses_polynomial_strings() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)
        .expect("ring should be valid");
    let polynomial = ring.parse("x^2 - y + 3/2*x*y").expect("should parse");

    assert_eq!(polynomial.nvars, 2);
    assert_eq!(polynomial.terms.len(), 3);
    assert_eq!(
        ring.format(&polynomial).expect("should format"),
        "x^2 + 3/2*x*y - y"
    );
}

#[test]
fn parses_implicit_multiplication_and_grouped_numbers() {
    let ring = PolynomialRing::<BigRational>::new(["d", "h", "w0", "z0"], MonomialOrder::Lex)
        .expect("ring should be valid");
    let polynomial = ring
        .parse("190 976 d h^2 w0 z0^3 - 4*d*h*w0*z0")
        .expect("should parse");

    assert_eq!(polynomial.terms.len(), 2);
    assert_eq!(
        ring.format(&polynomial).expect("should format"),
        "190976*d*h^2*w0*z0^3 - 4*d*h*w0*z0"
    );
}

#[test]
fn variable_order_controls_lexicographic_order() {
    let z1_first = PolynomialRing::<BigRational>::new(["z1", "z2", "z3"], MonomialOrder::Lex)
        .expect("ring should be valid");
    let z3_first = PolynomialRing::<BigRational>::new(["z3", "z1", "z2"], MonomialOrder::Lex)
        .expect("ring should be valid");

    let expression = "z1 + z3";
    assert_eq!(
        z1_first
            .format(&z1_first.parse(expression).expect("should parse"))
            .expect("should format"),
        "z1 + z3"
    );
    assert_eq!(
        z3_first
            .format(&z3_first.parse(expression).expect("should parse"))
            .expect("should format"),
        "z3 + z1"
    );
}

#[test]
fn rejects_unknown_variables() {
    let ring = PolynomialRing::<BigRational>::new(["z0", "z2"], MonomialOrder::Lex)
        .expect("ring should be valid");

    assert!(ring.parse("z011 + z2").is_err());
}

#[test]
fn parse_many_accepts_wrapped_polynomials() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)
        .expect("ring should be valid");
    let polynomials = ring
        .parse_many("x^2\n - y, x*y - 1")
        .expect("system should parse");

    assert_eq!(polynomials.len(), 2);
    assert_eq!(
        ring.format(&polynomials[0]).expect("should format"),
        "x^2 - y"
    );
    assert_eq!(
        ring.format(&polynomials[1]).expect("should format"),
        "x*y - 1"
    );
}
