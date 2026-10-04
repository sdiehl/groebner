#![allow(clippy::expect_used, clippy::unwrap_used)]

use groebner::{
    MonomialOrder, PolynomialRing, groebner_basis, groebner_basis_f4, is_groebner_basis,
};
use num_rational::BigRational;

#[test]
fn f4_over_rationals_matches_buchberger() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y", "z"], MonomialOrder::GRevLex).unwrap();
    let polys = ring
        .parse_many("x^2 + y^2 + z^2 - 1; x^2 + z^2 - y; x - z")
        .unwrap();
    let f4 = groebner_basis_f4(polys.clone(), true).unwrap();
    let buchberger = groebner_basis(polys, true).unwrap();
    assert_eq!(f4, buchberger);
    assert!(is_groebner_basis(&f4));
}

#[test]
fn f4_output_is_reduced() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y", "z"], MonomialOrder::Lex).unwrap();
    let polys = ring.parse_many("x - y^3; y^2 - z").unwrap();
    let basis = groebner_basis_f4(polys, true).unwrap();
    let text: Vec<_> = basis.iter().map(|p| ring.format(p).unwrap()).collect();
    assert_eq!(text, ["x - y*z", "y^2 - z"]);
}

#[test]
fn modular_matches_direct_across_orders_and_large_coefficients() {
    use groebner::f4::groebner_basis_f4_direct;
    for order in [
        MonomialOrder::Lex,
        MonomialOrder::GrLex,
        MonomialOrder::GRevLex,
        MonomialOrder::weighted(vec![2, 1], MonomialOrder::Lex),
        MonomialOrder::elimination(1, 1),
    ] {
        let ring = PolynomialRing::<BigRational>::new(["x", "y"], order).unwrap();
        let polys = ring
            .parse_many("(x + y)^2 - 1/7; 123456789012345678901234567890123456789*x - 19*y")
            .unwrap();
        let modular = groebner::groebner_basis_f4_rational(polys.clone(), true).unwrap();
        assert_eq!(modular, groebner_basis_f4_direct(polys, true).unwrap());
    }
}

#[test]
fn modular_validates_input_and_handles_the_unit_ideal() {
    use groebner::{GroebnerError, Polynomial};
    assert!(matches!(
        groebner_basis_f4::<BigRational>(vec![], true),
        Err(GroebnerError::EmptyInput)
    ));
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex).unwrap();
    let zero = Polynomial::<BigRational>::zero(2, MonomialOrder::Lex);
    assert!(matches!(
        groebner_basis_f4(vec![zero], true),
        Err(GroebnerError::EmptyInput)
    ));
    let mixed = vec![
        ring.parse("x").unwrap(),
        ring.with_order(MonomialOrder::GRevLex).parse("y").unwrap(),
    ];
    assert!(matches!(
        groebner_basis_f4(mixed, true),
        Err(GroebnerError::OrderMismatch)
    ));
    assert_eq!(
        groebner_basis_f4(ring.parse_many("x; x - 1").unwrap(), false).unwrap(),
        ring.parse_many("1").unwrap()
    );
}
