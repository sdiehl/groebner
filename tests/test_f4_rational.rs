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
    assert!(is_groebner_basis(&f4).unwrap());
}

#[test]
fn f4_output_is_reduced() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y", "z"], MonomialOrder::Lex).unwrap();
    let polys = ring.parse_many("x - y^3; y^2 - z").unwrap();
    let basis = groebner_basis_f4(polys, true).unwrap();
    let text: Vec<_> = basis.iter().map(|p| ring.format(p).unwrap()).collect();
    assert_eq!(text, ["x - y*z", "y^2 - z"]);
}
