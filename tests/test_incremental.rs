#![allow(clippy::expect_used)]

use groebner::{
    groebner_basis, groebner_basis_incremental, MonomialOrder, Polynomial, PolynomialRing,
};
use num_rational::BigRational;

fn basis_signature(polynomials: &[Polynomial<BigRational>]) -> Vec<String> {
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)
        .expect("ring should be valid");
    polynomials
        .iter()
        .map(|polynomial| ring.format(polynomial).expect("should format"))
        .collect()
}

#[test]
fn incremental_basis_matches_full_recompute() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)
        .expect("ring should be valid");
    let f1 = ring.parse("x^2 - y").expect("f1 should parse");
    let f2 = ring.parse("x*y - 1").expect("f2 should parse");
    let f3 = ring.parse("y^2 - x").expect("f3 should parse");

    let initial = groebner_basis(vec![f1.clone(), f2.clone()], ring.order(), true)
        .expect("initial basis should compute");
    let incremental = groebner_basis_incremental(initial, vec![f3.clone()], ring.order(), true)
        .expect("incremental basis should compute");
    let full =
        groebner_basis(vec![f1, f2, f3], ring.order(), true).expect("full basis should compute");

    assert_eq!(basis_signature(&incremental), basis_signature(&full));
}
