#![cfg(feature = "parallel")]

use groebner::{
    groebner_basis, groebner_basis_parallel, is_groebner_basis, is_groebner_basis_parallel,
    MonomialOrder, Polynomial, PolynomialRing,
};
use num_rational::BigRational;

fn leading_signature<F>(basis: &[Polynomial<F>]) -> Vec<Vec<u32>>
where
    F: groebner::Field,
{
    let mut signature = basis
        .iter()
        .filter_map(|polynomial| polynomial.leading_monomial())
        .map(|monomial| monomial.exponents.to_vec())
        .collect::<Vec<_>>();
    signature.sort();
    signature
}

#[test]
fn parallel_buchberger_matches_serial_leading_terms() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y", "z"], MonomialOrder::GrLex)
        .expect("ring should construct");
    let f1 = ring
        .parse("x^2 + y^2 + z^2 - 1")
        .expect("polynomial should parse");
    let f2 = ring.parse("x*y - z").expect("polynomial should parse");
    let f3 = ring.parse("y*z - x").expect("polynomial should parse");

    let serial = groebner_basis(vec![f1.clone(), f2.clone(), f3.clone()], true)
        .expect("serial basis should compute");
    let parallel =
        groebner_basis_parallel(vec![f1, f2, f3], true).expect("parallel basis should compute");

    assert_eq!(leading_signature(&parallel), leading_signature(&serial));
    assert!(is_groebner_basis(&parallel).expect("basis check should run"));
    assert!(is_groebner_basis_parallel(&parallel).expect("parallel basis check should run"));
}

#[test]
fn parallel_checker_detects_non_basis() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)
        .expect("ring should construct");
    let f1 = ring.parse("x^2 - y").expect("polynomial should parse");
    let f2 = ring.parse("x*y - 1").expect("polynomial should parse");

    assert!(!is_groebner_basis_parallel(&[f1, f2]).expect("parallel basis check should run"));
}
