#![allow(clippy::expect_used, clippy::unwrap_used)]

use groebner::{Ideal, MonomialOrder, PolynomialRing, PrimeField, fglm, groebner_basis_f4};
use num_rational::BigRational;

type F = PrimeField<32003>;

fn ring(vars: &[&str], order: MonomialOrder) -> PolynomialRing<F> {
    PolynomialRing::new(vars.iter().copied(), order).unwrap()
}

#[test]
fn membership_and_normal_form() {
    let r = ring(&["x", "y"], MonomialOrder::GRevLex);
    let ideal = Ideal::new(r.parse_many("x^2 - y; y^2 - x").unwrap()).unwrap();
    assert!(ideal.contains(&r.parse("x^4 - x").unwrap()).unwrap());
    assert!(!ideal.contains(&r.parse("x + 1").unwrap()).unwrap());
    let nf = ideal.normal_form(&r.parse("x^3 + y^3").unwrap()).unwrap();
    assert_eq!(r.format(&nf).unwrap(), "2*x*y");
    assert!(!ideal.is_trivial());
}

#[test]
fn trivial_ideal() {
    let r = ring(&["x", "y"], MonomialOrder::Lex);
    let ideal = Ideal::new(r.parse_many("x - 1; x - 2").unwrap()).unwrap();
    assert!(ideal.is_trivial());
    assert_eq!(ideal.vector_space_dimension(), Some(0));
}

#[test]
fn dimension_and_standard_monomials() {
    let r = ring(&["x", "y", "z"], MonomialOrder::GRevLex);
    let ideal = Ideal::new(r.parse_many("x^2 - 1; y^3 - 1; z - x*y").unwrap()).unwrap();
    assert!(ideal.is_zero_dimensional());
    assert_eq!(ideal.vector_space_dimension(), Some(6));
    let positive = Ideal::new(r.parse_many("x*y - 1").unwrap()).unwrap();
    assert!(!positive.is_zero_dimensional());
    assert_eq!(positive.standard_monomials(), None);
}

#[test]
fn elimination_ideal() {
    let r = ring(&["t", "x", "y"], MonomialOrder::GRevLex);
    let ideal = Ideal::new(r.parse_many("x - t^2; y - t^3").unwrap()).unwrap();
    let eliminated = ideal.eliminate(1).unwrap();
    assert_eq!(eliminated.basis().len(), 1);
    assert_eq!(r.format(&eliminated.basis()[0]).unwrap(), "x^3 + 32002*y^2");
}

#[test]
fn radical_membership() {
    let r = ring(&["x", "y"], MonomialOrder::GRevLex);
    let ideal = Ideal::new(r.parse_many("x^2; y^3").unwrap()).unwrap();
    assert!(ideal.radical_contains(&r.parse("x*y").unwrap()).unwrap());
    assert!(!ideal.contains(&r.parse("x*y").unwrap()).unwrap());
    assert!(!ideal.radical_contains(&r.parse("x + 1").unwrap()).unwrap());
}

#[test]
fn fglm_matches_direct_lex_basis() {
    let r = ring(&["x", "y", "z"], MonomialOrder::GRevLex);
    let polys = r
        .parse_many("x^2 + y^2 + z^2 - 1; x^2 + z^2 - y; x - z")
        .unwrap();
    let grevlex = groebner_basis_f4(polys.clone(), true).unwrap();
    let via_fglm = fglm(&grevlex, &MonomialOrder::Lex).unwrap();
    let lex_ring = r.with_order(MonomialOrder::Lex);
    let direct = groebner_basis_f4(
        polys
            .iter()
            .map(|p| p.reorder(MonomialOrder::Lex))
            .collect(),
        true,
    )
    .unwrap();
    assert_eq!(via_fglm, direct);
    assert!(lex_ring.format(&via_fglm[0]).unwrap().starts_with('x'));
}

#[test]
fn fglm_rejects_positive_dimension() {
    let r = ring(&["x", "y"], MonomialOrder::GRevLex);
    let basis = groebner_basis_f4(r.parse_many("x*y - 1").unwrap(), true).unwrap();
    assert!(fglm(&basis, &MonomialOrder::Lex).is_err());
}

#[test]
fn change_order_over_rationals() {
    let r = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::GRevLex).unwrap();
    let ideal = Ideal::new(r.parse_many("x^2 + y^2 - 1; x - y^3").unwrap()).unwrap();
    let lex = ideal.change_order(MonomialOrder::Lex).unwrap();
    let univariate = lex
        .basis()
        .iter()
        .find(|p| p.terms.iter().all(|t| t.monomial.exponents()[0] == 0))
        .expect("lex basis has a univariate polynomial in y");
    assert_eq!(r.format(univariate).unwrap(), "y^6 + y^2 - 1");
}

#[test]
fn multiplication_matrix_shape() {
    let r = ring(&["x", "y"], MonomialOrder::GRevLex);
    let ideal = Ideal::new(r.parse_many("x^2 - y; y^2 - 1").unwrap()).unwrap();
    let m = ideal.multiplication_matrix(0).unwrap().unwrap();
    assert_eq!(m.len(), 4);
    assert!(m.iter().all(|column| column.len() == 4));
}
