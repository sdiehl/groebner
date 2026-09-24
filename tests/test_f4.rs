#![allow(clippy::expect_used)]

use groebner::{
    groebner_basis, groebner_basis_f4, is_groebner_basis, MonomialOrder, Polynomial,
    PolynomialRing, PrimeField,
};

type F32003 = PrimeField<32003>;

fn ring(vars: &[&str], order: MonomialOrder) -> PolynomialRing<F32003> {
    PolynomialRing::new(vars.iter().copied(), order).expect("test ring should be valid")
}

fn parse_system(ring: &PolynomialRing<F32003>, expressions: &[&str]) -> Vec<Polynomial<F32003>> {
    expressions
        .iter()
        .map(|expression| {
            ring.parse(expression)
                .expect("test polynomial should parse")
        })
        .collect()
}

fn assert_f4_matches_buchberger(vars: &[&str], order: MonomialOrder, expressions: &[&str]) {
    let ring = ring(vars, order);
    let polynomials = parse_system(&ring, expressions);
    let buchberger =
        groebner_basis(polynomials.clone(), true).expect("Buchberger computation should succeed");
    let f4 = groebner_basis_f4(polynomials, true).expect("F4 computation should succeed");

    assert!(is_groebner_basis(&f4).expect("F4 output should be a Groebner basis"));
    assert_eq!(f4, buchberger);
}

#[test]
fn f4_matches_buchberger_for_two_variable_system() {
    assert_f4_matches_buchberger(&["x", "y"], MonomialOrder::Lex, &["x^2 - y", "x*y - 1"]);
}

#[test]
fn f4_matches_buchberger_for_katsura3() {
    assert_f4_matches_buchberger(
        &["x0", "x1", "x2"],
        MonomialOrder::Lex,
        &[
            "x0 + 2*x1 + 2*x2 - 1",
            "x0^2 + 2*x1^2 + 2*x2^2 - x0",
            "2*x0*x1 + 2*x1*x2 - x1",
        ],
    );
}

#[test]
fn f4_matches_buchberger_for_cyclic4_grlex() {
    assert_f4_matches_buchberger(
        &["x0", "x1", "x2", "x3"],
        MonomialOrder::GrLex,
        &[
            "x0 + x1 + x2 + x3",
            "x0*x1 + x1*x2 + x2*x3 + x0*x3",
            "x0*x1*x2 + x1*x2*x3 + x0*x1*x3 + x0*x2*x3",
            "x0*x1*x2*x3 - 1",
        ],
    );
}

#[test]
#[allow(deprecated)]
fn f4_rejects_mismatched_runtime_prime() {
    let ring = ring(&["x", "y"], MonomialOrder::Lex);
    let polynomials = parse_system(&ring, &["x^2 - y", "x*y - 1"]);

    assert!(groebner::groebner_basis_f4_mod(polynomials, 17, MonomialOrder::Lex).is_err());
}
