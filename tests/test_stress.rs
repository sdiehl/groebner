extern crate groebner;
use groebner::{groebner_basis, Monomial, MonomialOrder, Polynomial, Term};
use num_rational::BigRational;
use std::time::Instant;

const SEED: usize = 7919;

/// Helper to create a more varied polynomial with nvars variables, degree up to max_deg, and nterms terms
fn make_poly(
    nvars: usize,
    max_deg: u32,
    nterms: usize,
    order: MonomialOrder,
    offset: usize,
) -> Polynomial<BigRational> {
    let mut terms = Vec::new();
    for i in 0..nterms {
        let mut exps = vec![0u32; nvars];
        // Vary exponents: each term has a different pattern
        for (v, exp) in exps.iter_mut().enumerate() {
            *exp = (((i + v * offset + SEED * v) * (v + 1)) % (max_deg as usize + 1)) as u32;
        }
        // Vary coefficients: alternate sign, use offset
        let sign = if (i + offset) % 2 == 0 { 1 } else { -1 };
        let coeff = BigRational::new(
            (sign * ((i + 1 + offset) as i32)).into(),
            (1 + (i % 3) as i32).into(),
        );
        terms.push(Term::new(coeff, Monomial::new(exps)));
    }
    Polynomial::new(terms, nvars, order)
}

#[test]
fn stress_small() {
    // 2 polynomials, 2 variables, degree up to 2, 3 terms each
    let nvars = 2;
    let nterms = 3;
    let polys: Vec<_> = (0..2)
        .map(|i| make_poly(nvars, 2, nterms, MonomialOrder::Lex, i))
        .collect();
    let start = Instant::now();
    let basis = groebner_basis(polys, MonomialOrder::Lex, true).expect("Groebner failed");
    let elapsed = start.elapsed();
    println!(
        "[stress_small] Time: {}.{:03} seconds, basis size: {}",
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        basis.len()
    );
}

#[test]
fn stress_medium() {
    // 3 polynomials, 3 variables, degree up to 3, 4 terms each
    let nvars = 3;
    let nterms = 4;
    let polys: Vec<_> = (0..3)
        .map(|i| make_poly(nvars, 3, nterms, MonomialOrder::GrLex, i))
        .collect();
    let start = Instant::now();
    let basis = groebner_basis(polys, MonomialOrder::GrLex, true).expect("Groebner failed");
    let elapsed = start.elapsed();
    println!(
        "[stress_medium] Time: {}.{:03} seconds, basis size: {}",
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        basis.len()
    );
}

#[test]
fn stress_large() {
    // 4 polynomials, 4 variables, degree up to 4, 5 terms each
    let nvars = 4;
    let nterms = 5;
    let polys: Vec<_> = (0..4)
        .map(|i| make_poly(nvars, 3, nterms, MonomialOrder::Lex, i))
        .collect();
    let start = Instant::now();
    let basis = groebner_basis(polys, MonomialOrder::Lex, true).expect("Groebner failed");
    let elapsed = start.elapsed();
    println!(
        "[stress_large] Time: {}.{:03} seconds, basis size: {}",
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        basis.len()
    );
}

#[test]
fn stress_xlarge() {
    // 5 polynomials, 5 variables, degree up to 5, 8 terms each
    let nvars = 5;
    let nterms = 8;
    let polys: Vec<_> = (0..5)
        .map(|i| make_poly(nvars, 5, nterms, MonomialOrder::GrLex, i))
        .collect();
    let start = Instant::now();
    let basis = groebner_basis(polys, MonomialOrder::GrLex, true).expect("Groebner failed");
    let elapsed = start.elapsed();
    println!(
        "[stress_xlarge] Time: {}.{:03} seconds, basis size: {}",
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        basis.len()
    );
}
