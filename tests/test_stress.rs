extern crate groebner;
use groebner::{groebner_basis, MonomialOrder, Polynomial, PolynomialRing};
use num_rational::BigRational;
use std::time::Instant;

const SEED: usize = 7919;
const TEST_VARIABLES: [&str; 5] = ["x", "y", "z", "w", "u"];

/// Helper to create a more varied polynomial with nvars variables, degree up to max_deg, and nterms terms
#[allow(clippy::expect_used)]
fn make_poly(
    nvars: usize,
    max_deg: u32,
    nterms: usize,
    order: MonomialOrder,
    offset: usize,
) -> Polynomial<BigRational> {
    let ring = PolynomialRing::<BigRational>::new(TEST_VARIABLES[..nvars].iter().copied(), order)
        .expect("test ring should be valid");
    let mut expression = String::new();
    for i in 0..nterms {
        let mut exps = vec![0u32; nvars];
        // Vary exponents: each term has a different pattern
        for (v, exp) in exps.iter_mut().enumerate() {
            *exp = (((i + v * offset + SEED * v) * (v + 1)) % (max_deg as usize + 1)) as u32;
        }
        // Vary coefficients: alternate sign, use offset
        let sign = if (i + offset) % 2 == 0 { 1 } else { -1 };
        push_term(
            &mut expression,
            i,
            sign * ((i + 1 + offset) as i32),
            1 + (i % 3) as i32,
            &exps,
        );
    }
    ring.parse(&expression)
        .expect("generated stress polynomial should parse")
}

fn push_term(
    expression: &mut String,
    index: usize,
    numerator: i32,
    denominator: i32,
    exponents: &[u32],
) {
    let negative = numerator < 0;
    let abs_numerator = numerator.abs();
    if index == 0 {
        if negative {
            expression.push('-');
        }
    } else if negative {
        expression.push_str(" - ");
    } else {
        expression.push_str(" + ");
    }

    let is_constant = exponents.iter().all(|&exponent| exponent == 0);
    let mut factors = Vec::new();
    if abs_numerator != 1 || denominator != 1 || is_constant {
        if denominator == 1 {
            factors.push(abs_numerator.to_string());
        } else {
            factors.push(format!("{abs_numerator}/{denominator}"));
        }
    }
    for (variable, exponent) in TEST_VARIABLES.iter().zip(exponents) {
        match exponent {
            0 => {}
            1 => factors.push((*variable).to_string()),
            exponent => factors.push(format!("{variable}^{exponent}")),
        }
    }
    expression.push_str(&factors.join("*"));
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
    let basis = groebner_basis(polys, true).expect("Groebner failed");
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
    let basis = groebner_basis(polys, true).expect("Groebner failed");
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
    let basis = groebner_basis(polys, true).expect("Groebner failed");
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
    let basis = groebner_basis(polys, true).expect("Groebner failed");
    let elapsed = start.elapsed();
    println!(
        "[stress_xlarge] Time: {}.{:03} seconds, basis size: {}",
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        basis.len()
    );
}
