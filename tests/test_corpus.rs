//! Fast tier of the public benchmark corpus; the full corpus runs through
//! `cargo run --release --example corpus`.

mod common;

use common::corpus::{self, Algorithm, Field};
use std::path::Path;

fn smoke(field: Field, algorithm: Algorithm) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus");
    let Ok(list) = std::fs::read_to_string(root.join("smoke.txt")) else {
        return;
    };
    let mut failures = Vec::new();
    for name in list
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let system = match corpus::load(&root.join(format!("{name}.txt"))) {
            Ok(system) => system,
            Err(e) => {
                failures.push(format!("{name}: {e}"));
                continue;
            }
        };
        let (Some(reference), true) = (&system.reference, field.applies(&system)) else {
            continue;
        };
        let result = corpus::run(&system, field, algorithm)
            .and_then(|rows| corpus::compare(&rows, &reference.rows));
        if let Err(e) = result {
            failures.push(format!("{name}: {e}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{field:?} {algorithm:?} failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn corpus_f4_zp() {
    smoke(Field::Zp, Algorithm::F4);
}

#[test]
fn corpus_f4_gf32003() {
    smoke(Field::Gf32003, Algorithm::F4);
}

#[test]
fn corpus_f4_rational() {
    smoke(Field::Qq, Algorithm::F4);
}

#[test]
fn corpus_buchberger_zp() {
    smoke(Field::Zp, Algorithm::Buchberger);
}

#[test]
fn fingerprint_omits_terms_vanishing_at_reference_prime() {
    let ring = groebner::PolynomialRing::<num_rational::BigRational>::new(
        ["x", "y"],
        groebner::MonomialOrder::Lex,
    )
    .unwrap();
    let rational = ring.parse_many("x + 32003*y + 1").unwrap();
    let modular = rational
        .iter()
        .map(|p| {
            p.try_map(|c| groebner::Fp::from_rational(c, 32003))
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(
        corpus::fingerprint(&rational, |c| polycore::crt::reduce(c, 32003)),
        corpus::fingerprint(&modular, |c| Some(c.value()))
    );
}
