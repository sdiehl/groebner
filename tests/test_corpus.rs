//! Fast tier of the public benchmark corpus; the full corpus runs through
//! `cargo run --release --example corpus`.

mod common;

use common::corpus::{self, Algorithm, Field};
use std::path::Path;

fn smoke(field: Field, algorithm: Algorithm) -> Result<(), String> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus");
    // The published crate intentionally omits corpus data, but a checkout must
    // not silently skip smoke coverage if its manifest or a fixture disappears.
    if !root.exists() {
        return Ok(());
    }
    let systems = corpus::load_manifest(&root, &root.join("smoke.txt"))?;
    let mut failures = Vec::new();
    for system in systems {
        let name = &system.name;
        let (Some(reference), true) = (&system.reference, field.applies(&system)) else {
            continue;
        };
        let result = corpus::run(&system, field, algorithm)
            .and_then(|rows| corpus::compare(&rows, &reference.rows));
        if let Err(e) = result {
            failures.push(format!("{name}: {e}"));
        }
    }
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "{field:?} {algorithm:?} failures:\n{}",
            failures.join("\n")
        ))
    }
}

#[test]
fn corpus_f4_zp() -> Result<(), String> {
    smoke(Field::Zp, Algorithm::F4)
}

#[test]
fn corpus_f4_gf32003() -> Result<(), String> {
    smoke(Field::Gf32003, Algorithm::F4)
}

#[test]
fn corpus_f4_rational() -> Result<(), String> {
    smoke(Field::Qq, Algorithm::F4)
}

#[test]
fn corpus_buchberger_zp() -> Result<(), String> {
    smoke(Field::Zp, Algorithm::Buchberger)
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

#[test]
fn corpus_manifest_rejects_empty_duplicate_and_invalid_names() {
    assert!(corpus::manifest_names("# no cases\n").is_err());
    assert!(corpus::manifest_names("classic/cyclic3\nclassic/cyclic3").is_err());
    for invalid in ["/absolute", "../outside", "classic/../cyclic3", "classic/"] {
        assert!(corpus::manifest_names(invalid).is_err(), "{invalid}");
    }
    assert_eq!(
        corpus::manifest_names("# tier\nclassic/cyclic3 # small\n\nclassic/katsura3\n").unwrap(),
        ["classic/cyclic3", "classic/katsura3"]
    );
}

#[test]
fn ci_manifest_has_complete_inputs_and_references() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/corpus");
    // Corpus data is intentionally excluded from the crates.io package.
    if !root.exists() {
        return;
    }
    let systems = corpus::load_manifest(&root, &root.join("ci.txt")).unwrap();
    assert!(
        systems.len() >= 200,
        "CI must retain substantial corpus coverage"
    );
}
