#![allow(clippy::expect_used, clippy::unwrap_used)]

use groebner::{
    F4Field, Field, Ideal, LiftBasis, MonomialOrder, ParseCoefficient, Polynomial, PolynomialRing,
    PrimeField, groebner_basis, verify_lift,
};
use num_rational::BigRational;

type Zp = PrimeField<32003>;

const CYCLIC3: (&[&str], &str) = (&["x", "y", "z"], "x + y + z; x*y + y*z + z*x; x*y*z - 1");
const CYCLIC4: (&[&str], &str) = (
    &["a", "b", "c", "d"],
    "a + b + c + d; a*b + b*c + c*d + d*a; a*b*c + b*c*d + c*d*a + d*a*b; a*b*c*d - 1",
);
const KATSURA3: (&[&str], &str) = (
    &["x0", "x1", "x2", "x3"],
    "x0 + 2*x1 + 2*x2 + 2*x3 - 1; x0^2 + 2*x1^2 + 2*x2^2 + 2*x3^2 - x0; \
     2*x0*x1 + 2*x1*x2 + 2*x2*x3 - x1; x1^2 + 2*x0*x2 + 2*x1*x3 - x2",
);
// Redundant, non-monic generators and a zero: the basis collapses and interreduction does work.
const REDUNDANT: (&[&str], &str) = (
    &["x", "y"],
    "3*x^3 - 3*x*y; 0; 2*x^2*y - 2*y^2 + x; x^2 - y; x^4 - y^2",
);

fn ring<F: F4Field + ParseCoefficient>(vars: &[&str], order: MonomialOrder) -> PolynomialRing<F> {
    PolynomialRing::new(vars.iter().copied(), order).unwrap()
}

fn check_system<F: F4Field + ParseCoefficient>(
    (vars, src): (&[&str], &str),
    order: MonomialOrder,
    members: &str,
    combination: &str,
) {
    let r = ring::<F>(vars, order);
    let generators = r.parse_many(src).unwrap();
    let ideal = Ideal::new(generators.clone()).unwrap();
    let lifted = LiftBasis::new(&generators).unwrap();
    assert_eq!(lifted.basis, ideal.basis());
    assert_eq!(
        lifted.basis,
        groebner_basis(generators.clone(), true).unwrap()
    );
    for (g, row) in lifted.basis.iter().zip(&lifted.cofactors) {
        assert!(verify_lift(&generators, row, g));
    }
    let mut targets = r.parse_many(members).unwrap();
    let multipliers = r.parse_many(combination).unwrap();
    targets.push(multipliers.iter().zip(&generators).fold(
        Polynomial::zero(generators[0].nvars, r.order()),
        |acc, (m, g)| acc.add(&m.multiply(g)),
    ));
    targets.extend(ideal.basis().iter().cloned());
    for f in targets {
        let h = ideal.lift(&f).unwrap().expect("member");
        assert_eq!(h.len(), generators.len());
        assert!(verify_lift(ideal.generators(), &h, &f));
    }
}

fn check_all<F: F4Field + ParseCoefficient>() {
    for order in [MonomialOrder::Lex, MonomialOrder::GRevLex] {
        check_system::<F>(
            CYCLIC3,
            order.clone(),
            "x^3 - 1; y^3 - 1",
            "x*y; z^2 - 1; 3*x - y",
        );
        check_system::<F>(
            REDUNDANT,
            order,
            "x; y; x^3 - x*y; x^2*y - y^2",
            "y; x^5; 1; x*y; 0",
        );
    }
    check_system::<F>(
        KATSURA3,
        MonomialOrder::GRevLex,
        "0",
        "x3; x0*x1 - 2; 7; x2^2",
    );
    check_system::<F>(
        CYCLIC4,
        MonomialOrder::GRevLex,
        "a*b*c*d - 1",
        "a*b; c - d; 1; a^2 - b*c",
    );
}

#[test]
fn lift_over_rationals() {
    check_all::<BigRational>();
}

#[test]
fn lift_over_prime_field() {
    check_all::<Zp>();
}

#[test]
fn lift_katsura3_lex_prime_field() {
    check_system::<Zp>(KATSURA3, MonomialOrder::Lex, "0", "x1; 1; x3^2; 5");
}

#[test]
fn non_member_has_no_certificate() {
    let r = ring::<BigRational>(CYCLIC3.0, MonomialOrder::GRevLex);
    let ideal = Ideal::new(r.parse_many(CYCLIC3.1).unwrap()).unwrap();
    for src in ["x", "x^2 + 1", "y - 1"] {
        let f = r.parse(src).unwrap();
        assert!(!ideal.contains(&f).unwrap());
        assert_eq!(ideal.lift(&f).unwrap(), None);
    }
}

#[test]
fn trivial_ideal_certifies_one() {
    let r = ring::<Zp>(&["x", "y"], MonomialOrder::Lex);
    let ideal = Ideal::new(r.parse_many("x*y - 1; x^2; y + x").unwrap()).unwrap();
    let one = Polynomial::constant(Zp::one(), 2, MonomialOrder::Lex);
    let h = ideal.lift(&one).unwrap().unwrap();
    assert!(verify_lift(ideal.generators(), &h, &one));
}

#[test]
fn verify_rejects_bad_certificates() {
    let r = ring::<BigRational>(&["x", "y"], MonomialOrder::GRevLex);
    let generators = r.parse_many("x^2 - y; y^2 - x").unwrap();
    let f = r.parse("x^4 - x").unwrap();
    let mut h = Ideal::new(generators.clone())
        .unwrap()
        .lift(&f)
        .unwrap()
        .unwrap();
    assert!(verify_lift(&generators, &h, &f));
    assert!(!verify_lift(&generators, &h[..1], &f));
    h[0] = h[0].add(&r.parse("1").unwrap());
    assert!(!verify_lift(&generators, &h, &f));
}

#[test]
fn division_identity() {
    let r = ring::<BigRational>(&["x", "y", "z"], MonomialOrder::GrLex);
    let divisors = r.parse_many("x*y - z; 2*y^2 + x; 0; z^2 - 1").unwrap();
    let f = r.parse("x^3*y^2 + 3*x*y*z^3 - y^4 + 7").unwrap();
    let (quotients, remainder) = f.divide(&divisors).unwrap();
    assert_eq!(remainder, f.reduce(&divisors).unwrap());
    assert!(quotients[2].is_zero());
    let sum = quotients
        .iter()
        .zip(&divisors)
        .fold(remainder, |acc, (q, d)| acc.add(&q.multiply(d)));
    assert_eq!(sum, f);
}
