#![allow(clippy::expect_used, clippy::unwrap_used)]

use groebner::{
    groebner_basis, groebner_basis_f4, is_groebner_basis, Field, MonomialOrder, PolynomialRing,
    PrimeField, Zp,
};

const P64: u64 = 18_446_744_073_709_551_557;

#[test]
fn arithmetic_with_large_modulus() {
    let a = Zp::new(P64 - 1, P64);
    assert_eq!(a.multiply(&a), Zp::new(1, P64));
    assert_eq!(a.add(&Zp::one()), Zp::new(0, P64));
    assert_eq!(a.multiply(&a.inverse().unwrap()), Zp::one());
    assert_eq!(Zp::from_i64(-3, 7), Zp::new(4, 7));
}

#[test]
fn placeholders_adopt_modulus() {
    let one = Zp::one();
    assert_eq!(one.modulus(), 0);
    assert_eq!(one.negate().add(&Zp::new(2, 7)), Zp::new(1, 7));
    assert_eq!(one.negate().bind(7), Zp::new(6, 7));
}

#[test]
fn parse_requires_modulus() {
    let ring = PolynomialRing::<Zp>::new(["x"], MonomialOrder::Lex).unwrap();
    assert!(ring.parse("2*x").is_err());
    let ring = PolynomialRing::<Zp>::with_modulus(["x"], MonomialOrder::Lex, 7).unwrap();
    let p = ring.parse("9*x - 1/2").unwrap();
    assert_eq!(ring.format(&p).unwrap(), "2*x + 3");
    assert_eq!(ring.modulus(), Some(7));
}

#[test]
fn f4_matches_const_modulus_field() {
    let zp =
        PolynomialRing::<Zp>::with_modulus(["x", "y", "z"], MonomialOrder::GRevLex, 32003).unwrap();
    let pf =
        PolynomialRing::<PrimeField<32003>>::new(["x", "y", "z"], MonomialOrder::GRevLex).unwrap();
    let system = "x^2 + y^2 + z^2 - 1; x*y - z; x - y^2";
    let a = groebner_basis_f4(zp.parse_many(system).unwrap(), true).unwrap();
    let b = groebner_basis_f4(pf.parse_many(system).unwrap(), true).unwrap();
    let a_text: Vec<_> = a.iter().map(|p| zp.format(p).unwrap()).collect();
    let b_text: Vec<_> = b.iter().map(|p| pf.format(p).unwrap()).collect();
    assert_eq!(a_text, b_text);
    assert!(is_groebner_basis(&a).unwrap());
}

#[test]
fn f4_with_64_bit_modulus_matches_buchberger() {
    let ring = PolynomialRing::<Zp>::with_modulus(["x", "y"], MonomialOrder::Lex, P64).unwrap();
    let polys = ring.parse_many("x^2 - 3*y; x*y - 1/5").unwrap();
    let f4 = groebner_basis_f4(polys.clone(), true).unwrap();
    let buchberger = groebner_basis(polys, true).unwrap();
    assert_eq!(f4, buchberger);
    assert!(is_groebner_basis(&f4).unwrap());
}
