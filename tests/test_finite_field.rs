use groebner::{groebner_basis, Field, MonomialOrder, PolynomialRing, PrimeField};

type F32003 = PrimeField<32003>;

#[test]
fn prime_field_arithmetic_reduces_modulo_p() {
    let a = F32003::from(32005u32);
    let b = F32003::from(7u32);

    assert_eq!(a.value(), 2);
    assert_eq!(a.add(&b).value(), 9);
    assert_eq!(a.subtract(&b).value(), 31998);
    assert_eq!(a.multiply(&b).value(), 14);
    assert_eq!(
        b.multiply(&b.inverse().expect("nonzero has inverse"))
            .value(),
        1
    );
}

#[test]
fn polynomial_ring_parses_prime_field_coefficients() {
    let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::Lex)
        .expect("ring should be valid");
    let polynomial = ring
        .parse("-2*x^2 + 1/2*y - 32004")
        .expect("polynomial should parse");

    assert_eq!(
        ring.format(&polynomial).expect("should format"),
        "32001*x^2 + 16002*y + 32002"
    );
}

#[test]
fn computes_small_groebner_basis_over_prime_field() {
    let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::Lex)
        .expect("ring should be valid");
    let f1 = ring.parse("x^2 - y").expect("f1 should parse");
    let f2 = ring.parse("x*y - 1").expect("f2 should parse");

    let basis = groebner_basis(vec![f1, f2], ring.order(), true)
        .expect("Groebner basis should compute over GF(p)");

    assert!(!basis.is_empty());
}
