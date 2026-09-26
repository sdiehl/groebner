use groebner::{MonomialOrder, PolynomialRing, SelectionStrategy, groebner_basis_with_strategy};
use num_rational::BigRational;

#[test]
fn test_groebner_basis_gm_strategy() {
    // Example: x^2 - y, xy - 1
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)
        .expect("test ring should be valid");
    let f1 = ring.parse("x^2 - y").expect("f1 should parse");
    let f2 = ring.parse("x*y - 1").expect("f2 should parse");
    let polynomials = vec![f1, f2];
    let result = groebner_basis_with_strategy(polynomials, true, &SelectionStrategy::GebauerMoller);
    assert!(
        result.is_ok(),
        "Gebauer–Möller strategy Groebner basis computation failed"
    );
    let gm_basis =
        result.expect("Gebauer–Möller strategy Groebner basis computation should succeed");
    assert!(
        !gm_basis.is_empty(),
        "Gebauer–Möller strategy produced empty basis"
    );
    for poly in &gm_basis {
        assert!(
            !poly.terms.is_empty(),
            "Gebauer–Möller strategy basis contains empty polynomial"
        );
    }
}
