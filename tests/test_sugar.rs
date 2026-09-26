use groebner::groebner::{SelectionStrategy, groebner_basis_with_strategy};
use groebner::{MonomialOrder, PolynomialRing};
use num_rational::BigRational;

#[test]
fn test_groebner_basis_sugar_strategy() {
    // Example: x^2 - y, xy - 1
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)
        .expect("test ring should be valid");
    let f1 = ring.parse("x^2 - y").expect("f1 should parse");
    let f2 = ring.parse("x*y - 1").expect("f2 should parse");
    let polynomials = vec![f1, f2];
    let result = groebner_basis_with_strategy(polynomials.clone(), true, &SelectionStrategy::Sugar);
    assert!(
        result.is_ok(),
        "Sugar strategy Groebner basis computation failed"
    );
    let sugar_basis = result.expect("Sugar strategy Groebner basis computation should succeed");
    // Also check the default strategy for comparison
    let default_result =
        groebner_basis_with_strategy(polynomials, true, &SelectionStrategy::Degree);
    assert!(
        default_result.is_ok(),
        "Default strategy Groebner basis computation failed"
    );
    let default_basis =
        default_result.expect("Default strategy Groebner basis computation should succeed");
    // Both should produce non-empty bases
    assert!(
        !sugar_basis.is_empty(),
        "Sugar strategy produced empty basis"
    );
    assert!(
        !default_basis.is_empty(),
        "Default strategy produced empty basis"
    );
    // The bases should be valid (contain polynomials)
    for poly in &sugar_basis {
        assert!(
            !poly.terms.is_empty(),
            "Sugar strategy basis contains empty polynomial"
        );
    }
    for poly in &default_basis {
        assert!(
            !poly.terms.is_empty(),
            "Default strategy basis contains empty polynomial"
        );
    }
}
