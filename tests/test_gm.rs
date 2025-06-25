use groebner::{
    groebner_basis_with_strategy, Monomial, MonomialOrder, Polynomial, SelectionStrategy, Term,
};
use num_rational::BigRational;

#[test]
fn test_groebner_basis_gm_strategy() {
    // Example: x^2 - y, xy - 1
    let f1 = Polynomial::new(
        vec![
            Term::new(
                BigRational::new(1.into(), 1.into()),
                Monomial::new(vec![2, 0]),
            ),
            Term::new(
                BigRational::new((-1).into(), 1.into()),
                Monomial::new(vec![0, 1]),
            ),
        ],
        2,
        MonomialOrder::Lex,
    );
    let f2 = Polynomial::new(
        vec![
            Term::new(
                BigRational::new(1.into(), 1.into()),
                Monomial::new(vec![1, 1]),
            ),
            Term::new(
                BigRational::new((-1).into(), 1.into()),
                Monomial::new(vec![0, 0]),
            ),
        ],
        2,
        MonomialOrder::Lex,
    );
    let polynomials = vec![f1, f2];
    let result = groebner_basis_with_strategy(
        polynomials,
        MonomialOrder::Lex,
        true,
        &SelectionStrategy::GebauerMoller,
    );
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
