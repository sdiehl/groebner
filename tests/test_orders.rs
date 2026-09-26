#![allow(clippy::expect_used, clippy::unwrap_used)]

use groebner::{
    Monomial, MonomialOrder, PolynomialRing, PrimeField, groebner_basis, groebner_basis_f4,
    is_groebner_basis,
};
use std::cmp::Ordering;

#[test]
fn weighted_order_compares_by_weight_first() {
    let order = MonomialOrder::weighted(vec![1, 3], MonomialOrder::Lex);
    let a = Monomial::new(vec![2, 0]);
    let b = Monomial::new(vec![0, 1]);
    assert_eq!(a.compare(&b, &order), Ordering::Less);
    assert_eq!(a.compare(&b, &MonomialOrder::GrLex), Ordering::Greater);
}

#[test]
fn block_order_eliminates_first_block() {
    let order = MonomialOrder::elimination(1, 2);
    assert!(order.eliminates(1));
    assert!(!order.eliminates(2));
    assert!(MonomialOrder::Lex.eliminates(2));
    assert!(!MonomialOrder::GRevLex.eliminates(1));
    let a = Monomial::new(vec![1, 0, 0]);
    let b = Monomial::new(vec![0, 5, 5]);
    assert_eq!(a.compare(&b, &order), Ordering::Greater);
}

#[test]
fn elimination_order_groebner_basis_projects() {
    let ring =
        PolynomialRing::<PrimeField<32003>>::new(["t", "x", "y"], MonomialOrder::elimination(1, 2))
            .unwrap();
    let polys = ring.parse_many("x - t^2; y - t^3").unwrap();
    let f4 = groebner_basis_f4(polys.clone(), true).unwrap();
    let buchberger = groebner_basis(polys, true).unwrap();
    assert_eq!(f4, buchberger);
    assert!(is_groebner_basis(&f4).unwrap());
    let projected: Vec<_> = f4
        .iter()
        .filter(|p| p.terms.iter().all(|t| t.monomial.exponents()[0] == 0))
        .collect();
    assert_eq!(projected.len(), 1);
    assert_eq!(ring.format(projected[0]).unwrap(), "x^3 + 32002*y^2");
}

#[test]
fn weighted_order_groebner_basis_is_valid() {
    let order = MonomialOrder::weighted(vec![2, 3], MonomialOrder::GRevLex);
    let ring = PolynomialRing::<PrimeField<32003>>::new(["x", "y"], order).unwrap();
    let polys = ring.parse_many("y^2 - x^3; x*y - 1").unwrap();
    let basis = groebner_basis_f4(polys.clone(), true).unwrap();
    assert!(is_groebner_basis(&basis).unwrap());
    assert_eq!(basis, groebner_basis(polys, true).unwrap());
}
