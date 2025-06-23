extern crate groebner;
use groebner::{groebner_basis, Monomial, MonomialOrder, Polynomial, Term};
use num_rational::BigRational;

fn main() {
    // Example: Compute a Groebner basis for the ideal (x^2 + y^2 - 1, x - y)
    // Variables: x = 0, y = 1
    let nvars = 2;
    let order = MonomialOrder::Lex;

    // x^2 + y^2 - 1
    let f = Polynomial::new(
        vec![
            Term::new(
                BigRational::new(1.into(), 1.into()),
                Monomial::new(vec![2, 0]),
            ), // x^2
            Term::new(
                BigRational::new(1.into(), 1.into()),
                Monomial::new(vec![0, 2]),
            ), // y^2
            Term::new(
                BigRational::new((-1).into(), 1.into()),
                Monomial::new(vec![0, 0]),
            ), // -1
        ],
        nvars,
        order,
    );

    // x - y
    let g = Polynomial::new(
        vec![
            Term::new(
                BigRational::new(1.into(), 1.into()),
                Monomial::new(vec![1, 0]),
            ), // x
            Term::new(
                BigRational::new((-1).into(), 1.into()),
                Monomial::new(vec![0, 1]),
            ), // -y
        ],
        nvars,
        order,
    );

    // Compute the Groebner basis
    let basis = groebner_basis(vec![f, g], order, true);

    // Print the Groebner basis
    println!("Groebner basis:");
    for (i, poly) in basis.iter().enumerate() {
        println!("g{}: {:?}", i + 1, poly);
    }
}
