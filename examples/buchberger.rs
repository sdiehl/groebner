use groebner::{groebner_basis, MonomialOrder, PolynomialRing};
use num_rational::BigRational;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)?;
    let f = ring.parse("x^2 + y^2 - 1")?;
    let g = ring.parse("x - y")?;
    let basis = groebner_basis(vec![f, g], ring.order(), true)?;

    for (i, polynomial) in basis.iter().enumerate() {
        println!("g{}: {}", i + 1, ring.format(polynomial)?);
    }

    Ok(())
}
