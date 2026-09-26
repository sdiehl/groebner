use groebner::{MonomialOrder, PolynomialRing, PrimeField, groebner_basis_f4};

type F32003 = PrimeField<32003>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::Lex)?;
    let f = ring.parse("x^2 - y")?;
    let g = ring.parse("x*y - 1")?;
    let basis = groebner_basis_f4(vec![f, g], true)?;

    for (i, polynomial) in basis.iter().enumerate() {
        println!("g{}: {}", i + 1, ring.format(polynomial)?);
    }

    Ok(())
}
