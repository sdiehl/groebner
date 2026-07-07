use groebner::{groebner_basis_f4_mod, MonomialOrder, PolynomialRing, PrimeField};

type F32003 = PrimeField<32003>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::Lex)?;
    let f = ring.parse("x^2 - y")?;
    let g = ring.parse("x*y - 1")?;
    let basis = groebner_basis_f4_mod(vec![f, g], F32003::modulus(), ring.order())?;

    for (i, polynomial) in basis.iter().enumerate() {
        println!("g{}: {}", i + 1, ring.format(polynomial)?);
    }

    Ok(())
}
