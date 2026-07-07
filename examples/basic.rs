extern crate groebner;
use groebner::{groebner_basis, MonomialOrder, PolynomialRing};
use num_rational::BigRational;

fn main() {
    // Example: Compute a Groebner basis for the ideal (x^2 + y^2 - 1, x - y)
    let order = MonomialOrder::Lex;
    let ring = match PolynomialRing::<BigRational>::new(["x", "y"], order) {
        Ok(ring) => ring,
        Err(e) => {
            println!("Error creating polynomial ring: {e}");
            return;
        }
    };
    let f = match ring.parse("x^2 + y^2 - 1") {
        Ok(polynomial) => polynomial,
        Err(e) => {
            println!("Error parsing first polynomial: {e}");
            return;
        }
    };
    let g = match ring.parse("x - y") {
        Ok(polynomial) => polynomial,
        Err(e) => {
            println!("Error parsing second polynomial: {e}");
            return;
        }
    };

    // Compute the Groebner basis
    match groebner_basis(vec![f, g], order, true) {
        Ok(basis) => {
            println!("Groebner basis:");
            for (i, poly) in basis.iter().enumerate() {
                match ring.format(poly) {
                    Ok(formatted) => println!("g{}: {}", i + 1, formatted),
                    Err(e) => println!("Error formatting g{}: {e}", i + 1),
                }
            }
        }
        Err(e) => {
            println!("Error computing Groebner basis: {e}");
        }
    }
}
