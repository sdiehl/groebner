//! Groebner Basis
//!
//! This library computes Groebner bases of polynomial ideals, using the Buchberger and F4
//! algorithms.
//!
//! # Example
//! ```
//! use groebner::{groebner_basis, is_groebner_basis, MonomialOrder, PolynomialRing};
//! use num_rational::BigRational;
//!
//! let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::Lex)?;
//! let f1 = ring.parse("x^2 - y")?;
//! let f2 = ring.parse("x*y - 1")?;
//! let basis_result = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true);
//! match basis_result {
//!     Ok(basis) => {
//!         assert!(!basis.is_empty());
//!         match is_groebner_basis(&basis) {
//!             Ok(true) => {}
//!             Ok(false) => panic!("Basis is not a Groebner basis!"),
//!             Err(e) => panic!("Groebner basis check failed: {}", e),
//!         }
//!     }
//!     Err(e) => panic!("Groebner basis computation failed: {}", e),
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod field;
pub mod grebauer_moller;
pub mod groebner;
pub mod monomial;
pub mod polynomial;
pub mod ring;
pub mod sugar;

pub use field::Field;
pub use grebauer_moller::filter_gm_pairs;
pub use groebner::{
    groebner_basis, groebner_basis_with_strategy, is_groebner_basis, GroebnerError,
    SelectionStrategy,
};
pub use monomial::{Monomial, MonomialOrder};
pub use polynomial::{Polynomial, Term};
pub use ring::{ParsePolynomialError, PolynomialRing};
