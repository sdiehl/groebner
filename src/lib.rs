//! Groebner basis algorithms for multivariate polynomial ideals.
//!
//! The crate provides two public computation paths:
//!
//! - [`groebner_basis`] for the existing Buchberger implementation over any [`Field`].
//! - [`groebner_basis_f4_mod`] for a sparse F4-style implementation over [`PrimeField`].
//!
//! Polynomials can be built from strings using [`PolynomialRing`], where the variable list also
//! defines the lexicographic variable order.
//!
//! # Buchberger Example
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
//!
//! # F4 Example
//! ```
//! use groebner::{groebner_basis_f4_mod, MonomialOrder, PolynomialRing, PrimeField};
//!
//! type F32003 = PrimeField<32003>;
//!
//! let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::Lex)?;
//! let f1 = ring.parse("x^2 - y")?;
//! let f2 = ring.parse("x*y - 1")?;
//! let basis = groebner_basis_f4_mod(vec![f1, f2], F32003::modulus(), ring.order())?;
//!
//! assert!(!basis.is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod f4;
pub mod field;
pub mod finite_field;
pub mod grebauer_moller;
pub mod groebner;
pub mod monomial;
pub mod polynomial;
pub mod ring;
pub mod sugar;

pub use f4::groebner_basis_f4_mod;
pub use field::Field;
pub use finite_field::{PrimeField, PrimeFieldParseError};
pub use grebauer_moller::filter_gm_pairs;
pub use groebner::{
    groebner_basis, groebner_basis_with_strategy, is_groebner_basis, GroebnerError,
    SelectionStrategy,
};
pub use monomial::{Monomial, MonomialOrder};
pub use polynomial::{Polynomial, Term};
pub use ring::{ParsePolynomialError, PolynomialRing};
