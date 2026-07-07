//! Groebner basis algorithms for multivariate polynomial ideals.
//!
//! The crate provides two public computation paths:
//!
//! - [`groebner_basis`] for the existing Buchberger implementation over any [`Field`].
//! - [`groebner_basis_incremental`] for extending an existing Buchberger basis with new generators.
//! - [`groebner_basis_parallel`] for Rayon-backed Buchberger batches when the `parallel` feature is enabled.
//! - [`groebner_basis_f4_mod`] for a sparse F4-style implementation over [`PrimeField`].
//!
//! Polynomials can be built from strings using [`PolynomialRing`], where the variable list also
//! defines the lexicographic variable order. `PolynomialRing` also formats results with variable
//! names through [`PolynomialRing::format`] and [`PolynomialRing::format_latex`].
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
//! # Parallel Buchberger Example
//! ```
//! # #[cfg(feature = "parallel")]
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! use groebner::{groebner_basis_parallel, MonomialOrder, PolynomialRing};
//! use num_rational::BigRational;
//!
//! let ring = PolynomialRing::<BigRational>::new(["x", "y", "z"], MonomialOrder::GrLex)?;
//! let f1 = ring.parse("x^2 + y^2 + z^2 - 1")?;
//! let f2 = ring.parse("x*y - z")?;
//! let basis = groebner_basis_parallel(vec![f1, f2], ring.order(), true)?;
//!
//! assert!(!basis.is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! # }
//! # #[cfg(not(feature = "parallel"))]
//! # fn main() {}
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
    groebner_basis, groebner_basis_incremental, groebner_basis_with_strategy, is_groebner_basis,
    GroebnerError, SelectionStrategy,
};
#[cfg(feature = "parallel")]
pub use groebner::{groebner_basis_parallel, is_groebner_basis_parallel};
pub use monomial::{Monomial, MonomialOrder};
pub use polynomial::{Polynomial, Term};
pub use ring::{ParsePolynomialError, PolynomialRing};
