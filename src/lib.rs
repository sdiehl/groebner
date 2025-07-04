//! Groebner Basis
//!
//! This library computes Groebner bases of polynomial ideals, using the Buchberger and F4
//! algorithms.
//!
//! # Example
//! ```
//! use groebner::is_groebner_basis;
//! use groebner::{groebner_basis, Monomial, MonomialOrder, Polynomial, Term};
//! use num_rational::BigRational;
//! // x^2 - y, xy - 1
//! let f1 = Polynomial::new(
//!     vec![
//!         Term::new(
//!             BigRational::new(1.into(), 1.into()),
//!             Monomial::new(vec![2, 0]),
//!         ),
//!         Term::new(
//!             BigRational::new((-1).into(), 1.into()),
//!             Monomial::new(vec![0, 1]),
//!         ),
//!     ],
//!     2,
//!     MonomialOrder::Lex,
//! );
//! let f2 = Polynomial::new(
//!     vec![
//!         Term::new(
//!             BigRational::new(1.into(), 1.into()),
//!             Monomial::new(vec![1, 1]),
//!         ),
//!         Term::new(
//!             BigRational::new((-1).into(), 1.into()),
//!             Monomial::new(vec![0, 0]),
//!         ),
//!     ],
//!     2,
//!     MonomialOrder::Lex,
//! );
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
//! ```

pub mod field;
pub mod grebauer_moller;
pub mod groebner;
pub mod monomial;
pub mod polynomial;
pub mod sugar;

pub use field::Field;
pub use grebauer_moller::filter_gm_pairs;
pub use groebner::{
    groebner_basis, groebner_basis_with_strategy, is_groebner_basis, GroebnerError,
    SelectionStrategy,
};
pub use monomial::{Monomial, MonomialOrder};
pub use polynomial::{Polynomial, Term};
