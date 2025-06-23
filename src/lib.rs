//! Groebner Basis Library
//!
//! This library implements algorithms for computing Groebner bases of polynomial ideals,
//! including the Buchberger algorithm and F4 algorithm optimizations.
//!
//! # Example
//! ```
//! use num_rational::BigRational;
//! use groebner::{groebner_basis, Polynomial, MonomialOrder, Term, Monomial};
//! // x^2 - y, xy - 1
//! let f1 = Polynomial::new(
//!     vec![
//!         Term::new(BigRational::new(1.into(), 1.into()), Monomial::new(vec![2, 0])),
//!         Term::new(BigRational::new((-1).into(), 1.into()), Monomial::new(vec![0, 1]))
//!     ],
//!     2,
//!     MonomialOrder::Lex
//! );
//! let f2 = Polynomial::new(
//!     vec![
//!         Term::new(BigRational::new(1.into(), 1.into()), Monomial::new(vec![1, 1])),
//!         Term::new(BigRational::new((-1).into(), 1.into()), Monomial::new(vec![0, 0]))
//!     ],
//!     2,
//!     MonomialOrder::Lex
//! );
//! let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true);
//! assert!(!basis.is_empty());
//! ```

pub mod field;
pub mod groebner;
pub mod monomial;
pub mod polynomial;

pub use field::Field;
pub use groebner::{groebner_basis, is_groebner_basis};
pub use monomial::{Monomial, MonomialOrder};
pub use polynomial::{Polynomial, Term};

// (All other code has been moved to modules)
