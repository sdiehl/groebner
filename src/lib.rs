//! Groebner bases for multivariate polynomial ideals.
//!
//! - [`groebner_basis_f4`]: F4 over any [`F4Field`], with a fast dense modular path for
//!   [`PrimeField`] and [`Zp`].
//! - [`groebner_basis`], [`groebner_basis_incremental`], [`groebner_basis_with_strategy`] and
//!   [`groebner_basis_parallel`]: Buchberger over any [`Field`].
//! - [`Ideal`]: normal forms, membership, elimination, dimension, radical membership and
//!   change of order via [`fglm()`], and membership certificates via [`Ideal::lift`].
//! - [`PolynomialRing`]: parsing and formatting with named variables under any [`MonomialOrder`].
//! - [`RationalFunction`]: coefficients in Q(a) for systems with one symbolic parameter.
//!
//! # F4 over a prime field
//! ```
//! use groebner::{groebner_basis_f4, is_groebner_basis, MonomialOrder, PolynomialRing, PrimeField};
//!
//! let ring = PolynomialRing::<PrimeField<32003>>::new(["x", "y"], MonomialOrder::GRevLex)?;
//! let polys = ring.parse_many("x^2 - y; x*y - 1")?;
//! let basis = groebner_basis_f4(polys, true)?;
//! assert!(is_groebner_basis(&basis)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Runtime modulus
//! ```
//! use groebner::{groebner_basis_f4, MonomialOrder, PolynomialRing, Zp};
//!
//! let ring = PolynomialRing::<Zp>::with_modulus(["x", "y"], MonomialOrder::Lex, 1_000_003)?;
//! let basis = groebner_basis_f4(ring.parse_many("x^2 - y; x*y - 1")?, true)?;
//! assert_eq!(basis.len(), 2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Ideals and change of order
//! ```
//! use groebner::{Ideal, MonomialOrder, PolynomialRing};
//! use num_rational::BigRational;
//!
//! let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::GRevLex)?;
//! let ideal = Ideal::new(ring.parse_many("x^2 + y^2 - 1; x - y")?)?;
//! assert_eq!(ideal.vector_space_dimension(), Some(2));
//! let lex = ideal.change_order(MonomialOrder::Lex)?;
//! assert!(lex.contains(&ring.parse("2*y^2 - 1")?)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

pub mod f4;
pub mod modular;
mod par;
pub use modular::{RationalOptions, groebner_basis_f4_rational, groebner_basis_f4_rational_with};
pub mod fglm;
pub use polycore::Modular as ModularField;
pub use polycore::{Field, Fp, Gf, Modular};
pub mod finite_field;
pub mod grebauer_moller;
pub mod groebner;
pub mod ideal;
pub mod lift;
pub mod monomial;
#[cfg(feature = "parameters")]
pub mod parameters;
pub mod polynomial;
pub mod rational_function;
pub mod ring;
pub mod sugar;

pub use f4::{F4Field, Reducers, SparseRow, groebner_basis_f4};
pub use fglm::{fglm, is_zero_dimensional, standard_monomials};
pub use finite_field::{PrimeField, PrimeFieldParseError, Zp};
pub use grebauer_moller::filter_gm_pairs;
pub use groebner::{
    GroebnerError, SelectionStrategy, groebner_basis, groebner_basis_incremental,
    groebner_basis_with_strategy, is_groebner_basis,
};
#[cfg(feature = "parallel")]
pub use groebner::{groebner_basis_parallel, is_groebner_basis_parallel};
pub use ideal::Ideal;
pub use lift::{LiftBasis, verify_lift};
pub use monomial::{Monomial, MonomialExt, MonomialOrder};
#[cfg(feature = "parameters")]
pub use parameters::Frac;
pub use polycore::{Order, Poly, Ring};
pub use polynomial::{Polynomial, PolynomialExt, Term};
pub use rational_function::{RationalFunction, specialize};
pub use ring::{ParseCoefficient, ParsePolynomialError, PolynomialRing};
