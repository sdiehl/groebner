//! Field trait and Rational number implementation
//!
//! This module defines the `Field` trait, which abstracts the algebraic concept of a field,
//! and provides a rational number implementation (`Rational`). Fields are used as coefficient
//! domains for polynomials in Groebner basis computations.
//!
//! # Extending
//! To use your own field type, implement the `Field` trait for your type.
//!
//! # Example
//! ```
//! use groebner::Field;
//! use num_rational::BigRational;
//! let a = BigRational::new(1.into(), 2.into());
//! let b = BigRational::new(1.into(), 3.into());
//! let sum = a.add(&b);
//! assert_eq!(sum, BigRational::new(5.into(), 6.into()));
//! ```

use num_bigint::BigInt;
use num_rational::BigRational;
use std::fmt;

pub trait Field: Clone + PartialEq + fmt::Debug + fmt::Display {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    #[must_use]
    fn add(&self, other: &Self) -> Self;
    #[must_use]
    fn subtract(&self, other: &Self) -> Self;
    #[must_use]
    fn multiply(&self, other: &Self) -> Self;
    #[must_use]
    fn negate(&self) -> Self;
    fn inverse(&self) -> Option<Self>;
    fn divide(&self, other: &Self) -> Option<Self> {
        other.inverse().map(|inv| self.multiply(&inv))
    }
}

impl Field for BigRational {
    fn zero() -> Self {
        BigRational::from_integer(BigInt::from(0))
    }
    fn one() -> Self {
        BigRational::from_integer(BigInt::from(1))
    }
    fn is_zero(&self) -> bool {
        <num_rational::Ratio<BigInt> as num_traits::Zero>::is_zero(self)
    }
    fn is_one(&self) -> bool {
        <num_rational::Ratio<BigInt> as num_traits::One>::is_one(self)
    }
    fn add(&self, other: &Self) -> Self {
        self + other
    }
    fn subtract(&self, other: &Self) -> Self {
        self - other
    }
    fn multiply(&self, other: &Self) -> Self {
        self * other
    }
    fn negate(&self) -> Self {
        -self
    }
    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(self.recip())
        }
    }
}
