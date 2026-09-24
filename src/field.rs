//! The [`Field`] trait and its rational implementation.
//!
//! ```
//! use groebner::Field;
//! use num_rational::BigRational;
//! let a = BigRational::new(1.into(), 2.into());
//! let b = BigRational::new(1.into(), 3.into());
//! assert_eq!(a.add(&b), BigRational::new(5.into(), 6.into()));
//! ```

use num_rational::BigRational;
use std::fmt;

/// Coefficient field for polynomials. Implement this to use your own coefficient type.
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

/// A prime field whose elements fit in a machine word, used by the fast F4 linear algebra.
pub trait ModularField: Field + Copy + Send + Sync {
    /// The prime, or 0 for a placeholder constant that has not been bound to a modulus yet.
    fn modulus(&self) -> u64;
    fn residue(&self) -> u64;
    fn from_residue(residue: u64, modulus: u64) -> Self;
    /// Residue in `[0, modulus)`, mapping an unbound placeholder as a signed integer.
    fn residue_mod(&self, modulus: u64) -> u64 {
        if self.modulus() == 0 {
            i128::from(self.residue() as i64).rem_euclid(i128::from(modulus)) as u64
        } else {
            self.residue()
        }
    }
}

impl Field for BigRational {
    fn zero() -> Self {
        <Self as num_traits::Zero>::zero()
    }
    fn one() -> Self {
        <Self as num_traits::One>::one()
    }
    fn is_zero(&self) -> bool {
        <Self as num_traits::Zero>::is_zero(self)
    }
    fn is_one(&self) -> bool {
        <Self as num_traits::One>::is_one(self)
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
