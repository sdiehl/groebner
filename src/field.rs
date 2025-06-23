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
//! use groebner::{Rational, Field};
//! let a = Rational::new(1, 2);
//! let b = Rational::new(1, 3);
//! let sum = a.add(&b);
//! assert_eq!(sum, Rational::new(5, 6));
//! ```

use std::fmt;

pub trait Field: Clone + PartialEq + fmt::Debug + fmt::Display {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn add(&self, other: &Self) -> Self;
    fn subtract(&self, other: &Self) -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn negate(&self) -> Self;
    fn inverse(&self) -> Self;
    fn divide(&self, other: &Self) -> Self {
        self.multiply(&other.inverse())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Rational {
    pub numerator: i64,
    pub denominator: i64,
}

impl Rational {
    pub fn new(num: i64, den: i64) -> Self {
        if den == 0 {
            panic!("Division by zero");
        }
        let g = gcd(num.abs(), den.abs());
        let (num, den) = if den < 0 { (-num, -den) } else { (num, den) };
        Self {
            numerator: num / g,
            denominator: den / g,
        }
    }
    pub fn from_integer(n: i64) -> Self {
        Self::new(n, 1)
    }
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

impl Field for Rational {
    fn zero() -> Self {
        Self::new(0, 1)
    }
    fn one() -> Self {
        Self::new(1, 1)
    }
    fn is_zero(&self) -> bool {
        self.numerator == 0
    }
    fn is_one(&self) -> bool {
        self.numerator == 1 && self.denominator == 1
    }
    fn add(&self, other: &Self) -> Self {
        Self::new(
            self.numerator * other.denominator + other.numerator * self.denominator,
            self.denominator * other.denominator,
        )
    }
    fn subtract(&self, other: &Self) -> Self {
        Self::new(
            self.numerator * other.denominator - other.numerator * self.denominator,
            self.denominator * other.denominator,
        )
    }
    fn multiply(&self, other: &Self) -> Self {
        Self::new(
            self.numerator * other.numerator,
            self.denominator * other.denominator,
        )
    }
    fn negate(&self) -> Self {
        Self::new(-self.numerator, self.denominator)
    }
    fn inverse(&self) -> Self {
        if self.is_zero() {
            panic!("Division by zero");
        }
        Self::new(self.denominator, self.numerator)
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}
