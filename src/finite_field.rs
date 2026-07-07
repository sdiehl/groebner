//! Prime finite fields for modular Groebner basis computations.
//!
//! [`PrimeField<P>`] is a small const-generic field type for arithmetic modulo a machine prime.
//! It is the coefficient type used by the sparse F4 API.
//!
//! # Example
//! ```
//! use groebner::{Field, PrimeField};
//!
//! type F7 = PrimeField<7>;
//!
//! let a = F7::from(10_u32);
//! let b = F7::from(5_u32);
//! assert_eq!(a.value(), 3);
//! assert_eq!(a.multiply(&b).value(), 1);
//! ```

use crate::field::Field;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PrimeField<const P: u32> {
    value: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PrimeFieldParseError {
    InvalidModulus,
    InvalidInteger,
}

impl fmt::Display for PrimeFieldParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PrimeFieldParseError::InvalidModulus => {
                write!(f, "field modulus must be greater than 1")
            }
            PrimeFieldParseError::InvalidInteger => write!(f, "invalid finite field integer"),
        }
    }
}

impl std::error::Error for PrimeFieldParseError {}

impl<const P: u32> PrimeField<P> {
    pub const fn modulus() -> u32 {
        P
    }

    pub fn new(value: u32) -> Self {
        Self {
            value: Self::normalize(value),
        }
    }

    pub fn from_i64(value: i64) -> Self {
        if P <= 1 {
            return Self { value: 0 };
        }
        let modulus = i64::from(P);
        let normalized = value.rem_euclid(modulus);
        Self {
            value: normalized as u32,
        }
    }

    pub fn value(self) -> u32 {
        self.value
    }

    pub(crate) fn parse_digits_mod(digits: &str) -> Result<Self, PrimeFieldParseError> {
        if P <= 1 {
            return Err(PrimeFieldParseError::InvalidModulus);
        }
        if digits.is_empty() {
            return Err(PrimeFieldParseError::InvalidInteger);
        }

        let mut value = 0u64;
        let modulus = u64::from(P);
        for digit in digits.bytes() {
            if !digit.is_ascii_digit() {
                return Err(PrimeFieldParseError::InvalidInteger);
            }
            value = (value * 10 + u64::from(digit - b'0')) % modulus;
        }
        Ok(Self {
            value: value as u32,
        })
    }

    const fn normalize(value: u32) -> u32 {
        if P <= 1 {
            0
        } else {
            value % P
        }
    }
}

impl<const P: u32> Field for PrimeField<P> {
    fn zero() -> Self {
        Self { value: 0 }
    }

    fn one() -> Self {
        Self::new(1)
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }

    fn is_one(&self) -> bool {
        self.value == Self::one().value
    }

    fn add(&self, other: &Self) -> Self {
        if P <= 1 {
            return Self::zero();
        }
        let value = (u64::from(self.value) + u64::from(other.value)) % u64::from(P);
        Self {
            value: value as u32,
        }
    }

    fn subtract(&self, other: &Self) -> Self {
        if P <= 1 {
            return Self::zero();
        }
        let modulus = u64::from(P);
        let value = (modulus + u64::from(self.value) - u64::from(other.value)) % modulus;
        Self {
            value: value as u32,
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        if P <= 1 {
            return Self::zero();
        }
        let value = (u64::from(self.value) * u64::from(other.value)) % u64::from(P);
        Self {
            value: value as u32,
        }
    }

    fn negate(&self) -> Self {
        if self.is_zero() || P <= 1 {
            Self::zero()
        } else {
            Self {
                value: P - self.value,
            }
        }
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() || P <= 1 {
            return None;
        }

        let mut t = 0i64;
        let mut new_t = 1i64;
        let mut r = i64::from(P);
        let mut new_r = i64::from(self.value);

        while new_r != 0 {
            let quotient = r / new_r;
            let next_t = t - quotient * new_t;
            t = new_t;
            new_t = next_t;

            let next_r = r - quotient * new_r;
            r = new_r;
            new_r = next_r;
        }

        if r != 1 {
            return None;
        }

        Some(Self::from_i64(t))
    }
}

impl<const P: u32> fmt::Display for PrimeField<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

impl<const P: u32> From<u32> for PrimeField<P> {
    fn from(value: u32) -> Self {
        Self::new(value)
    }
}

impl<const P: u32> From<i32> for PrimeField<P> {
    fn from(value: i32) -> Self {
        Self::from_i64(i64::from(value))
    }
}

impl<const P: u32> FromStr for PrimeField<P> {
    type Err = PrimeFieldParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        let input = input.trim();
        let (negative, digits) = match input.strip_prefix('-') {
            Some(digits) => (true, digits),
            None => (false, input.strip_prefix('+').unwrap_or(input)),
        };
        let value = Self::parse_digits_mod(digits)?;
        if negative {
            Ok(value.negate())
        } else {
            Ok(value)
        }
    }
}
