//! Prime fields: [`PrimeField<P>`] with a const modulus and [`Zp`] with a runtime modulus.
//!
//! ```
//! use groebner::{Field, PrimeField, Zp};
//!
//! type F7 = PrimeField<7>;
//! let a = F7::from(10_u32);
//! assert_eq!(a.multiply(&F7::from(5_u32)).value(), 1);
//!
//! let p = 18446744073709551557;
//! let b = Zp::new(3, p);
//! assert_eq!(b.multiply(&b.inverse().unwrap()), Zp::new(1, p));
//! ```

use crate::field::{Field, ModularField};
use std::fmt;
use std::hash::{Hash, Hasher};
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

fn parse_digits_mod(digits: &str, modulus: u64) -> Result<u64, PrimeFieldParseError> {
    if modulus <= 1 {
        return Err(PrimeFieldParseError::InvalidModulus);
    }
    if digits.is_empty() {
        return Err(PrimeFieldParseError::InvalidInteger);
    }
    let mut value = 0u128;
    for digit in digits.bytes() {
        if !digit.is_ascii_digit() {
            return Err(PrimeFieldParseError::InvalidInteger);
        }
        value = (value * 10 + u128::from(digit - b'0')) % u128::from(modulus);
    }
    Ok(value as u64)
}

fn mod_inverse(value: u64, modulus: u64) -> Option<u64> {
    if value == 0 || modulus <= 1 {
        return None;
    }
    let (mut t, mut new_t) = (0i128, 1i128);
    let (mut r, mut new_r) = (i128::from(modulus), i128::from(value));
    while new_r != 0 {
        let q = r / new_r;
        (t, new_t) = (new_t, t - q * new_t);
        (r, new_r) = (new_r, r - q * new_r);
    }
    (r == 1).then(|| t.rem_euclid(i128::from(modulus)) as u64)
}

fn mul_mod(a: u64, b: u64, modulus: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(modulus)) as u64
}

impl<const P: u32> PrimeField<P> {
    pub const fn modulus() -> u32 {
        P
    }

    pub fn new(value: u32) -> Self {
        Self {
            value: if P <= 1 { 0 } else { value % P },
        }
    }

    pub fn from_i64(value: i64) -> Self {
        if P <= 1 {
            return Self { value: 0 };
        }
        Self {
            value: value.rem_euclid(i64::from(P)) as u32,
        }
    }

    pub fn value(self) -> u32 {
        self.value
    }

    pub(crate) fn parse_digits(digits: &str) -> Result<Self, PrimeFieldParseError> {
        parse_digits_mod(digits, u64::from(P)).map(|value| Self {
            value: value as u32,
        })
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
        let sum = u64::from(self.value) + u64::from(other.value);
        Self {
            value: (sum % u64::from(P)) as u32,
        }
    }
    fn subtract(&self, other: &Self) -> Self {
        self.add(&other.negate())
    }
    fn multiply(&self, other: &Self) -> Self {
        if P <= 1 {
            return Self::zero();
        }
        Self {
            value: mul_mod(u64::from(self.value), u64::from(other.value), u64::from(P)) as u32,
        }
    }
    fn negate(&self) -> Self {
        if self.value == 0 || P <= 1 {
            Self::zero()
        } else {
            Self {
                value: P - self.value,
            }
        }
    }
    fn inverse(&self) -> Option<Self> {
        mod_inverse(u64::from(self.value), u64::from(P)).map(|v| Self { value: v as u32 })
    }
}

impl<const P: u32> ModularField for PrimeField<P> {
    fn modulus(&self) -> u64 {
        u64::from(P)
    }
    fn residue(&self) -> u64 {
        u64::from(self.value)
    }
    fn from_residue(residue: u64, _modulus: u64) -> Self {
        Self {
            value: (residue % u64::from(P)) as u32,
        }
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
        let (negative, digits) = split_sign(input);
        let value = Self::parse_digits(digits)?;
        Ok(if negative { value.negate() } else { value })
    }
}

fn split_sign(input: &str) -> (bool, &str) {
    let input = input.trim();
    match input.strip_prefix('-') {
        Some(digits) => (true, digits),
        None => (false, input.strip_prefix('+').unwrap_or(input)),
    }
}

/// Element of `GF(p)` for a runtime prime `p` up to 64 bits.
///
/// [`Field::zero`] and [`Field::one`] produce placeholder constants with modulus 0 that adopt the
/// modulus of whatever they are combined with, so the nullary `Field` constructors still work.
#[derive(Debug, Clone, Copy)]
pub struct Zp {
    value: u64,
    modulus: u64,
}

impl Zp {
    pub fn new(value: u64, modulus: u64) -> Self {
        Self {
            value: if modulus <= 1 { 0 } else { value % modulus },
            modulus,
        }
    }

    pub fn from_i64(value: i64, modulus: u64) -> Self {
        Self {
            value: if modulus <= 1 {
                0
            } else {
                i128::from(value).rem_euclid(i128::from(modulus)) as u64
            },
            modulus,
        }
    }

    pub fn value(self) -> u64 {
        self.value
    }

    pub fn modulus(self) -> u64 {
        self.modulus
    }

    /// Bind a placeholder constant to `modulus`; bound elements are returned unchanged.
    pub fn bind(self, modulus: u64) -> Self {
        if self.modulus != 0 {
            self
        } else {
            Self::from_i64(self.value as i64, modulus)
        }
    }

    pub(crate) fn parse_digits(digits: &str, modulus: u64) -> Result<Self, PrimeFieldParseError> {
        parse_digits_mod(digits, modulus).map(|value| Self { value, modulus })
    }

    fn unify(self, other: Self) -> (Self, Self, u64) {
        let modulus = self.modulus.max(other.modulus);
        (self.bind(modulus), other.bind(modulus), modulus)
    }
}

impl PartialEq for Zp {
    fn eq(&self, other: &Self) -> bool {
        let (a, b, _) = self.unify(*other);
        a.value == b.value
    }
}

impl Eq for Zp {}

impl Hash for Zp {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.value.hash(state);
    }
}

impl Field for Zp {
    fn zero() -> Self {
        Self {
            value: 0,
            modulus: 0,
        }
    }
    fn one() -> Self {
        Self {
            value: 1,
            modulus: 0,
        }
    }
    fn is_zero(&self) -> bool {
        self.value == 0
    }
    fn is_one(&self) -> bool {
        self.value == 1
    }
    fn add(&self, other: &Self) -> Self {
        let (a, b, modulus) = self.unify(*other);
        if modulus == 0 {
            return Self {
                value: a.value.wrapping_add(b.value),
                modulus,
            };
        }
        let value = if a.value >= modulus - b.value {
            a.value - (modulus - b.value)
        } else {
            a.value + b.value
        };
        Self { value, modulus }
    }
    fn subtract(&self, other: &Self) -> Self {
        self.add(&other.negate())
    }
    fn multiply(&self, other: &Self) -> Self {
        let (a, b, modulus) = self.unify(*other);
        if modulus == 0 {
            return Self {
                value: a.value.wrapping_mul(b.value),
                modulus,
            };
        }
        Self {
            value: mul_mod(a.value, b.value, modulus),
            modulus,
        }
    }
    fn negate(&self) -> Self {
        if self.modulus == 0 {
            return Self {
                value: 0u64.wrapping_sub(self.value),
                modulus: 0,
            };
        }
        Self {
            value: if self.value == 0 {
                0
            } else {
                self.modulus - self.value
            },
            modulus: self.modulus,
        }
    }
    fn inverse(&self) -> Option<Self> {
        if self.modulus == 0 {
            return (self.value == 1 || self.value == u64::MAX).then_some(*self);
        }
        mod_inverse(self.value, self.modulus).map(|value| Self {
            value,
            modulus: self.modulus,
        })
    }
}

impl ModularField for Zp {
    fn modulus(&self) -> u64 {
        self.modulus
    }
    fn residue(&self) -> u64 {
        self.value
    }
    fn from_residue(residue: u64, modulus: u64) -> Self {
        Self::new(residue, modulus)
    }
}

impl fmt::Display for Zp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.modulus == 0 {
            write!(f, "{}", self.value as i64)
        } else {
            write!(f, "{}", self.value)
        }
    }
}
