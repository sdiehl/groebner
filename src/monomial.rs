//! Monomial types and orderings for Groebner basis computations
//!
//! This module provides the `Monomial` struct and `MonomialOrder` enum, which are used to
//! represent and compare monomials in multivariate polynomial rings. Monomial orderings
//! are essential for defining leading terms and for the correctness of Groebner basis algorithms.
//!
//! # Example
//! ```
//! use groebner::{Monomial, MonomialOrder};
//! let m1 = Monomial::new(vec![2, 1]); // x0^2 * x1
//! let m2 = Monomial::new(vec![1, 2]); // x0 * x1^2
//! assert_eq!(
//!     m1.compare(&m2, MonomialOrder::Lex),
//!     std::cmp::Ordering::Greater
//! );
//! ```

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt;
use std::sync::{Arc, Mutex, OnceLock, Weak};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MonomialOrder {
    Lex,
    GrLex,
    GRevLex,
}

impl MonomialOrder {
    pub fn compare(&self, a: &[u32], b: &[u32]) -> Ordering {
        match self {
            MonomialOrder::Lex => {
                for (ai, bi) in a.iter().zip(b.iter()) {
                    match ai.cmp(bi) {
                        Ordering::Equal => {}
                        other => return other,
                    }
                }
                Ordering::Equal
            }
            MonomialOrder::GrLex => {
                let deg_a: u32 = a.iter().sum();
                let deg_b: u32 = b.iter().sum();
                match deg_a.cmp(&deg_b) {
                    Ordering::Equal => {
                        for (ai, bi) in a.iter().zip(b.iter()) {
                            match ai.cmp(bi) {
                                Ordering::Equal => {}
                                other => return other,
                            }
                        }
                        Ordering::Equal
                    }
                    other => other,
                }
            }
            MonomialOrder::GRevLex => {
                let deg_a: u32 = a.iter().sum();
                let deg_b: u32 = b.iter().sum();
                match deg_a.cmp(&deg_b) {
                    Ordering::Equal => {
                        for (ai, bi) in a.iter().rev().zip(b.iter().rev()) {
                            match ai.cmp(bi) {
                                Ordering::Equal => {}
                                other => return other,
                            }
                        }
                        Ordering::Equal
                    }
                    other => other.reverse(),
                }
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Monomial {
    pub exponents: Arc<[u32]>,
}

impl Monomial {
    pub fn new(exponents: Vec<u32>) -> Self {
        Self {
            exponents: intern_exponents(exponents),
        }
    }
    pub fn one(nvars: usize) -> Self {
        Self::new(vec![0; nvars])
    }
    pub fn nvars(&self) -> usize {
        self.exponents.len()
    }
    pub fn degree(&self) -> u32 {
        self.exponents.iter().sum()
    }
    pub fn multiply(&self, other: &Self) -> Self {
        assert_eq!(self.nvars(), other.nvars());
        let exponents = self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .map(|(a, b)| a + b)
            .collect();
        Self::new(exponents)
    }
    pub fn divides(&self, other: &Self) -> bool {
        assert_eq!(self.nvars(), other.nvars());
        self.exponents
            .iter()
            .zip(other.exponents.iter())
            .all(|(a, b)| a <= b)
    }
    pub fn divide(&self, other: &Self) -> Option<Self> {
        assert_eq!(self.nvars(), other.nvars());
        if !other.divides(self) {
            return None;
        }
        let exponents = self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .map(|(a, b)| a - b)
            .collect();
        Some(Self::new(exponents))
    }
    pub fn lcm(&self, other: &Self) -> Self {
        assert_eq!(self.nvars(), other.nvars());
        let exponents = self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .map(|(a, b)| (*a).max(*b))
            .collect();
        Self::new(exponents)
    }
    pub fn compare(&self, other: &Self, order: MonomialOrder) -> Ordering {
        order.compare(&self.exponents, &other.exponents)
    }
}

fn intern_exponents(exponents: Vec<u32>) -> Arc<[u32]> {
    type ExponentInterner = Mutex<HashMap<Vec<u32>, Weak<[u32]>>>;

    static INTERNER: OnceLock<ExponentInterner> = OnceLock::new();

    let interner = INTERNER.get_or_init(|| Mutex::new(HashMap::new()));
    let Ok(mut interner) = interner.lock() else {
        return Arc::from(exponents);
    };

    if let Some(existing) = interner.get(&exponents) {
        if let Some(existing) = existing.upgrade() {
            return existing;
        }
    }

    let interned: Arc<[u32]> = Arc::from(exponents.clone());
    interner.insert(exponents, Arc::downgrade(&interned));
    interned
}

impl fmt::Display for Monomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.exponents.iter().all(|&e| e == 0) {
            return write!(f, "1");
        }
        let mut first = true;
        for (i, &exp) in self.exponents.iter().enumerate() {
            if exp > 0 {
                if !first {
                    write!(f, "*")?;
                }
                first = false;
                if exp == 1 {
                    write!(f, "x{i}")?;
                } else {
                    write!(f, "x{i}^{exp}")?;
                }
            }
        }
        Ok(())
    }
}
