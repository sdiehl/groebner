//! Monomials and monomial orders.
//!
//! ```
//! use groebner::{Monomial, MonomialOrder};
//! let m1 = Monomial::new(vec![2, 1]);
//! let m2 = Monomial::new(vec![1, 2]);
//! assert_eq!(m1.compare(&m2, &MonomialOrder::Lex), std::cmp::Ordering::Greater);
//! ```

use std::cmp::Ordering;
use std::fmt;
use std::sync::Arc;

/// A monomial order on exponent vectors.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MonomialOrder {
    Lex,
    GrLex,
    GRevLex,
    /// Weighted degree first, ties broken by `tie_break`.
    Weighted {
        weights: Arc<[u32]>,
        tie_break: Arc<MonomialOrder>,
    },
    /// Product order: each block of `usize` variables is compared with its own order, left to right.
    Block(Arc<[(MonomialOrder, usize)]>),
}

impl MonomialOrder {
    /// Weighted degree order with the given weights and tie-breaking order.
    pub fn weighted(weights: impl Into<Arc<[u32]>>, tie_break: MonomialOrder) -> Self {
        Self::Weighted {
            weights: weights.into(),
            tie_break: Arc::new(tie_break),
        }
    }

    /// Product order from a list of `(order, block size)` pairs.
    pub fn block(blocks: impl Into<Arc<[(MonomialOrder, usize)]>>) -> Self {
        Self::Block(blocks.into())
    }

    /// Elimination order for the first `first` variables: GRevLex on those, then GRevLex on the remaining `rest`.
    pub fn elimination(first: usize, rest: usize) -> Self {
        Self::block(vec![(Self::GRevLex, first), (Self::GRevLex, rest)])
    }

    /// True when the order eliminates the first `k` variables, so a Groebner basis restricts to an elimination ideal.
    pub fn eliminates(&self, k: usize) -> bool {
        match self {
            Self::Lex => true,
            Self::Block(blocks) => {
                let mut covered = 0;
                for (order, size) in blocks.iter() {
                    if covered == k {
                        return true;
                    }
                    if covered + size > k {
                        return matches!(order, Self::Lex);
                    }
                    covered += size;
                }
                covered == k
            }
            _ => k == 0,
        }
    }

    pub fn compare(&self, a: &[u32], b: &[u32]) -> Ordering {
        match self {
            Self::Lex => lex(a, b),
            Self::GrLex => degree(a).cmp(&degree(b)).then_with(|| lex(a, b)),
            Self::GRevLex => degree(a).cmp(&degree(b)).then_with(|| revlex(a, b)),
            Self::Weighted { weights, tie_break } => weighted(a, weights)
                .cmp(&weighted(b, weights))
                .then_with(|| tie_break.compare(a, b)),
            Self::Block(blocks) => {
                let mut start = 0;
                for (order, size) in blocks.iter() {
                    let end = (start + size).min(a.len());
                    match order.compare(&a[start..end], &b[start..end]) {
                        Ordering::Equal => start = end,
                        other => return other,
                    }
                }
                Ordering::Equal
            }
        }
    }
}

fn degree(a: &[u32]) -> u32 {
    a.iter().sum()
}

fn weighted(a: &[u32], weights: &[u32]) -> u64 {
    a.iter()
        .zip(weights)
        .map(|(e, w)| u64::from(*e) * u64::from(*w))
        .sum()
}

fn lex(a: &[u32], b: &[u32]) -> Ordering {
    a.iter()
        .zip(b)
        .map(|(x, y)| x.cmp(y))
        .find(|o| *o != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
}

fn revlex(a: &[u32], b: &[u32]) -> Ordering {
    a.iter()
        .rev()
        .zip(b.iter().rev())
        .map(|(x, y)| y.cmp(x))
        .find(|o| *o != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
}

/// An exponent vector with a cached total degree.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Monomial {
    pub exponents: Arc<[u32]>,
    degree: u32,
}

impl Monomial {
    pub fn new(exponents: Vec<u32>) -> Self {
        let degree = degree(&exponents);
        Self {
            exponents: exponents.into(),
            degree,
        }
    }
    pub fn one(nvars: usize) -> Self {
        Self::new(vec![0; nvars])
    }
    pub fn variable(index: usize, nvars: usize) -> Self {
        let mut exponents = vec![0; nvars];
        exponents[index] = 1;
        Self::new(exponents)
    }
    pub fn exponents(&self) -> &[u32] {
        &self.exponents
    }
    pub fn nvars(&self) -> usize {
        self.exponents.len()
    }
    pub fn degree(&self) -> u32 {
        self.degree
    }
    pub fn is_one(&self) -> bool {
        self.degree == 0
    }
    pub fn multiply(&self, other: &Self) -> Self {
        debug_assert_eq!(self.nvars(), other.nvars());
        Self::new(self.zip(other, |a, b| a + b))
    }
    pub fn divides(&self, other: &Self) -> bool {
        debug_assert_eq!(self.nvars(), other.nvars());
        self.degree <= other.degree
            && self
                .exponents
                .iter()
                .zip(other.exponents.iter())
                .all(|(a, b)| a <= b)
    }
    pub fn divide(&self, other: &Self) -> Option<Self> {
        other
            .divides(self)
            .then(|| Self::new(self.zip(other, |a, b| a - b)))
    }
    pub fn lcm(&self, other: &Self) -> Self {
        Self::new(self.zip(other, |a, b| a.max(b)))
    }
    pub fn is_coprime(&self, other: &Self) -> bool {
        self.exponents
            .iter()
            .zip(other.exponents.iter())
            .all(|(a, b)| *a == 0 || *b == 0)
    }
    pub fn compare(&self, other: &Self, order: &MonomialOrder) -> Ordering {
        order.compare(&self.exponents, &other.exponents)
    }

    fn zip(&self, other: &Self, f: impl Fn(u32, u32) -> u32) -> Vec<u32> {
        self.exponents
            .iter()
            .zip(other.exponents.iter())
            .map(|(a, b)| f(*a, *b))
            .collect()
    }
}

impl fmt::Display for Monomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_one() {
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
