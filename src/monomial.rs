//! Monomials and orders supplied by polycore.
//!
//! Orders are `Lex`, `GrLex`, `GRevLex`, `MonomialOrder::weighted(weights, tie_break)`,
//! and product orders via `MonomialOrder::block` or `MonomialOrder::elimination(k, rest)`.
pub use polycore::{Monomial, Order as MonomialOrder};
use std::cmp::Ordering;
/// Compatibility helpers for the 0.3 monomial API.
pub trait MonomialExt {
    fn variable(index: usize, nvars: usize) -> Self;
    fn exponents(&self) -> &[u32];
    fn multiply(&self, other: &Self) -> Self;
    fn divide(&self, other: &Self) -> Option<Self>
    where
        Self: Sized;
    fn compare(&self, other: &Self, order: &MonomialOrder) -> Ordering;
}
impl MonomialExt for Monomial {
    fn variable(index: usize, nvars: usize) -> Self {
        Self::var(index, nvars)
    }
    fn exponents(&self) -> &[u32] {
        self.exps()
    }
    fn multiply(&self, other: &Self) -> Self {
        self * other
    }
    fn divide(&self, other: &Self) -> Option<Self> {
        self.quo(other)
    }
    fn compare(&self, other: &Self, order: &MonomialOrder) -> Ordering {
        order.compare(self, other)
    }
}

// A necessary divisibility condition: each variable gets four cumulative degree buckets.
// Collisions above sixteen variables only weaken the filter, never reject a divisor.
pub(crate) fn divisibility_mask(m: &Monomial) -> u64 {
    m.exps().iter().enumerate().fold(0, |mask, (i, &e)| {
        let bits = (1u64 << e.min(4)) - 1;
        mask | (bits << ((i % 16) * 4))
    })
}
