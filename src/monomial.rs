//! Monomials and orders supplied by polycore.
//!
//! Orders are `Lex`, `GrLex`, `GRevLex`, `MonomialOrder::weighted(weights, tie_break)`,
//! and product orders via `MonomialOrder::block` or `MonomialOrder::elimination(k, rest)`.
pub use polycore::{Monomial, Order as MonomialOrder};

// A necessary divisibility condition: each variable gets four cumulative degree buckets.
// Collisions above sixteen variables only weaken the filter, never reject a divisor.
pub(crate) fn divisibility_mask(m: &Monomial) -> u64 {
    m.exps().iter().enumerate().fold(0, |mask, (i, &e)| {
        let bits = (1u64 << e.min(4)) - 1;
        mask | (bits << ((i % 16) * 4))
    })
}
