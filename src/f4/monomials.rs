//! Monomials packed into 7-bit lanes of a `u128`, and the tables that intern them as
//! matrix columns.

use crate::monomial::{Monomial, MonomialOrder};
use crate::par;
use rustc_hash::FxHashMap as HashMap;
use std::borrow::Borrow;
use std::hash::{Hash, Hasher};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, PoisonError};

// Keyed by exponent slice, so a product is looked up before it is allocated.
struct Key(Monomial);
impl PartialEq for Key {
    fn eq(&self, other: &Self) -> bool {
        self.0.exps() == other.0.exps()
    }
}
impl Eq for Key {}
impl Hash for Key {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.0.exps().hash(state);
    }
}
impl Borrow<[u32]> for Key {
    fn borrow(&self) -> &[u32] {
        self.0.exps()
    }
}

// Exponents of up to 16 variables in 7-bit lanes of one integer, so a product is one
// addition that cannot carry, and a lane with its high bit set flags an overflow.
pub(super) const LANES: u128 = 0x8080_8080_8080_8080_8080_8080_8080_8080;

// 2^64 over the golden ratio, the Fibonacci hashing multiplier: the top bits of a
// product with it depend on every bit of the key.
const GOLDEN_GAMMA: u64 = 0x9e37_79b9_7f4a_7c15;
// The first splitmix64 finalizer multiplier, remixing a key already hashed once.
const SPLITMIX_MULTIPLIER: u64 = 0xbf58_476d_1ce4_e5b9;

pub(super) fn pack(exps: &[u32]) -> Option<u128> {
    if exps.len() > 16 {
        return None;
    }
    exps.iter()
        .rev()
        .try_fold(0u128, |acc, &e| (e < 128).then(|| acc << 8 | u128::from(e)))
}

pub(super) fn unpack(k: u128, nvars: usize) -> Monomial {
    let exps: Vec<u32> = (0..nvars).map(|v| (k >> (8 * v)) as u32 & 0xff).collect();
    Monomial::new(exps.as_slice())
}

// Lanewise maximum of packed monomials: the lcm.
pub(super) fn lane_max(a: u128, b: u128) -> u128 {
    let ge = ((a | LANES) - b) & LANES;
    let mask = (ge >> 7) * 0xff;
    (a & mask) | (b & !mask)
}

pub(super) fn lane_divides(a: u128, b: u128) -> bool {
    ((b | LANES) - a) & LANES == LANES
}

pub(super) fn lane_sum(a: u128) -> u32 {
    const BYTES: u128 = 0x00ff_00ff_00ff_00ff_00ff_00ff_00ff_00ff;
    const ONES: u128 = 0x0001_0001_0001_0001_0001_0001_0001_0001;
    let pairs = (a & BYTES) + ((a >> 8) & BYTES);
    (pairs.wrapping_mul(ONES) >> 112) as u32
}

// Ids of packed monomials, probed linearly with each key beside its id, so a lookup
// usually reads one cache line.
#[derive(Default)]
pub(super) struct PackedIds {
    slots: Vec<([u64; 2], u32)>,
    len: usize,
}
impl PackedIds {
    const EMPTY: u32 = u32::MAX;

    // The top bits of the product depend on every key bit, so nearby monomials spread.
    // The length is a power of two; shift by its bit width, not usize's, for 32-bit targets.
    fn start(&self, k: u128) -> usize {
        let h = ((k as u64).wrapping_mul(GOLDEN_GAMMA) ^ (k >> 64) as u64)
            .wrapping_mul(SPLITMIX_MULTIPLIER);
        (h >> (u64::BITS - self.slots.len().trailing_zeros())) as usize
    }

    pub(super) fn get(&self, k: u128) -> Option<u32> {
        if self.slots.is_empty() {
            return None;
        }
        let key = [k as u64, (k >> 64) as u64];
        let mask = self.slots.len() - 1;
        let mut i = self.start(k);
        loop {
            let (s, id) = self.slots[i];
            if id == Self::EMPTY {
                return None;
            }
            if s == key {
                return Some(id);
            }
            i = (i + 1) & mask;
        }
    }

    // `k` must be absent.
    fn insert(&mut self, k: u128, id: u32) {
        if 2 * (self.len + 1) > self.slots.len() {
            let size = (2 * self.slots.len()).max(1024);
            let old = std::mem::replace(&mut self.slots, vec![([0; 2], Self::EMPTY); size]);
            for (s, id) in old.into_iter().filter(|&(_, id)| id != Self::EMPTY) {
                self.place(u128::from(s[0]) | u128::from(s[1]) << 64, id);
            }
        }
        self.place(k, id);
        self.len += 1;
    }

    fn place(&mut self, k: u128, id: u32) {
        let mask = self.slots.len() - 1;
        let mut i = self.start(k);
        while self.slots[i].1 != Self::EMPTY {
            i = (i + 1) & mask;
        }
        self.slots[i] = ([k as u64, (k >> 64) as u64], id);
    }
}

#[derive(Default)]
pub(super) struct MonomialTable {
    pub(super) packed: PackedIds,
    ids: HashMap<Key, u32>,
    pub(super) monomials: Vec<Monomial>,
    scratch: Vec<u32>,
}
impl MonomialTable {
    pub(super) fn intern(&mut self, m: &Monomial) -> (u32, bool) {
        let found = match pack(m.exps()) {
            Some(k) => self.packed.get(k),
            None => self.ids.get(m.exps()).copied(),
        };
        match found {
            Some(id) => (id, false),
            None => (self.insert(m.clone()), true),
        }
    }

    // Packed key of a product, when it fits.
    pub(super) fn product_key(a: &Monomial, b: &Monomial) -> Option<u128> {
        let k = pack(a.exps())? + pack(b.exps())?;
        (k & LANES == 0).then_some(k)
    }

    pub(super) fn intern_product(&mut self, a: &Monomial, b: &Monomial) -> (u32, bool) {
        if let Some(k) = Self::product_key(a, b) {
            if let Some(id) = self.packed.get(k) {
                return (id, false);
            }
            return (self.insert(a * b), true);
        }
        self.scratch.clear();
        self.scratch
            .extend(a.exps().iter().zip(b.exps()).map(|(x, y)| x + y));
        match self.ids.get(self.scratch.as_slice()) {
            Some(&id) => (id, false),
            None => (self.insert(Monomial::new(self.scratch.as_slice())), true),
        }
    }

    // Ids in column order, by descending monomial. Distinct monomials never tie.
    pub(super) fn columns(&self, order: &MonomialOrder) -> Vec<usize> {
        let mut ids: Vec<usize> = (0..self.monomials.len()).collect();
        let cmp = |&a: &usize, &b: &usize| order.compare(&self.monomials[b], &self.monomials[a]);
        par::sort_unstable_by(&mut ids, cmp);
        ids
    }

    #[allow(clippy::expect_used)] // More than 2^32 matrix columns cannot fit in practical memory.
    pub(super) fn insert(&mut self, m: Monomial) -> u32 {
        let id = u32::try_from(self.monomials.len()).expect("F4 matrix exceeds u32 columns");
        match pack(m.exps()) {
            Some(k) => self.packed.insert(k, id),
            None => {
                self.ids.insert(Key(m.clone()), id);
            }
        }
        self.monomials.push(m);
        id
    }
}

// Monomials first met during one level of symbolic preprocessing, interned by
// concurrent workers into shards and numbered after the table's columns.
pub(super) struct Fresh {
    next: AtomicU32,
    shards: Vec<Mutex<HashMap<u128, u32>>>,
}
impl Fresh {
    const SHARD_BITS: u32 = 6;
    const SHARDS: usize = 1 << Self::SHARD_BITS;

    pub(super) fn new(next: usize) -> Self {
        Self {
            next: AtomicU32::new(next as u32),
            shards: (0..Self::SHARDS).map(|_| Mutex::default()).collect(),
        }
    }

    pub(super) fn intern(&self, k: u128) -> u32 {
        let h = ((k >> 64) as u64 ^ k as u64).wrapping_mul(GOLDEN_GAMMA);
        let shard = &self.shards[(h >> (u64::BITS - Self::SHARD_BITS)) as usize];
        let mut map = shard.lock().unwrap_or_else(PoisonError::into_inner);
        *map.entry(k)
            .or_insert_with(|| self.next.fetch_add(1, Ordering::Relaxed))
    }

    // Append the new monomials to the table in id order, returning their ids.
    pub(super) fn drain(self, table: &mut MonomialTable, nvars: usize) -> Vec<u32> {
        let mut new: Vec<(u32, u128)> = self
            .shards
            .into_iter()
            .flat_map(|m| m.into_inner().unwrap_or_else(PoisonError::into_inner))
            .map(|(k, id)| (id, k))
            .collect();
        new.sort_unstable();
        new.into_iter()
            .map(|(_, k)| table.insert(unpack(k, nvars)))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::monomial::{Monomial, divisibility_mask};

    #[test]
    fn mask_never_rejects_a_divisor() {
        for n in [2, 16, 20] {
            for seed in 0..100 {
                let small =
                    Monomial::new((0..n).map(|i| ((seed + i) % 9) as u32).collect::<Vec<_>>());
                let large = Monomial::new(small.exps().iter().map(|e| e + 3).collect::<Vec<_>>());
                assert_eq!(divisibility_mask(&small) & !divisibility_mask(&large), 0);
            }
        }
    }
}
