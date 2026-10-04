//! Sparse row elimination against known pivots: generic field kernels, delayed
//! reduction over word-size residues, and a Monte Carlo echelon form.

use super::{Reducers, SparseRow};
use crate::{Field, ModularField, par};
use polycore::modp::{MulBy, add as add_mod, inv as inv_mod, mul as mul_mod};
use polycore::sample::Rng;
use std::sync::OnceLock;

// A borrowed row, so known pivots can share coefficients.
pub(super) struct Row<'a, F> {
    pub(super) columns: &'a [u32],
    pub(super) coefficients: &'a [F],
}
impl<F> Clone for Row<'_, F> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<F> Copy for Row<'_, F> {}
impl<F> SparseRow<F> {
    fn row(&self) -> Row<'_, F> {
        Row {
            columns: &self.columns,
            coefficients: &self.coefficients,
        }
    }
}

/// Pivot rows indexed by leading column.
pub(super) trait Pivots<F>: Sync {
    fn get(&self, col: u32) -> Option<Row<'_, F>>;
}

pub(super) struct Table<'a, F> {
    rows: Vec<Row<'a, F>>,
    index: Vec<u32>,
}
impl<'a, F> Table<'a, F> {
    pub(super) fn new(rows: impl IntoIterator<Item = Row<'a, F>>, ncols: usize) -> Self {
        let rows: Vec<_> = rows.into_iter().collect();
        let mut index = vec![u32::MAX; ncols];
        for (k, row) in rows.iter().enumerate() {
            if let Some(&c) = row.columns.first() {
                index[c as usize] = k as u32;
            }
        }
        Self { rows, index }
    }
}
impl<F: Sync> Pivots<F> for Table<'_, F> {
    fn get(&self, col: u32) -> Option<Row<'_, F>> {
        self.rows.get(self.index[col as usize] as usize).copied()
    }
}

// Known pivots plus one slot per column for rows claimed during elimination.
struct Claimed<'a, F> {
    known: &'a Table<'a, F>,
    slots: Vec<OnceLock<(usize, SparseRow<F>)>>,
}
impl<F: Send + Sync> Pivots<F> for Claimed<'_, F> {
    fn get(&self, col: u32) -> Option<Row<'_, F>> {
        self.known
            .get(col)
            .or_else(|| self.slots[col as usize].get().map(|(_, r)| r.row()))
    }
}

/// Row reduction and normalization for one coefficient representation.
pub(super) trait Kernel<F>: Sync {
    type Scratch;
    fn scratch(&self) -> Self::Scratch;
    fn reduce(
        &self,
        row: &SparseRow<F>,
        pivots: &impl Pivots<F>,
        scratch: &mut Self::Scratch,
    ) -> SparseRow<F>;
    fn normalize(&self, row: &mut SparseRow<F>);
    // Rows per parallel block when mapping `n` rows.
    fn block(&self, _n: usize) -> usize {
        32
    }
}

struct Generic(usize);
impl<F: Field + Send + Sync> Kernel<F> for Generic {
    type Scratch = ();
    fn scratch(&self) {}
    fn reduce(&self, row: &SparseRow<F>, pivots: &impl Pivots<F>, (): &mut ()) -> SparseRow<F> {
        reduce_generic(row, pivots, self.0)
    }
    fn normalize(&self, row: &mut SparseRow<F>) {
        let inv = row.coefficients[0]
            .inverse()
            .unwrap_or_else(|| unreachable!("nonzero coefficient"));
        for v in &mut row.coefficients {
            *v = v.clone() * inv.clone();
        }
    }
}

/// Stored residues in `[0, p)`: `u32` whenever the prime fits, halving row traffic.
trait Residue: Copy + Send + Sync + Into<u64> {
    fn new(v: u64) -> Self;
}
impl Residue for u16 {
    fn new(v: u64) -> Self {
        v as u16
    }
}
impl Residue for u32 {
    fn new(v: u64) -> Self {
        v as u32
    }
}
impl Residue for u64 {
    fn new(v: u64) -> Self {
        v
    }
}

struct Dense {
    p: u64,
    ncols: usize,
}
impl<C: Residue> Kernel<C> for Dense {
    type Scratch = Vec<u64>;
    fn scratch(&self) -> Vec<u64> {
        vec![0; self.ncols]
    }
    fn reduce(
        &self,
        row: &SparseRow<C>,
        pivots: &impl Pivots<C>,
        buf: &mut Vec<u64>,
    ) -> SparseRow<C> {
        reduce_dense(row, pivots, buf, self.p)
    }
    fn normalize(&self, row: &mut SparseRow<C>) {
        normalize_residues(&mut row.coefficients, self.p);
    }
}

// Reduce `rows` by the monic `pivots` and by each other, then back substitute. Workers
// claim free leading columns concurrently, and a row that loses a claim keeps reducing
// by the winner. The reduced echelon form is unique, so the result does not depend on
// scheduling, though which rows are reported live (independent) may.
pub(super) fn echelon<F: Clone + Send + Sync, K: Kernel<F>>(
    kernel: &K,
    known: &Table<'_, F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    let claimed = Claimed {
        known,
        slots: (0..ncols).map(|_| OnceLock::new()).collect(),
    };
    let order = sorted(rows);
    par::map_blocks(
        &order,
        kernel.block(order.len()),
        || kernel.scratch(),
        |&i, s| {
            let r = kernel.reduce(&rows[i], &claimed, s);
            claim(kernel, &claimed, i, r, s);
        },
    );
    back_substitute(kernel, claimed, ncols)
}

// Nonempty row indices by leading column, then length.
fn sorted<F>(rows: &[SparseRow<F>]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..rows.len())
        .filter(|&i| !rows[i].columns.is_empty())
        .collect();
    order.sort_by_key(|&i| (rows[i].columns[0], rows[i].columns.len()));
    order
}

// Claim the leading column of a reduced row, reducing by the winner after each lost
// race. Returns whether a new pivot was stored.
fn claim<F: Send + Sync, K: Kernel<F>>(
    kernel: &K,
    claimed: &Claimed<'_, F>,
    i: usize,
    mut r: SparseRow<F>,
    s: &mut K::Scratch,
) -> bool {
    while let Some(&lead) = r.columns.first() {
        kernel.normalize(&mut r);
        match claimed.slots[lead as usize].set((i, r)) {
            Ok(()) => return true,
            Err((_, lost)) => r = kernel.reduce(&lost, claimed, s),
        }
    }
    false
}

fn back_substitute<F: Clone + Send + Sync, K: Kernel<F>>(
    kernel: &K,
    claimed: Claimed<'_, F>,
    ncols: usize,
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    let (live, new): (Vec<usize>, Vec<SparseRow<F>>) = claimed
        .slots
        .into_iter()
        .filter_map(OnceLock::into_inner)
        .unzip();
    // The column sweep clears every pivot column it meets, fill-in included, so each
    // tail reduces independently against the unreduced echelon form.
    let table = Table::new(new.iter().map(SparseRow::row), ncols);
    let reduced = par::map_blocks(
        &new,
        kernel.block(new.len()),
        || kernel.scratch(),
        |r, s| {
            if !r.columns[1..].iter().any(|&c| table.get(c).is_some()) {
                return r.clone();
            }
            let tail = SparseRow {
                columns: r.columns[1..].to_vec(),
                coefficients: r.coefficients[1..].to_vec(),
            };
            let mut out = kernel.reduce(&tail, &table, s);
            out.columns.insert(0, r.columns[0]);
            out.coefficients.insert(0, r.coefficients[0].clone());
            out
        },
    );
    (reduced, live)
}

pub(super) fn echelonize_generic<F: Field + Send + Sync>(
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    let known = Table::new(pivots.view(&pivots.coefficients), ncols);
    echelon(&Generic(ncols), &known, rows, ncols).0
}

// Span of a row's columns, as buffer indices.
pub(super) fn span(columns: &[u32]) -> (usize, usize) {
    match (columns.first(), columns.last()) {
        (Some(&lo), Some(&hi)) => (lo as usize, hi as usize),
        _ => (usize::MAX, 0),
    }
}

pub(super) fn reduce_generic<F: Field>(
    row: &SparseRow<F>,
    pivots: &impl Pivots<F>,
    ncols: usize,
) -> SparseRow<F> {
    let mut buf: Vec<F> = vec![F::zero(); ncols];
    for (&c, v) in row.columns.iter().zip(&row.coefficients) {
        buf[c as usize] = v.clone();
    }
    let (lo, mut hi) = span(&row.columns);
    let mut j = lo;
    while j <= hi && j < ncols {
        if !buf[j].is_zero()
            && let Some(piv) = pivots.get(j as u32)
        {
            let c = std::mem::replace(&mut buf[j], F::zero());
            for (&pc, pv) in piv.columns.iter().zip(piv.coefficients).skip(1) {
                let pc = pc as usize;
                buf[pc] = buf[pc].clone() - c.clone() * pv.clone();
            }
            hi = hi.max(span(piv.columns).1);
        }
        j += 1;
    }
    compress(buf, lo)
}

fn compress<F: Field>(buf: Vec<F>, lo: usize) -> SparseRow<F> {
    let mut row = SparseRow {
        columns: Vec::new(),
        coefficients: Vec::new(),
    };
    for (c, v) in buf.into_iter().enumerate().skip(lo) {
        if !v.is_zero() {
            row.columns.push(c as u32);
            row.coefficients.push(v);
        }
    }
    row
}

pub(super) fn echelonize_modular<F: ModularField + Send + Sync>(
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    // Tiny fields leave too few multipliers for random combinations to be reliable.
    match modulus(pivots, rows) {
        None => echelonize_generic(pivots, rows, ncols),
        Some(p) if p < 1 << 12 => echelonize_modular_traced(pivots, rows, ncols).0,
        Some(p) if p <= u64::from(u16::MAX) => random_residues::<F, u16>(p, pivots, rows, ncols),
        Some(p) if p <= u64::from(u32::MAX) => random_residues::<F, u32>(p, pivots, rows, ncols),
        Some(p) => random_residues::<F, u64>(p, pivots, rows, ncols),
    }
}

fn random_residues<F: ModularField + Send + Sync, C: Residue>(
    p: u64,
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    let (known, pending) = residues::<F, C>(p, pivots, rows);
    let known = Table::new(pivots.view(&known), ncols);
    lift(echelon_random(p, &known, &pending, ncols), p)
}

// Monte Carlo echelon form. Sorted rows are split into about sqrt(n / 3) blocks, and
// each block is replaced by random linear combinations of its rows, reduced and
// claimed like ordinary rows, until `zeros` consecutive combinations vanish. A block
// that still has an independent row escapes with probability at most p^-zeros, so
// only about rank + blocks dense reductions are needed instead of one per row.
fn echelon_random<C: Residue>(
    p: u64,
    known: &Table<'_, C>,
    rows: &[SparseRow<C>],
    ncols: usize,
) -> Vec<SparseRow<C>> {
    let kernel = Dense { p, ncols };
    let claimed = Claimed {
        known,
        slots: (0..ncols).map(|_| OnceLock::new()).collect(),
    };
    let order = sorted(rows);
    let nblocks = (order.len() as f64 / 3.0).sqrt() as usize + 1;
    let blocks: Vec<&[usize]> = order.chunks(order.len().div_ceil(nblocks).max(1)).collect();
    let zeros = 40u32.div_ceil(p.ilog2()).max(2) as usize;
    let block = |(buf, lanes): &mut (Vec<u64>, Vec<[u64; COMBOS]>), b: usize, block: &&[usize]| {
        let mut rng = Rng::new(b as u64);
        let (mut found, mut run) = (0, 0);
        // Several combinations share each pass over the pivots when every cell can defer
        // its reductions, and are claimed in turn as if swept one after another.
        if p < (1 << 31) && fits(ncols + block.len(), p) {
            lanes.resize(ncols, [0; COMBOS]);
            while found < block.len() && run < zeros {
                let (mut lo, mut hi) = (usize::MAX, 0);
                for &i in *block {
                    let m = std::array::from_fn(|_| rng.nonzero(p) as u32);
                    let (cols, coefs) = (&rows[i].columns, &rows[i].coefficients);
                    accumulate_combos(lanes, cols, coefs, m);
                    let (l, h) = span(cols);
                    (lo, hi) = (lo.min(l), hi.max(h));
                }
                for r in sweep_combos(lanes, (lo, hi), &claimed, p) {
                    if found == block.len() || run == zeros {
                        break;
                    }
                    if claim(&kernel, &claimed, b, r, buf) {
                        (found, run) = (found + 1, 0);
                    } else {
                        run += 1;
                    }
                }
            }
            return;
        }
        // Unreduced sums only when the sweep can also defer every reduction.
        let lazy = p < (1 << 31) && fits(ncols + block.len(), p);
        while found < block.len() && run < zeros {
            let (mut lo, mut hi) = (usize::MAX, 0);
            for &i in *block {
                let m = rng.nonzero(p);
                let by = MulBy::new(m, p);
                let (cols, coefs) = (&rows[i].columns, &rows[i].coefficients);
                for (&c, &v) in cols.iter().zip(coefs) {
                    let cell = &mut buf[c as usize];
                    *cell = if lazy {
                        *cell + m * v.into()
                    } else {
                        add_mod(*cell, by.mul(v.into(), p), p)
                    };
                }
                let (l, h) = span(cols);
                (lo, hi) = (lo.min(l), hi.max(h));
            }
            let terms = if lazy { block.len() } else { 0 };
            let r = sweep(buf, (lo, hi), terms, &claimed, p);
            if claim(&kernel, &claimed, b, r, buf) {
                (found, run) = (found + 1, 0);
            } else {
                run += 1;
            }
        }
    };
    let scratch = || (Kernel::<C>::scratch(&kernel), Vec::new());
    par::for_each_init(&blocks, scratch, block);
    back_substitute(&kernel, claimed, ncols).0
}

// Whether `terms` products of residues plus one residue fit in 64 bits.
pub(super) fn fits(terms: usize, p: u64) -> bool {
    (terms as u128) * u128::from(p - 1).pow(2) + u128::from(p - 1) <= u128::from(u64::MAX)
}

pub(super) fn echelonize_modular_traced<F: ModularField + Send + Sync>(
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    match modulus(pivots, rows) {
        None => (
            echelonize_generic(pivots, rows, ncols),
            (0..rows.len()).collect(),
        ),
        Some(p) if p <= u64::from(u16::MAX) => echelon_residues::<F, u16>(p, pivots, rows, ncols),
        Some(p) if p <= u64::from(u32::MAX) => echelon_residues::<F, u32>(p, pivots, rows, ncols),
        Some(p) => echelon_residues::<F, u64>(p, pivots, rows, ncols),
    }
}

fn echelon_residues<F: ModularField + Send + Sync, C: Residue>(
    p: u64,
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    let (known, pending) = residues::<F, C>(p, pivots, rows);
    let known = Table::new(pivots.view(&known), ncols);
    let (reduced, live) = echelon(&Dense { p, ncols }, &known, &pending, ncols);
    (lift(reduced, p), live)
}

pub(super) fn reduce_rows_modular<F: ModularField + Send + Sync>(
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    match modulus(pivots, rows) {
        None => rows.to_vec(),
        Some(p) if p <= u64::from(u16::MAX) => reduce_residues::<F, u16>(p, pivots, rows, ncols),
        Some(p) if p <= u64::from(u32::MAX) => reduce_residues::<F, u32>(p, pivots, rows, ncols),
        Some(p) => reduce_residues::<F, u64>(p, pivots, rows, ncols),
    }
}

fn reduce_residues<F: ModularField + Send + Sync, C: Residue>(
    p: u64,
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    let (known, pending) = residues::<F, C>(p, pivots, rows);
    let kernel = Dense { p, ncols };
    let table = Table::new(pivots.view(&known), ncols);
    let reduced = par::map_rows(
        &pending,
        || Kernel::<C>::scratch(&kernel),
        |r, s| kernel.reduce(r, &table, s),
    );
    lift(reduced, p)
}

// The modulus, or `None` when every coefficient is zero and so none can be recovered.
fn modulus<F: ModularField>(pivots: &Reducers<F>, rows: &[SparseRow<F>]) -> Option<u64> {
    pivots
        .coefficients
        .iter()
        .flatten()
        .chain(rows.iter().flat_map(|r| &r.coefficients))
        .map(ModularField::modulus)
        .find(|&m| m != 0)
}

type Residues<C> = (Vec<Vec<C>>, Vec<SparseRow<C>>);

// Residues of monic-normalized pivot coefficients and of rows.
fn residues<F: ModularField, C: Residue>(
    p: u64,
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
) -> Residues<C> {
    let convert = |cs: &[F]| -> Vec<C> { cs.iter().map(|v| C::new(v.residue_mod(p))).collect() };
    let mut known: Vec<_> = pivots.coefficients.iter().map(|cs| convert(cs)).collect();
    for cs in &mut known {
        normalize_residues(cs, p);
    }
    let rows = rows
        .iter()
        .map(|r| SparseRow {
            columns: r.columns.clone(),
            coefficients: convert(&r.coefficients),
        })
        .collect();
    (known, rows)
}

fn lift<F: ModularField, C: Residue>(rows: Vec<SparseRow<C>>, p: u64) -> Vec<SparseRow<F>> {
    rows.into_iter()
        .map(|r| SparseRow {
            columns: r.columns,
            coefficients: r
                .coefficients
                .into_iter()
                .map(|v| F::from_residue(v.into(), p))
                .collect(),
        })
        .collect()
}

fn normalize_residues<C: Residue>(coefficients: &mut [C], p: u64) {
    let Some(&lead) = coefficients.first() else {
        return;
    };
    if lead.into() == 1 {
        return;
    }
    let inv = inv_mod(lead.into(), p);
    for v in coefficients {
        *v = C::new(mul_mod((*v).into(), inv, p));
    }
}

// Monagan-Pearce dense reduction: sums are reduced lazily while they fit in 64 bits.
fn reduce_dense<C: Residue>(
    row: &SparseRow<C>,
    pivots: &impl Pivots<C>,
    buf: &mut [u64],
    p: u64,
) -> SparseRow<C> {
    // `buf` is all zero on entry and is cleared again before returning.
    for (&c, &v) in row.columns.iter().zip(&row.coefficients) {
        buf[c as usize] = v.into();
    }
    let (lo, hi) = span(&row.columns);
    sweep(buf, (lo, hi), 0, pivots, p)
}

// Random combinations of a block reduced together in `echelon_random`.
const COMBOS: usize = 4;

// `sweep` of `COMBOS` rows at once with deferred reductions, reading each pivot once.
// Returns the reduced rows, possibly empty, and clears the cells it read.
fn sweep_combos<C: Residue>(
    buf: &mut [[u64; COMBOS]],
    (lo, mut hi): (usize, usize),
    pivots: &impl Pivots<C>,
    p: u64,
) -> [SparseRow<C>; COMBOS] {
    let mut out: [SparseRow<C>; COMBOS] = std::array::from_fn(|_| SparseRow {
        columns: Vec::new(),
        coefficients: Vec::new(),
    });
    let mut j = lo;
    while j <= hi && j < buf.len() {
        if buf[j] == [0; COMBOS] {
            j += 1;
            continue;
        }
        let v = std::mem::take(&mut buf[j]).map(|x| x % p);
        let Some(piv) = pivots.get(j as u32) else {
            for (row, &x) in out.iter_mut().zip(&v) {
                if x != 0 {
                    row.columns.push(j as u32);
                    row.coefficients.push(C::new(x));
                }
            }
            j += 1;
            continue;
        };
        if v != [0; COMBOS] {
            let c = v.map(|x| ((p - x) % p) as u32);
            hi = hi.max(span(piv.columns).1);
            accumulate_combos(buf, &piv.columns[1..], &piv.coefficients[1..], c);
        }
        j += 1;
    }
    out
}

// `buf[column] += c * coefficient` in every lane, without reduction, for `p < 2^31`.
#[cfg(not(target_arch = "aarch64"))]
fn accumulate_combos<C: Residue>(
    buf: &mut [[u64; COMBOS]],
    columns: &[u32],
    coefficients: &[C],
    c: [u32; COMBOS],
) {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 support was just detected.
        #[allow(unsafe_code)]
        return unsafe { accumulate_combos_avx2(buf, columns, coefficients, c) };
    }
    for (&pc, &pv) in columns.iter().zip(coefficients) {
        let cell = &mut buf[pc as usize];
        let pv: u64 = pv.into();
        for l in 0..COMBOS {
            cell[l] += u64::from(c[l]) * pv;
        }
    }
}

// One widening multiply by the broadcast coefficient and an add cover the four lanes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(unsafe_code)]
unsafe fn accumulate_combos_avx2<C: Residue>(
    buf: &mut [[u64; COMBOS]],
    columns: &[u32],
    coefficients: &[C],
    c: [u32; COMBOS],
) {
    use std::arch::x86_64::{
        __m128i, __m256i, _mm_loadu_si128, _mm256_add_epi64, _mm256_cvtepu32_epi64,
        _mm256_loadu_si256, _mm256_mul_epu32, _mm256_set1_epi64x, _mm256_storeu_si256,
    };
    const { assert!(COMBOS == 4) };
    // SAFETY: each pointer addresses a whole array of four `u32` or four `u64` lanes.
    unsafe {
        let m = _mm256_cvtepu32_epi64(_mm_loadu_si128(c.as_ptr().cast::<__m128i>()));
        for (&pc, &pv) in columns.iter().zip(coefficients) {
            let cell = buf[pc as usize].as_mut_ptr().cast::<__m256i>();
            let v = _mm256_set1_epi64x(pv.into() as i64);
            let acc = _mm256_add_epi64(_mm256_loadu_si256(cell), _mm256_mul_epu32(v, m));
            _mm256_storeu_si256(cell, acc);
        }
    }
}

// Two widening multiply-accumulates by the scalar coefficient cover the four lanes.
#[cfg(target_arch = "aarch64")]
#[allow(unsafe_code)]
fn accumulate_combos<C: Residue>(
    buf: &mut [[u64; COMBOS]],
    columns: &[u32],
    coefficients: &[C],
    c: [u32; COMBOS],
) {
    use std::arch::aarch64::{
        vget_low_u32, vld1q_u32, vld1q_u64, vmlal_high_n_u32, vmlal_n_u32, vst1q_u64,
    };
    const { assert!(COMBOS == 4) };
    // SAFETY: NEON is part of the aarch64 baseline, and each pointer addresses a whole
    // array of four `u32` or four `u64` lanes.
    unsafe {
        let m = vld1q_u32(c.as_ptr());
        for (&pc, &pv) in columns.iter().zip(coefficients) {
            let cell = buf[pc as usize].as_mut_ptr();
            let pv = pv.into() as u32;
            let lo = vmlal_n_u32(vld1q_u64(cell), vget_low_u32(m), pv);
            let hi = vmlal_high_n_u32(vld1q_u64(cell.add(2)), m, pv);
            vst1q_u64(cell, lo);
            vst1q_u64(cell.add(2), hi);
        }
    }
}

// Primes below which a cell in [0, 4p) fits in 64 bits, for Shoup reduction in `sweep`.
const SHOUP_LIMIT: u64 = 1 << 62;

// Eliminate pivot columns from `buf[lo..=hi]` and drain it into a row. Each cell holds
// a residue plus at most `terms` unreduced products of two residues.
fn sweep<C: Residue>(
    buf: &mut [u64],
    (lo, mut hi): (usize, usize),
    terms: usize,
    pivots: &impl Pivots<C>,
    p: u64,
) -> SparseRow<C> {
    let small = p < (1 << 31);
    let square = if small { p * p } else { 0 };
    let ncols = buf.len();
    // At most ncols triangular pivot eliminations can contribute to any cell.
    let deferred = small && fits(ncols + terms, p);
    // Cells behind the sweep are final, so leftovers are reduced in place and counted,
    // and a row that reduces to zero needs no second pass.
    let (mut first, mut left) = (usize::MAX, 0);
    let mut j = lo;
    while j <= hi && j < ncols {
        if buf[j] == 0 {
            j += 1;
            continue;
        }
        let Some(piv) = pivots.get(j as u32) else {
            buf[j] %= p;
            if buf[j] != 0 {
                first = first.min(j);
                left += 1;
            }
            j += 1;
            continue;
        };
        let v = std::mem::take(&mut buf[j]) % p;
        if v != 0 {
            let c = p - v;
            hi = hi.max(span(piv.columns).1);
            let (cols, coefs) = (&piv.columns[1..], &piv.coefficients[1..]);
            let tail = cols.iter().zip(coefs);
            if deferred {
                // Fixed-width chunks let the independent scatter updates overlap.
                let (kc, kr) = cols.as_chunks::<8>();
                let (vc, vr) = coefs.as_chunks::<8>();
                for (ks, vs) in kc.iter().zip(vc) {
                    for u in 0..8 {
                        buf[ks[u] as usize] += c * vs[u].into();
                    }
                }
                for (&pc, &pv) in kr.iter().zip(vr) {
                    buf[pc as usize] += c * pv.into();
                }
            } else if small {
                // Cells stay below p^2 < 2^62, so one product never overflows, and
                // the conditional subtraction compiles to a select rather than a branch.
                let (kc, kr) = cols.as_chunks::<8>();
                let (vc, vr) = coefs.as_chunks::<8>();
                let update = |cell: &mut u64, v: C| {
                    let acc = *cell + c * v.into();
                    *cell = if acc >= square { acc - square } else { acc };
                };
                for (ks, vs) in kc.iter().zip(vc) {
                    for u in 0..8 {
                        update(&mut buf[ks[u] as usize], vs[u]);
                    }
                }
                for (&pc, &pv) in kr.iter().zip(vr) {
                    update(&mut buf[pc as usize], pv);
                }
            } else if p < SHOUP_LIMIT {
                // Shoup's precomputed quotient makes `c * v - q * p` exact in [0, 2p)
                // without a division per cell, and cells stay below 2p until read.
                let cq = ((u128::from(c) << 64) / u128::from(p)) as u64;
                let twice = 2 * p;
                for (&pc, &pv) in tail {
                    let v: u64 = pv.into();
                    let q = ((u128::from(cq) * u128::from(v)) >> 64) as u64;
                    let cell = &mut buf[pc as usize];
                    let acc = *cell + c.wrapping_mul(v).wrapping_sub(q.wrapping_mul(p));
                    *cell = if acc >= twice { acc - twice } else { acc };
                }
            } else {
                let by = MulBy::new(c, p);
                for (&pc, &pv) in tail {
                    buf[pc as usize] = add_mod(buf[pc as usize], by.mul(pv.into(), p), p);
                }
            }
        }
        j += 1;
    }
    let mut out = SparseRow {
        columns: Vec::with_capacity(left),
        coefficients: Vec::with_capacity(left),
    };
    for (c, cell) in buf.iter_mut().enumerate().skip(first) {
        if out.columns.len() == left {
            break;
        }
        let v = std::mem::take(cell);
        if v != 0 {
            out.columns.push(c as u32);
            out.coefficients.push(C::new(v));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_field::Zp;

    // Knuth's MMIX linear congruential multiplier, for a reproducible coefficient stream.
    const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;

    #[test]
    fn modular_row_blocks_match_generic_elimination_at_word_boundaries() {
        for p in [
            2,
            32003,
            2_147_483_647,
            2_147_483_659,
            4_611_686_018_427_387_847,
            18_446_744_073_709_551_557,
        ] {
            let ncols = 80;
            let mut seed = 42u64;
            let mut next = || {
                seed = seed.wrapping_mul(LCG_MULTIPLIER).wrapping_add(1);
                Zp::new(seed, p)
            };
            let pivots: Vec<_> = (0..ncols)
                .step_by(3)
                .map(|i| {
                    let columns: Vec<_> = (i as u32..ncols as u32).collect();
                    let mut coefficients: Vec<_> = columns.iter().map(|_| next()).collect();
                    coefficients[0] = Zp::new(1, p);
                    SparseRow {
                        columns,
                        coefficients,
                    }
                })
                .collect();
            let rows: Vec<_> = (0..96)
                .map(|_| SparseRow {
                    columns: (0..ncols as u32).collect(),
                    coefficients: (0..ncols).map(|_| next()).collect(),
                })
                .collect();
            let pivots = Reducers::from(pivots);
            let expected = echelonize_generic(&pivots, &rows, ncols);
            let actual = echelonize_modular(&pivots, &rows, ncols);
            assert_eq!(actual.len(), expected.len());
            for (a, b) in actual.iter().zip(&expected) {
                assert_eq!(a.columns, b.columns, "modulus {p}");
                assert_eq!(a.coefficients, b.coefficients, "modulus {p}");
            }
        }
    }
}
