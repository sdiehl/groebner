//! F4 Groebner basis computation (Faugere) with Gebauer-Moller pair selection and dense
//! delayed-reduction linear algebra over machine-word prime fields.
//!
//! ```
//! use groebner::{groebner_basis_f4, MonomialOrder, PolynomialRing, PrimeField};
//!
//! type F32003 = PrimeField<32003>;
//!
//! let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::Lex)?;
//! let f1 = ring.parse("x^2 - y")?;
//! let f2 = ring.parse("x*y - 1")?;
//! let basis = groebner_basis_f4(vec![f1, f2], true)?;
//!
//! assert!(!basis.is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::PolynomialExt;
use crate::finite_field::{PrimeField, Zp};
use crate::groebner::{GroebnerError, finish_basis, prepare_input};
use crate::monomial::{Monomial, MonomialOrder, divisibility_mask};
use crate::polynomial::Polynomial;
use crate::polynomial::term;
use crate::{Field, ModularField};
use num_rational::BigRational;
use polycore::modp::{add as add_mod, inv as inv_mod, mul as mul_mod};
use polycore::sample::Rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::borrow::{Borrow, Cow};
use std::hash::{Hash, Hasher};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

/// A sparse matrix row with columns in ascending index order (descending monomial order).
#[derive(Debug, Clone)]
pub struct SparseRow<F> {
    pub columns: Vec<u32>,
    pub coefficients: Vec<F>,
}

/// Known pivot rows. Each is a monomial multiple of one polynomial, so all multiples share
/// that polynomial's coefficients and differ only in their columns.
#[derive(Debug, Clone)]
pub struct Reducers<F> {
    /// Coefficients of each polynomial, leading coefficient first.
    pub coefficients: Vec<Vec<F>>,
    /// Each row's index into `coefficients`, and its columns.
    pub rows: Vec<(usize, Vec<u32>)>,
}

impl<F> From<Vec<SparseRow<F>>> for Reducers<F> {
    fn from(rows: Vec<SparseRow<F>>) -> Self {
        let (rows, coefficients) = rows
            .into_iter()
            .enumerate()
            .map(|(k, r)| ((k, r.columns), r.coefficients))
            .unzip();
        Self { coefficients, rows }
    }
}

impl<F> Reducers<F> {
    // Rows over `coefficients`, either `self.coefficients` or a converted copy of them.
    fn view<'a, C>(&'a self, coefficients: &'a [Vec<C>]) -> impl Iterator<Item = Row<'a, C>> {
        self.rows.iter().map(|(k, columns)| Row {
            columns,
            coefficients: &coefficients[*k],
        })
    }
}

// Known pivots as multiples of `n` basis elements, given their coefficients, with the
// given columns.
fn reducers<F>(
    n: usize,
    coefficients_of: impl Fn(usize) -> Vec<F>,
    rows: impl IntoIterator<Item = (usize, Vec<u32>)>,
) -> Reducers<F> {
    let mut shared = vec![usize::MAX; n];
    let mut coefficients = Vec::new();
    let rows = rows
        .into_iter()
        .map(|(i, columns)| {
            if shared[i] == usize::MAX {
                shared[i] = coefficients.len();
                coefficients.push(coefficients_of(i));
            }
            (shared[i], columns)
        })
        .collect();
    Reducers { coefficients, rows }
}

/// Coefficient fields usable by [`groebner_basis_f4`]. The default method is a generic dense
/// eliminator; prime fields override it with a delayed-reduction fast path.
pub trait F4Field: Field + Send + Sync {
    /// Optional multi-modular backend. Return `None` to use direct F4.
    fn modular(
        _polynomials: &[Polynomial<Self>],
        _canonicalize: bool,
    ) -> Option<Result<Vec<Polynomial<Self>>, GroebnerError>> {
        None
    }

    /// Echelon form and indices of independent input rows, used when learning a trace.
    fn echelonize_traced(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> (Vec<SparseRow<Self>>, Vec<usize>) {
        (
            Self::echelonize(pivots, rows, ncols),
            (0..rows.len()).collect(),
        )
    }

    /// Reduce `rows` modulo the monic `pivots` (distinct leading columns), then echelonize and
    /// interreduce the remainders, returning the nonzero monic rows.
    fn echelonize(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_generic(pivots, rows, ncols)
    }

    /// Fully reduce each row modulo the monic `pivots` (distinct leading columns).
    fn reduce_rows(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        let table = Table::new(pivots.view(&pivots.coefficients), ncols);
        rows.iter()
            .map(|r| reduce_generic(r, &table, ncols))
            .collect()
    }
}

impl F4Field for BigRational {
    fn modular(
        polynomials: &[Polynomial<Self>],
        _canonicalize: bool,
    ) -> Option<Result<Vec<Polynomial<Self>>, GroebnerError>> {
        Some(crate::modular::groebner_basis_f4_rational(
            polynomials.to_vec(),
            false,
        ))
    }
}

impl<const P: u64> F4Field for PrimeField<P> {
    fn echelonize_traced(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> (Vec<SparseRow<Self>>, Vec<usize>) {
        echelonize_modular_traced(pivots, rows, ncols)
    }

    fn echelonize(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_modular(pivots, rows, ncols)
    }

    fn reduce_rows(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        reduce_rows_modular(pivots, rows, ncols)
    }
}

impl F4Field for Zp {
    fn echelonize_traced(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> (Vec<SparseRow<Self>>, Vec<usize>) {
        echelonize_modular_traced(pivots, rows, ncols)
    }

    fn echelonize(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_modular(pivots, rows, ncols)
    }

    fn reduce_rows(
        pivots: &Reducers<Self>,
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        reduce_rows_modular(pivots, rows, ncols)
    }
}

/// Compute a Groebner basis with F4 under the polynomials' own order.
///
/// Over `BigRational`, uses multi-modular reconstruction and returns a reduced basis
/// even when `canonicalize` is false. Acceptance uses a fresh prime; see
/// [`crate::groebner_basis_f4_rational`] for optional exact checks.
pub fn groebner_basis_f4<F: F4Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    if let Some(result) = F::modular(&polynomials, canonicalize) {
        return result;
    }
    groebner_basis_f4_direct(polynomials, canonicalize)
}

/// Direct F4, bypassing modular reconstruction (useful for exact comparisons).
pub fn groebner_basis_f4_direct<F: F4Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    compute(polynomials, canonicalize, None)
}

// Runs on a pool thread, so parallel steps split work by stealing instead of each
// waking the pool from outside and sleeping until it finishes.
fn compute<F: F4Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
    trace: Option<&mut F4Trace>,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    #[cfg(feature = "parallel")]
    if rayon::current_thread_index().is_none() {
        return rayon::scope(|_| compute_on(polynomials, canonicalize, trace));
    }
    compute_on(polynomials, canonicalize, trace)
}

fn compute_on<F: F4Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
    mut trace: Option<&mut F4Trace>,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    let input = prepare_input(polynomials)?;
    let order = input[0].order.clone();
    let nvars = input[0].nvars;
    let traced = trace.is_some();
    let mut state = State {
        basis: Vec::new(),
        active: Vec::new(),
        by_size: Vec::new(),
        pairs: Vec::new(),
        masks: Vec::new(),
        keys: Vec::new(),
        nvars,
    };
    for poly in input {
        state.update(Element::new(&poly, traced));
    }
    if let Some(t) = trace.as_deref_mut() {
        t.input_leads = state.basis.iter().map(|g| g.lm.clone()).collect();
    }
    while !state.pairs.is_empty() {
        let selected = state.select();
        let mut plan = state.symbolic_preprocessing(&selected, &order, traced);
        let (new_rows, live) = {
            let matrix = state.matrix(&mut plan, traced);
            if traced {
                F::echelonize_traced(&matrix.pivots, &matrix.rows, plan.columns.len())
            } else {
                (
                    F::echelonize(&matrix.pivots, &matrix.rows, plan.columns.len()),
                    Vec::new(),
                )
            }
        };
        let keys: Vec<Option<u128>> = plan.columns.iter().map(|m| pack(m.exps())).collect();
        let mut new: Vec<Element<F>> = new_rows
            .into_iter()
            .filter(|row| !row.columns.is_empty())
            .map(|row| Element::decode(row, &plan.columns, &keys, traced))
            .collect();
        new.sort_by(|a, b| order.compare(&a.lm, &b.lm));
        if let Some(t) = trace.as_deref_mut() {
            let leads = new.iter().map(|g| g.lm.clone()).collect();
            t.rounds.push(TraceRound { plan, leads, live });
        }
        for g in new {
            state.update(g);
        }
    }
    if let Some(t) = trace {
        t.active = state.active.clone();
    }
    // Terms share one monomial per distinct exponent vector, as decoded rows would.
    let mut cache = HashMap::default();
    let basis = state
        .basis
        .into_iter()
        .zip(state.active)
        .filter_map(|(g, active)| active.then(|| g.polynomial(nvars, &order, &mut cache)))
        .collect();
    finish(basis, canonicalize)
}

/// Coefficient-independent F4 plans, learned at a prime and reusable at other primes.
#[derive(Default)]
pub(crate) struct F4Trace {
    input_leads: Vec<Monomial>,
    rounds: Vec<TraceRound>,
    active: Vec<bool>,
}
struct TraceRound {
    plan: MatrixPlan,
    leads: Vec<Monomial>,
    live: Vec<usize>,
}

pub(crate) fn learn(
    polynomials: Vec<Polynomial<Zp>>,
) -> Result<(Vec<Polynomial<Zp>>, F4Trace), GroebnerError> {
    let mut trace = F4Trace::default();
    let basis = compute(polynomials, true, Some(&mut trace))?;
    Ok((basis, trace))
}

impl F4Trace {
    /// Replay only independent rows. The caller must validate the reconstructed basis
    /// with a full run independent of this trace before accepting it.
    pub(crate) fn replay(&self, polynomials: Vec<Polynomial<Zp>>) -> Option<Vec<Polynomial<Zp>>> {
        let mut basis = prepare_input(polynomials).ok()?;
        if basis
            .iter()
            .filter_map(|p| p.lm())
            .ne(self.input_leads.iter())
        {
            return None;
        }
        let (nvars, order) = (basis[0].nvars, basis[0].order.clone());
        for round in &self.rounds {
            if round.live.is_empty() {
                continue;
            }
            let matrix = round
                .plan
                .execute_rows(&basis, round.live.iter().copied(), |i| {
                    basis[i].terms.iter().map(|(_, c)| *c).collect()
                })?;
            let rows =
                Zp::echelonize_traced(&matrix.pivots, &matrix.rows, round.plan.columns.len()).0;
            let mut polys: Vec<_> = rows
                .iter()
                .map(|r| decode(r, &round.plan.columns, nvars, &order))
                .collect();
            polys.sort_by(|a, b| crate::groebner::compare_leading(a, b, &order));
            if polys.iter().filter_map(|p| p.lm()).ne(round.leads.iter()) {
                return None;
            }
            basis.extend(polys);
        }
        if basis.len() != self.active.len() {
            return None;
        }
        finish(
            basis
                .into_iter()
                .zip(&self.active)
                .filter_map(|(p, &active)| active.then_some(p))
                .collect(),
            true,
        )
        .ok()
    }

    /// [`F4Trace::replay`] at [`PRIMES`] primes in one pass: the images share every
    /// matrix structure, so one sweep eliminates all their residues side by side. `None`
    /// when any replay would fail or a coefficient vanishes at only some of the primes,
    /// which the caller resolves by replaying the primes one at a time.
    pub(crate) fn replay_lanes(
        &self,
        inputs: Vec<Vec<Polynomial<Zp>>>,
        primes: [u64; PRIMES],
    ) -> Option<Vec<Vec<Polynomial<Zp>>>> {
        let inputs = inputs
            .into_iter()
            .map(|p| prepare_input(p).ok())
            .collect::<Option<Vec<_>>>()?;
        let mut basis = inputs[0].clone();
        let same_support = |a: &Polynomial<Zp>, b: &Polynomial<Zp>| {
            a.terms.len() == b.terms.len() && a.terms.iter().zip(&b.terms).all(|(x, y)| x.0 == y.0)
        };
        if inputs.len() != PRIMES
            || basis
                .iter()
                .filter_map(|p| p.lm())
                .ne(self.input_leads.iter())
            || inputs[1..].iter().any(|other| {
                other.len() != basis.len()
                    || !other.iter().zip(&basis).all(|(a, b)| same_support(a, b))
            })
        {
            return None;
        }
        let mut coefficients: Vec<Vec<Lanes>> = (0..basis.len())
            .map(|i| {
                (0..basis[i].terms.len())
                    .map(|k| {
                        std::array::from_fn(|l| {
                            inputs[l][i].terms[k].1.residue_mod(primes[l]) as u32
                        })
                    })
                    .collect()
            })
            .collect();
        let (nvars, order) = (basis[0].nvars, basis[0].order.clone());
        for round in &self.rounds {
            if round.live.is_empty() {
                continue;
            }
            let matrix = round
                .plan
                .execute_rows(&basis, round.live.iter().copied(), |i| {
                    coefficients[i].clone()
                })?;
            let ncols = round.plan.columns.len();
            let kernel = LaneKernel::new(primes, ncols);
            let known = Table::new(matrix.pivots.view(&matrix.pivots.coefficients), ncols);
            let (rows, _) = echelon(&kernel, &known, &matrix.rows, ncols);
            if kernel.failed.into_inner() {
                return None;
            }
            let mut new: Vec<_> = rows
                .into_iter()
                .map(|r| {
                    let first = SparseRow {
                        columns: r.columns,
                        coefficients: r
                            .coefficients
                            .iter()
                            .map(|v| Zp::from_residue(v[0].into(), primes[0]))
                            .collect(),
                    };
                    (
                        decode(&first, &round.plan.columns, nvars, &order),
                        r.coefficients,
                    )
                })
                .collect();
            new.sort_by(|a, b| crate::groebner::compare_leading(&a.0, &b.0, &order));
            if new
                .iter()
                .filter_map(|(p, _)| p.lm())
                .ne(round.leads.iter())
            {
                return None;
            }
            for (p, c) in new {
                basis.push(p);
                coefficients.push(c);
            }
        }
        if basis.len() != self.active.len() {
            return None;
        }
        (0..PRIMES)
            .map(|l| {
                let polys = basis
                    .iter()
                    .zip(&coefficients)
                    .zip(&self.active)
                    .filter(|&(_, &active)| active)
                    .map(|((p, cs), _)| Polynomial {
                        terms: p
                            .terms
                            .iter()
                            .zip(cs)
                            .filter(|(_, v)| v[l] != 0)
                            .map(|((m, _), v)| {
                                term(Zp::from_residue(v[l].into(), primes[l]), m.clone())
                            })
                            .collect(),
                        nvars,
                        order: order.clone(),
                    })
                    .collect();
                finish(polys, true).ok()
            })
            .collect()
    }
}

/// Number of primes [`F4Trace::replay_lanes`] eliminates together.
pub(crate) const PRIMES: usize = 4;
type Lanes = [u32; PRIMES];

// Dense reduction of lane residues. A leading coefficient that vanishes at only some
// primes splits the images' structure, which is recorded in `failed`.
struct LaneKernel {
    p: [u64; PRIMES],
    ncols: usize,
    deferred: bool,
    failed: std::sync::atomic::AtomicBool,
}
impl LaneKernel {
    fn new(p: [u64; PRIMES], ncols: usize) -> Self {
        let max = p.into_iter().max().unwrap_or(0);
        Self {
            p,
            ncols,
            deferred: max < (1 << 31) && fits(ncols, max),
            failed: false.into(),
        }
    }
}
impl Kernel<Lanes> for LaneKernel {
    type Scratch = Vec<[u64; PRIMES]>;
    fn scratch(&self) -> Self::Scratch {
        vec![[0; PRIMES]; self.ncols]
    }
    fn reduce(
        &self,
        row: &SparseRow<Lanes>,
        pivots: &impl Pivots<Lanes>,
        buf: &mut Self::Scratch,
    ) -> SparseRow<Lanes> {
        for (&c, v) in row.columns.iter().zip(&row.coefficients) {
            buf[c as usize] = v.map(u64::from);
        }
        sweep_lanes(buf, span(&row.columns), pivots, self).unwrap_or_else(|| {
            self.failed.store(true, Ordering::Relaxed);
            SparseRow {
                columns: Vec::new(),
                coefficients: Vec::new(),
            }
        })
    }
    fn normalize(&self, row: &mut SparseRow<Lanes>) {
        let lead = row.coefficients[0];
        if lead == [1; PRIMES] {
            return;
        }
        let inv: [u64; PRIMES] = std::array::from_fn(|l| inv_mod(lead[l].into(), self.p[l]));
        for v in &mut row.coefficients {
            *v = std::array::from_fn(|l| mul_mod(v[l].into(), inv[l], self.p[l]) as u32);
        }
    }
}

// `sweep` over lanes of residues, or `None` if the leading cell vanishes at only some
// primes. Other cells keep the union of the supports, with zeros where a coefficient
// vanishes. `buf` is cleared either way.
fn sweep_lanes(
    buf: &mut [[u64; PRIMES]],
    (lo, mut hi): (usize, usize),
    pivots: &impl Pivots<Lanes>,
    kernel: &LaneKernel,
) -> Option<SparseRow<Lanes>> {
    let p = kernel.p;
    let square = p.map(|p| if p < (1 << 31) { p * p } else { 0 });
    let (mut first, mut left) = (usize::MAX, 0);
    let mut j = lo;
    while j <= hi && j < buf.len() {
        if buf[j] == [0; PRIMES] {
            j += 1;
            continue;
        }
        let Some(piv) = pivots.get(j as u32) else {
            let cell = &mut buf[j];
            for l in 0..PRIMES {
                cell[l] %= p[l];
            }
            if *cell != [0; PRIMES] {
                first = first.min(j);
                left += 1;
            }
            j += 1;
            continue;
        };
        let v = std::mem::take(&mut buf[j]);
        // Multipliers as `u32`, so each product is one widening multiply per lane.
        let c: [u32; PRIMES] = std::array::from_fn(|l| match v[l] % p[l] {
            0 => 0,
            r => (p[l] - r) as u32,
        });
        if c != [0; PRIMES] {
            hi = hi.max(span(piv.columns).1);
            let tail = piv.columns[1..].iter().zip(&piv.coefficients[1..]);
            if kernel.deferred {
                accumulate_lanes(buf, tail, c);
            } else {
                for (&pc, pv) in tail {
                    let cell = &mut buf[pc as usize];
                    for l in 0..PRIMES {
                        let acc = cell[l] + u64::from(c[l]) * u64::from(pv[l]);
                        cell[l] = if acc >= square[l] {
                            acc - square[l]
                        } else {
                            acc
                        };
                    }
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
        if v != [0; PRIMES] {
            out.columns.push(c as u32);
            out.coefficients.push(v.map(|x| x as u32));
        }
    }
    let split = out.coefficients.first().is_some_and(|v| v.contains(&0));
    (!split).then_some(out)
}

// `buf[column] += c * coefficient` in every lane, without reduction.
#[cfg(not(target_arch = "aarch64"))]
fn accumulate_lanes<'a>(
    buf: &mut [[u64; PRIMES]],
    tail: impl Iterator<Item = (&'a u32, &'a Lanes)>,
    c: Lanes,
) {
    for (&pc, pv) in tail {
        let cell = &mut buf[pc as usize];
        for l in 0..PRIMES {
            cell[l] += u64::from(c[l]) * u64::from(pv[l]);
        }
    }
}

// Two widening multiply-accumulates cover the four lanes, where scalar code needs four.
#[cfg(target_arch = "aarch64")]
#[allow(unsafe_code)]
fn accumulate_lanes<'a>(
    buf: &mut [[u64; PRIMES]],
    tail: impl Iterator<Item = (&'a u32, &'a Lanes)>,
    c: Lanes,
) {
    use std::arch::aarch64::{
        vget_low_u32, vld1q_u32, vld1q_u64, vmlal_high_u32, vmlal_u32, vst1q_u64,
    };
    const { assert!(PRIMES == 4) };
    // SAFETY: NEON is part of the aarch64 baseline, and each pointer addresses a whole
    // array of four `u32` or four `u64` lanes.
    unsafe {
        let m = vld1q_u32(c.as_ptr());
        for (&pc, pv) in tail {
            let cell = buf[pc as usize].as_mut_ptr();
            let v = vld1q_u32(pv.as_ptr());
            let lo = vmlal_u32(vld1q_u64(cell), vget_low_u32(v), vget_low_u32(m));
            let hi = vmlal_high_u32(vld1q_u64(cell.add(2)), v, m);
            vst1q_u64(cell, lo);
            vst1q_u64(cell.add(2), hi);
        }
    }
}

#[derive(Debug, Clone)]
struct Pair {
    i: usize,
    j: usize,
    lcm: Monomial,
    degree: u32,
    mask: u64,
}

// A basis element: its coefficients in term order, the packed monomials of its terms when
// all fit, and otherwise, or when a trace records the support, the monomials themselves.
struct Element<F> {
    lm: Monomial,
    coefficients: Vec<F>,
    keys: Option<Box<[u128]>>,
    monomials: Vec<Monomial>,
}

impl<F: Field> Element<F> {
    fn new(p: &Polynomial<F>, traced: bool) -> Self {
        let keys: Option<Box<[u128]>> = p.terms.iter().map(|(m, _)| pack(m.exps())).collect();
        Self {
            lm: p.terms[0].0.clone(),
            coefficients: p.terms.iter().map(|(_, c)| c.clone()).collect(),
            monomials: if keys.is_none() || traced {
                p.support()
            } else {
                Vec::new()
            },
            keys,
        }
    }

    fn decode(
        row: SparseRow<F>,
        columns: &[Monomial],
        keys: &[Option<u128>],
        traced: bool,
    ) -> Self {
        let packed: Option<Box<[u128]>> = row.columns.iter().map(|&c| keys[c as usize]).collect();
        Self {
            lm: columns[row.columns[0] as usize].clone(),
            monomials: if packed.is_none() || traced {
                row.columns
                    .iter()
                    .map(|&c| columns[c as usize].clone())
                    .collect()
            } else {
                Vec::new()
            },
            keys: packed,
            coefficients: row.coefficients,
        }
    }

    fn len(&self) -> usize {
        self.coefficients.len()
    }

    fn term(&self, k: usize, nvars: usize) -> Cow<'_, Monomial> {
        match (self.monomials.get(k), &self.keys) {
            (Some(m), _) => Cow::Borrowed(m),
            (None, Some(keys)) => Cow::Owned(unpack(keys[k], nvars)),
            (None, None) => unreachable!("an element keeps its monomials or their keys"),
        }
    }

    fn polynomial(
        self,
        nvars: usize,
        order: &MonomialOrder,
        cache: &mut HashMap<u128, Monomial>,
    ) -> Polynomial<F> {
        let terms = match (self.monomials.is_empty(), &self.keys) {
            (true, Some(keys)) => keys
                .iter()
                .zip(self.coefficients)
                .map(|(&k, c)| {
                    let m = cache.entry(k).or_insert_with(|| unpack(k, nvars));
                    term(c, m.clone())
                })
                .collect(),
            _ => self
                .monomials
                .into_iter()
                .zip(self.coefficients)
                .map(|(m, c)| term(c, m))
                .collect(),
        };
        Polynomial {
            terms,
            nvars,
            order: order.clone(),
        }
    }
}

struct State<F> {
    basis: Vec<Element<F>>,
    active: Vec<bool>,
    // Active indices by (term count, index), so the first divisor found is the shortest.
    by_size: Vec<usize>,
    pairs: Vec<Pair>,
    masks: Vec<u64>,
    // Packed leading monomials, when they fit.
    keys: Vec<Option<u128>>,
    nvars: usize,
}

struct Matrix<F> {
    pivots: Reducers<F>,
    rows: Vec<SparseRow<F>>,
}

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
const LANES: u128 = 0x8080_8080_8080_8080_8080_8080_8080_8080;

fn pack(exps: &[u32]) -> Option<u128> {
    if exps.len() > 16 {
        return None;
    }
    exps.iter()
        .rev()
        .try_fold(0u128, |acc, &e| (e < 128).then(|| acc << 8 | u128::from(e)))
}

fn unpack(k: u128, nvars: usize) -> Monomial {
    let exps: Vec<u32> = (0..nvars).map(|v| (k >> (8 * v)) as u32 & 0xff).collect();
    Monomial::new(exps.as_slice())
}

// Lanewise maximum of packed monomials: the lcm.
fn lane_max(a: u128, b: u128) -> u128 {
    let ge = ((a | LANES) - b) & LANES;
    let mask = (ge >> 7) * 0xff;
    (a & mask) | (b & !mask)
}

fn lane_divides(a: u128, b: u128) -> bool {
    ((b | LANES) - a) & LANES == LANES
}

fn lane_sum(a: u128) -> u32 {
    const BYTES: u128 = 0x00ff_00ff_00ff_00ff_00ff_00ff_00ff_00ff;
    const ONES: u128 = 0x0001_0001_0001_0001_0001_0001_0001_0001;
    let pairs = (a & BYTES) + ((a >> 8) & BYTES);
    (pairs.wrapping_mul(ONES) >> 112) as u32
}

// Indices of candidates, as lcm, degree, coprime flag and mask, that no other candidate
// eliminates. A candidate goes when an lcm of lower degree divides it, or an equal lcm
// comes later or is coprime. Coprime candidates go too.
fn survivors<T: PartialEq + Clone + Send + Sync>(
    cands: &[(T, u32, bool, u64)],
    divides: impl Fn(&T, &T) -> bool + Sync,
) -> Vec<usize> {
    // Counting sort by degree, stable so equal lcms stay in candidate order.
    let top = cands.iter().map(|c| c.1 as usize).max().unwrap_or(0);
    let mut starts = vec![0; top + 2];
    for c in cands {
        starts[c.1 as usize + 1] += 1;
    }
    for d in 0..=top {
        starts[d + 1] += starts[d];
    }
    let mut next = starts.clone();
    let mut order = vec![0; cands.len()];
    for (p, c) in cands.iter().enumerate() {
        order[next[c.1 as usize]] = p;
        next[c.1 as usize] += 1;
    }
    let sorted: Vec<(T, u64, usize, bool)> = order
        .into_iter()
        .map(|p| (cands[p].0.clone(), cands[p].3, p, cands[p].2))
        .collect();
    let free = map_big(&sorted, |(lp, mp, p, coprime)| {
        let d = cands[*p].1 as usize;
        !coprime
            && !sorted[..starts[d]]
                .iter()
                .any(|(lq, mq, _, _)| mq & !mp == 0 && divides(lq, lp))
            && !sorted[starts[d]..starts[d + 1]]
                .iter()
                .any(|(lq, _, q, cq)| q != p && lq == lp && (q > p || *cq))
    });
    let mut kept: Vec<usize> = sorted
        .iter()
        .zip(free)
        .filter_map(|(c, free)| free.then_some(c.2))
        .collect();
    kept.sort_unstable();
    kept
}

// Ids of packed monomials, probed linearly with each key beside its id, so a lookup
// usually reads one cache line.
#[derive(Default)]
struct PackedIds {
    slots: Vec<([u64; 2], u32)>,
    len: usize,
}
impl PackedIds {
    const EMPTY: u32 = u32::MAX;

    fn start(&self, k: u128) -> usize {
        let h = ((k >> 64) as u64 ^ (k as u64).rotate_left(29)).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        (h >> 32) as usize & (self.slots.len() - 1)
    }

    fn get(&self, k: u128) -> Option<u32> {
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
struct MonomialTable {
    packed: PackedIds,
    ids: HashMap<Key, u32>,
    monomials: Vec<Monomial>,
    scratch: Vec<u32>,
}
impl MonomialTable {
    fn intern(&mut self, m: &Monomial) -> (u32, bool) {
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
    fn product_key(a: &Monomial, b: &Monomial) -> Option<u128> {
        let k = pack(a.exps())? + pack(b.exps())?;
        (k & LANES == 0).then_some(k)
    }

    fn intern_product(&mut self, a: &Monomial, b: &Monomial) -> (u32, bool) {
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
    fn columns(&self, order: &MonomialOrder) -> Vec<usize> {
        let mut ids: Vec<usize> = (0..self.monomials.len()).collect();
        let cmp = |&a: &usize, &b: &usize| order.compare(&self.monomials[b], &self.monomials[a]);
        #[cfg(feature = "parallel")]
        ids.par_sort_unstable_by(cmp);
        #[cfg(not(feature = "parallel"))]
        ids.sort_unstable_by(cmp);
        ids
    }

    #[allow(clippy::expect_used)] // More than 2^32 matrix columns cannot fit in practical memory.
    fn insert(&mut self, m: Monomial) -> u32 {
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
struct Fresh {
    next: AtomicU32,
    shards: Vec<Mutex<HashMap<u128, u32>>>,
}
impl Fresh {
    const SHARDS: usize = 64;

    fn new(next: usize) -> Self {
        Self {
            next: AtomicU32::new(next as u32),
            shards: (0..Self::SHARDS).map(|_| Mutex::default()).collect(),
        }
    }

    fn intern(&self, k: u128) -> u32 {
        let h = ((k >> 64) as u64 ^ k as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let shard = &self.shards[(h >> 58) as usize];
        let mut map = shard.lock().unwrap_or_else(PoisonError::into_inner);
        *map.entry(k)
            .or_insert_with(|| self.next.fetch_add(1, Ordering::Relaxed))
    }

    // Append the new monomials to the table in id order, returning their ids.
    fn drain(self, table: &mut MonomialTable, nvars: usize) -> Vec<u32> {
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

struct PlannedRow {
    basis: usize,
    // Support at planning time, kept only when recording a trace.
    terms: Vec<Monomial>,
    // Taken by the matrix unless the plan is traced.
    columns: Vec<u32>,
}
struct MatrixPlan {
    columns: Vec<Monomial>,
    pivots: Vec<PlannedRow>,
    rows: Vec<PlannedRow>,
}
impl MatrixPlan {
    // Some rows of a traced plan, over a replayed basis whose coefficients, in term
    // order, are `coefficients(i)`.
    fn execute_rows<F: Field, C>(
        &self,
        basis: &[Polynomial<F>],
        rows: impl Iterator<Item = usize>,
        coefficients: impl Fn(usize) -> Vec<C>,
    ) -> Option<Matrix<C>> {
        // Rows of one element share its planned support: an element that kept all of
        // it, as most do, takes the planned columns without matching monomials again.
        let mut kept = vec![None; basis.len()];
        // Columns of the basis polynomial, which may have lost terms since planning.
        let mut columns = |row: &PlannedRow| {
            let p = basis.get(row.basis)?;
            let same =
                *kept[row.basis].get_or_insert_with(|| p.terms.iter().map(|t| &t.0).eq(&row.terms));
            if same {
                return Some(row.columns.clone());
            }
            let mut columns = Vec::with_capacity(p.terms.len());
            let mut k = 0;
            for (m, _) in &p.terms {
                while row.terms.get(k).is_some_and(|n| n != m) {
                    k += 1;
                }
                row.terms.get(k)?;
                columns.push(row.columns[k]);
            }
            // Every known pivot must retain its learned leading column.
            (columns.first() == row.columns.first()).then_some(columns)
        };
        let pivots = self
            .pivots
            .iter()
            .map(|row| Some((row.basis, columns(row)?)))
            .collect::<Option<Vec<_>>>()?;
        let rows = rows
            .map(|i| {
                let row = &self.rows[i];
                Some(SparseRow {
                    columns: columns(row)?,
                    coefficients: coefficients(row.basis),
                })
            })
            .collect::<Option<_>>()?;
        Some(Matrix {
            pivots: reducers(basis.len(), coefficients, pivots),
            rows,
        })
    }
}

impl<F: F4Field> State<F> {
    fn lm(&self, i: usize) -> &Monomial {
        &self.basis[i].lm
    }

    // The matrix of a plan made from this basis. Unless the plan is traced, each row's
    // columns are freed as they are copied, so only one row is held twice, and the copies
    // lie together, which elimination reads faster than the plan's scattered vectors.
    fn matrix(&self, plan: &mut MatrixPlan, traced: bool) -> Matrix<F> {
        let columns = |row: &mut PlannedRow| {
            if traced {
                row.columns.clone()
            } else {
                std::mem::take(&mut row.columns).as_slice().to_vec()
            }
        };
        let pivots: Vec<_> = plan
            .pivots
            .iter_mut()
            .map(|row| (row.basis, columns(row)))
            .collect();
        let rows = plan
            .rows
            .iter_mut()
            .map(|row| SparseRow {
                columns: columns(row),
                coefficients: self.basis[row.basis].coefficients.clone(),
            })
            .collect();
        Matrix {
            pivots: reducers(
                self.basis.len(),
                |i| self.basis[i].coefficients.clone(),
                pivots,
            ),
            rows,
        }
    }

    fn pair(&self, i: usize, j: usize) -> Pair {
        let lcm = self.lm(i).lcm(self.lm(j));
        Pair {
            i,
            j,
            degree: lcm.degree(),
            mask: divisibility_mask(&lcm),
            lcm,
        }
    }

    // Gebauer-Moller update (Becker and Weispfenning, section 5.5).
    fn update(&mut self, mut h: Element<F>) {
        if !h.coefficients[0].is_one()
            && let Some(inv) = h.coefficients[0].inverse()
        {
            for c in &mut h.coefficients {
                *c = c.clone() * inv.clone();
            }
        }
        let lm_h = h.lm.clone();
        let t = self.basis.len();
        let mask_h = divisibility_mask(&lm_h);
        self.masks.push(mask_h);
        self.keys.push(pack(lm_h.exps()));
        self.basis.push(h);
        self.active.push(true);
        let kept = self.new_pairs(t);

        let basis = &self.basis;
        let lm = |i: usize| &basis[i].lm;
        let lcm_is = |a: &Monomial, l: &Monomial| {
            a.exps()
                .iter()
                .zip(lm_h.exps())
                .zip(l.exps())
                .all(|((x, y), z)| x.max(y) == z)
        };
        self.pairs.retain(|p| {
            mask_h & !p.mask != 0
                || !lm_h.divides(&p.lcm)
                || lcm_is(lm(p.i), &p.lcm)
                || lcm_is(lm(p.j), &p.lcm)
        });
        self.pairs.extend(kept);
        for g in 0..t {
            if self.active[g] && mask_h & !self.masks[g] == 0 && lm_h.divides(lm(g)) {
                self.active[g] = false;
            }
        }
        let active = &self.active;
        self.by_size.retain(|&g| active[g]);
        let size = |g: usize| (basis[g].len(), g);
        let at = self.by_size.partition_point(|&g| size(g) < size(t));
        self.by_size.insert(at, t);
    }

    // Pairs of `t` with the active basis that survive the chain criterion. A pair goes
    // when another's lcm properly divides its own, or equals it and either comes later
    // or has coprime leading monomials. Coprime pairs go too, after serving as divisors.
    fn new_pairs(&self, t: usize) -> Vec<Pair> {
        let lm_h = self.lm(t);
        let others: Vec<usize> = (0..t).filter(|&g| self.active[g]).collect();
        // The mask of an lcm is the union of the masks of its arguments.
        let mask = |g: usize| self.masks[g] | self.masks[t];
        let kept = match self.keys[t] {
            Some(h) if others.iter().all(|&g| self.keys[g].is_some()) => {
                let cands = map_big(&others, |&g| {
                    let k = self.keys[g].unwrap_or_default();
                    let lcm = lane_max(k, h);
                    (lcm, lane_sum(lcm), lcm == k + h, mask(g))
                });
                survivors(&cands, |&a, &b| lane_divides(a, b))
            }
            _ => {
                let cands = map_big(&others, |&g| {
                    let lcm = self.lm(g).lcm(lm_h);
                    let degree = lcm.degree();
                    (lcm, degree, self.lm(g).is_coprime(lm_h), mask(g))
                });
                survivors(&cands, Monomial::divides)
            }
        };
        kept.into_iter().map(|k| self.pair(others[k], t)).collect()
    }

    fn select(&mut self) -> Vec<Pair> {
        let min_degree = self
            .pairs
            .iter()
            .map(|p| p.degree)
            .min()
            .unwrap_or_default();
        let (selected, rest) = std::mem::take(&mut self.pairs)
            .into_iter()
            .partition(|p| p.degree == min_degree);
        self.pairs = rest;
        selected
    }

    fn reducer(&self, m: &Monomial) -> Option<usize> {
        let mask = divisibility_mask(m);
        let key = pack(m.exps());
        self.by_size.iter().copied().find(|&i| {
            self.masks[i] & !mask == 0
                && match (self.keys[i], key) {
                    (Some(a), Some(b)) => lane_divides(a, b),
                    _ => self.lm(i).divides(m),
                }
        })
    }

    fn symbolic_preprocessing(
        &self,
        pairs: &[Pair],
        order: &MonomialOrder,
        traced: bool,
    ) -> MatrixPlan {
        let mut table = MonomialTable::default();
        let mut multiples: HashSet<(usize, Monomial)> = HashSet::default();
        let mut by_lcm: HashMap<u32, Vec<(usize, Monomial)>> = HashMap::default();
        for p in pairs {
            let (id, _) = table.intern(&p.lcm);
            for i in [p.i, p.j] {
                if let Some(mult) = p.lcm.quo(self.lm(i))
                    && multiples.insert((i, mult.clone()))
                {
                    by_lcm.entry(id).or_default().push((i, mult));
                }
            }
        }
        let mut pivot_leads: HashSet<u32> = by_lcm.keys().copied().collect();
        let mut groups: Vec<_> = by_lcm.into_iter().collect();
        groups.sort_by(|(a, _), (b, _)| {
            order.compare(&table.monomials[*b as usize], &table.monomials[*a as usize])
        });
        // Rows to build, flagged when they are pivots: the shortest multiple per lcm.
        let mut level: Vec<(usize, Monomial, bool)> = Vec::new();
        for (_, mut group) in groups {
            group.sort_by_key(|(i, _)| (self.basis[*i].len(), *i));
            level.extend(
                group
                    .into_iter()
                    .enumerate()
                    .map(|(k, (i, m))| (i, m, k == 0)),
            );
        }
        let nvars = self.nvars;
        let mut pivots = Vec::new();
        let mut rows = Vec::new();
        // Each level interns its products in parallel and adds a reducer for each new
        // monomial as the next level.
        while !level.is_empty() {
            let fresh = Fresh::new(table.monomials.len());
            let found = map_rows(
                &level,
                // Each worker caches the fresh ids it has seen, sparing the shard locks.
                HashMap::<u128, u32>::default,
                |(i, mult, _), seen| {
                    let mut column = |k: Option<u128>| match k.filter(|k| k & LANES == 0) {
                        None => u32::MAX,
                        Some(k) => match table.packed.get(k) {
                            Some(id) => id,
                            None => *seen.entry(k).or_insert_with(|| fresh.intern(k)),
                        },
                    };
                    let g = &self.basis[*i];
                    let columns: Vec<u32> = match (&g.keys, pack(mult.exps())) {
                        (Some(keys), Some(m)) => keys.iter().map(|k| column(Some(k + m))).collect(),
                        _ => (0..g.len())
                            .map(|k| column(MonomialTable::product_key(&g.term(k, nvars), mult)))
                            .collect(),
                    };
                    let unpacked = columns.contains(&u32::MAX);
                    (columns, unpacked)
                },
            );
            let mut fresh = fresh.drain(&mut table, nvars);
            for ((i, mult, pivot), (mut columns, unpacked)) in level.into_iter().zip(found) {
                // Products too large to pack are interned here.
                let g = &self.basis[i];
                for (k, c) in columns.iter_mut().enumerate().filter(|_| unpacked) {
                    if *c == u32::MAX {
                        let (id, new) = table.intern_product(&g.term(k, nvars), &mult);
                        if new {
                            fresh.push(id);
                        }
                        *c = id;
                    }
                }
                let row = PlannedRow {
                    basis: i,
                    terms: if traced {
                        g.monomials.clone()
                    } else {
                        Vec::new()
                    },
                    columns,
                };
                if pivot {
                    pivots.push(row)
                } else {
                    rows.push(row)
                }
            }
            fresh.retain(|&id| pivot_leads.insert(id));
            let found = map_rows(
                &fresh,
                || (),
                |&id, ()| {
                    let m = &table.monomials[id as usize];
                    let i = self.reducer(m)?;
                    Some((i, m.quo(self.lm(i))?, true))
                },
            );
            level = found.into_iter().flatten().collect();
        }
        let ids = table.columns(order);
        let index = column_index(&ids);
        let columns = ids
            .into_iter()
            .map(|id| table.monomials[id].clone())
            .collect();
        let reindex = |row: &mut PlannedRow| {
            for c in &mut row.columns {
                *c = index[*c as usize];
            }
        };
        for_rows_mut(&mut pivots, reindex);
        for_rows_mut(&mut rows, reindex);
        pivots.sort_by_key(|r| r.columns[0]);
        MatrixPlan {
            columns,
            pivots,
            rows,
        }
    }
}

/// Minimize and, when `canonicalize` is set, interreduce with one Macaulay matrix: the tails
/// of the minimal basis are the rows, and symbolic preprocessing supplies the reducers.
fn finish<F: F4Field>(
    basis: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    let basis = finish_basis(basis, false)?;
    if !canonicalize {
        return Ok(basis);
    }
    let first = basis.first().ok_or(GroebnerError::EmptyInput)?;
    let (nvars, order) = (first.nvars, first.order.clone());
    let basis: Vec<_> = basis
        .into_iter()
        .map(|g| {
            if g.leading_coefficient().is_some_and(num_traits::One::is_one) {
                g
            } else {
                g.make_monic()
            }
        })
        .collect();
    let masks: Vec<u64> = basis
        .iter()
        .map(|g| divisibility_mask(&g.terms[0].0))
        .collect();
    let reducer = |m: &Monomial| {
        let mask = divisibility_mask(m);
        (0..basis.len())
            .filter(|&i| masks[i] & !mask == 0 && basis[i].terms[0].0.divides(m))
            .min_by_key(|&i| basis[i].terms.len())
    };
    let mut table = MonomialTable::default();
    let mut queue = Vec::new();
    let fresh = |(id, fresh): (u32, bool), queue: &mut Vec<u32>| {
        if fresh {
            queue.push(id);
        }
        id
    };
    let tails: Vec<Vec<u32>> = basis
        .iter()
        .map(|g| {
            g.terms[1..]
                .iter()
                .map(|(m, _)| fresh(table.intern(m), &mut queue))
                .collect()
        })
        .collect();
    let mut pivots: Vec<(usize, Vec<u32>)> = Vec::new();
    while let Some(id) = queue.pop() {
        let m = table.monomials[id as usize].clone();
        let Some(i) = reducer(&m) else { continue };
        let Some(mult) = m.quo(&basis[i].terms[0].0) else {
            continue;
        };
        let mut columns = vec![id];
        for (t, _) in &basis[i].terms[1..] {
            columns.push(fresh(table.intern_product(t, &mult), &mut queue));
        }
        pivots.push((i, columns));
    }
    let ids = table.columns(&order);
    let index = column_index(&ids);
    let columns: Vec<Monomial> = ids.iter().map(|&id| table.monomials[id].clone()).collect();
    let reindex = |mut cols: Vec<u32>| {
        cols.iter_mut().for_each(|c| *c = index[*c as usize]);
        cols
    };
    let pivot_rows = reducers(
        basis.len(),
        |i| basis[i].terms.iter().map(|(_, c)| c.clone()).collect(),
        pivots.into_iter().map(|(i, cols)| (i, reindex(cols))),
    );
    let tail_rows: Vec<_> = basis
        .iter()
        .zip(tails)
        .map(|(g, cols)| SparseRow {
            columns: reindex(cols),
            coefficients: g.terms[1..].iter().map(|(_, c)| c.clone()).collect(),
        })
        .collect();
    // The rows hold everything else, so only the leading terms are kept.
    let leads: Vec<_> = basis
        .into_iter()
        .map(|mut g| g.terms.swap_remove(0))
        .collect();
    let reduced = F::reduce_rows(&pivot_rows, &tail_rows, columns.len());
    drop((pivot_rows, tail_rows));
    let mut out: Vec<Polynomial<F>> = leads
        .into_iter()
        .zip(reduced)
        .map(|(lead, tail)| {
            let mut p = decode(&tail, &columns, nvars, &order);
            p.terms.insert(0, lead);
            p
        })
        .collect();
    out.sort_by(|a, b| crate::groebner::compare_leading(b, a, &order));
    Ok(out)
}

// Column of each interned monomial, given the ids in column order.
fn column_index(ids: &[usize]) -> Vec<u32> {
    let mut index = vec![0; ids.len()];
    for (col, &id) in ids.iter().enumerate() {
        index[id] = col as u32;
    }
    index
}

fn decode<F: Field>(
    row: &SparseRow<F>,
    columns: &[Monomial],
    nvars: usize,
    order: &MonomialOrder,
) -> Polynomial<F> {
    let terms = row
        .columns
        .iter()
        .zip(&row.coefficients)
        .map(|(&c, v)| term(v.clone(), columns[c as usize].clone()))
        .collect();
    Polynomial {
        terms,
        nvars,
        order: order.clone(),
    }
}

// A borrowed row, so known pivots can share coefficients.
struct Row<'a, F> {
    columns: &'a [u32],
    coefficients: &'a [F],
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
trait Pivots<F>: Sync {
    fn get(&self, col: u32) -> Option<Row<'_, F>>;
}

struct Table<'a, F> {
    rows: Vec<Row<'a, F>>,
    index: Vec<u32>,
}
impl<'a, F> Table<'a, F> {
    fn new(rows: impl IntoIterator<Item = Row<'a, F>>, ncols: usize) -> Self {
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
trait Kernel<F>: Sync {
    type Scratch;
    fn scratch(&self) -> Self::Scratch;
    fn reduce(
        &self,
        row: &SparseRow<F>,
        pivots: &impl Pivots<F>,
        scratch: &mut Self::Scratch,
    ) -> SparseRow<F>;
    fn normalize(&self, row: &mut SparseRow<F>);
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

// Map rows in parallel blocks, each worker reusing one scratch value.
fn map_rows<T: Sync, R: Send, S>(
    rows: &[T],
    scratch: impl Fn() -> S + Sync,
    f: impl Fn(&T, &mut S) -> R + Sync,
) -> Vec<R> {
    let block = |block: &[T]| {
        let mut s = scratch();
        block.iter().map(|r| f(r, &mut s)).collect::<Vec<_>>()
    };
    #[cfg(feature = "parallel")]
    if rows.len() >= 64 {
        return rows.par_chunks(32).flat_map_iter(block).collect();
    }
    block(rows)
}

fn for_rows_mut<T: Send>(rows: &mut [T], f: impl Fn(&mut T) + Sync) {
    #[cfg(feature = "parallel")]
    if rows.len() >= 64 {
        return rows
            .par_chunks_mut(32)
            .for_each(|block| block.iter_mut().for_each(&f));
    }
    rows.iter_mut().for_each(f);
}

// `map_rows` for cheap per-row work, which only pays to split when there is a lot.
fn map_big<T: Sync, R: Send>(rows: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    if rows.len() < 1024 {
        return rows.iter().map(f).collect();
    }
    map_rows(rows, || (), |r, ()| f(r))
}

// Reduce `rows` by the monic `pivots` and by each other, then back substitute. Workers
// claim free leading columns concurrently, and a row that loses a claim keeps reducing
// by the winner. The reduced echelon form is unique, so the result does not depend on
// scheduling, though which rows are reported live (independent) may.
fn echelon<F: Clone + Send + Sync, K: Kernel<F>>(
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
    map_rows(
        &order,
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
    let reduced = map_rows(
        &new,
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

fn echelonize_generic<F: Field + Send + Sync>(
    pivots: &Reducers<F>,
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    let known = Table::new(pivots.view(&pivots.coefficients), ncols);
    echelon(&Generic(ncols), &known, rows, ncols).0
}

// Span of a row's columns, as buffer indices.
fn span(columns: &[u32]) -> (usize, usize) {
    match (columns.first(), columns.last()) {
        (Some(&lo), Some(&hi)) => (lo as usize, hi as usize),
        _ => (usize::MAX, 0),
    }
}

fn reduce_generic<F: Field>(
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

fn echelonize_modular<F: ModularField + Send + Sync>(
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
    let zeros = 40u32.div_ceil(p.ilog2()).max(2);
    let block = |(b, block): (usize, &&[usize]), buf: &mut Vec<u64>| {
        // Unreduced sums only when the sweep can also defer every reduction.
        let lazy = p < (1 << 31) && fits(ncols + block.len(), p);
        let mut rng = Rng::new(b as u64);
        let (mut found, mut run) = (0, 0);
        while found < block.len() && run < zeros {
            let (mut lo, mut hi) = (usize::MAX, 0);
            for &i in *block {
                let m = rng.nonzero(p);
                let (cols, coefs) = (&rows[i].columns, &rows[i].coefficients);
                for (&c, &v) in cols.iter().zip(coefs) {
                    let cell = &mut buf[c as usize];
                    *cell = if lazy {
                        *cell + m * v.into()
                    } else {
                        add_mod(*cell, mul_mod(m, v.into(), p), p)
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
    let scratch = || Kernel::<C>::scratch(&kernel);
    #[cfg(feature = "parallel")]
    blocks
        .par_iter()
        .enumerate()
        .for_each_init(scratch, |buf, item| block(item, buf));
    #[cfg(not(feature = "parallel"))]
    {
        let mut buf = scratch();
        blocks
            .iter()
            .enumerate()
            .for_each(|item| block(item, &mut buf));
    }
    back_substitute(&kernel, claimed, ncols).0
}

// Whether `terms` products of residues plus one residue fit in 64 bits.
fn fits(terms: usize, p: u64) -> bool {
    (terms as u128) * u128::from(p - 1).pow(2) + u128::from(p - 1) <= u128::from(u64::MAX)
}

fn echelonize_modular_traced<F: ModularField + Send + Sync>(
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

fn reduce_rows_modular<F: ModularField + Send + Sync>(
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
    let reduced = map_rows(
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
            } else {
                for (&pc, &pv) in tail {
                    buf[pc as usize] = add_mod(buf[pc as usize], mul_mod(c, pv.into(), p), p);
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
mod trace_tests {
    use super::*;
    use crate::PolynomialRing;
    fn input(p: u64, system: &str) -> Vec<Polynomial<Zp>> {
        PolynomialRing::with_modulus(["x", "y"], MonomialOrder::GRevLex, p)
            .expect("ring")
            .parse_many(system)
            .expect("system")
    }
    #[test]
    fn trace_replays_at_another_prime() {
        let system = "x^3 - y; x*y - 1; y^3 - x";
        let (_, trace) = learn(input(32003, system)).expect("learn");
        assert!(!trace.rounds.is_empty());
        let polys = input(32009, system);
        assert_eq!(
            trace.replay(polys.clone()).expect("replay"),
            groebner_basis_f4_direct(polys, true).expect("full")
        );
    }
    #[test]
    fn omitted_zero_rows_require_independent_validation() {
        let system = "x + y; x - y";
        let (_, trace) = learn(input(2, system)).expect("learn");
        assert!(
            trace
                .rounds
                .iter()
                .any(|r| r.live.len() < r.plan.rows.len())
        );
        let replay = trace
            .replay(input(3, system))
            .expect("structurally valid trace");
        let full = groebner_basis_f4_direct(input(3, system), true).expect("full");
        assert_ne!(
            replay, full,
            "independent validation must reject this trace"
        );
    }
    #[test]
    fn replay_rejects_a_missing_pivot() {
        let system = "x + y; x - y";
        let (_, trace) = learn(input(3, system)).expect("learn");
        assert!(trace.replay(input(2, system)).is_none());
        let primes = [3, 5, 7, 2];
        let inputs = primes.iter().map(|&p| input(p, system)).collect();
        assert!(trace.replay_lanes(inputs, primes).is_none());
    }
    #[test]
    fn lane_replay_matches_scalar_replay() {
        let system = "x^3 - 2*y + 5; x*y^2 - 3*x + 1; y^3 - x^2 + 7*y";
        let (_, trace) = learn(input(32003, system)).expect("learn");
        let primes = [8_388_593, 8_388_587, 8_388_581, 2_147_483_647];
        let inputs = primes.iter().map(|&p| input(p, system)).collect();
        let lanes = trace.replay_lanes(inputs, primes).expect("lanes");
        for (basis, p) in lanes.into_iter().zip(primes) {
            assert_eq!(Some(basis), trace.replay(input(p, system)));
        }
    }
    #[test]
    fn modular_row_blocks_match_generic_elimination_at_word_boundaries() {
        for p in [2, 32003, 2_147_483_647, 18_446_744_073_709_551_557] {
            let ncols = 80;
            let mut seed = 42u64;
            let mut next = || {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
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
