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
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::borrow::Borrow;
use std::collections::VecDeque;
use std::hash::{Hash, Hasher};
use std::sync::OnceLock;

/// A sparse matrix row with columns in ascending index order (descending monomial order).
#[derive(Debug, Clone)]
pub struct SparseRow<F> {
    pub columns: Vec<u32>,
    pub coefficients: Vec<F>,
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
        pivots: &[SparseRow<Self>],
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
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_generic(pivots, rows, ncols)
    }

    /// Fully reduce each row modulo the monic `pivots` (distinct leading columns).
    fn reduce_rows(
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        let table = Table::new(pivots, ncols);
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
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> (Vec<SparseRow<Self>>, Vec<usize>) {
        echelonize_modular_traced(pivots, rows, ncols)
    }

    fn echelonize(
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_modular(pivots, rows, ncols)
    }

    fn reduce_rows(
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        reduce_rows_modular(pivots, rows, ncols)
    }
}

impl F4Field for Zp {
    fn echelonize_traced(
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> (Vec<SparseRow<Self>>, Vec<usize>) {
        echelonize_modular_traced(pivots, rows, ncols)
    }

    fn echelonize(
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_modular(pivots, rows, ncols)
    }

    fn reduce_rows(
        pivots: &[SparseRow<Self>],
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

fn compute<F: F4Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
    mut trace: Option<&mut F4Trace>,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    let input = prepare_input(polynomials)?;
    let order = input[0].order.clone();
    let nvars = input[0].nvars;
    let mut state = State {
        basis: Vec::new(),
        active: Vec::new(),
        by_size: Vec::new(),
        pairs: Vec::new(),
        masks: Vec::new(),
    };
    for poly in input {
        state.update(poly);
    }
    if let Some(t) = trace.as_deref_mut() {
        t.input_leads = state.basis.iter().filter_map(|p| p.lm().cloned()).collect();
    }
    while !state.pairs.is_empty() {
        let selected = state.select();
        let plan = state.symbolic_preprocessing(&selected, &order, trace.is_some());
        let matrix = plan
            .execute(&state.basis)
            .unwrap_or_else(|| unreachable!("fresh matrix plan"));
        let (new_rows, live) = if trace.is_some() {
            F::echelonize_traced(&matrix.pivots, &matrix.rows, plan.columns.len())
        } else {
            (
                F::echelonize(&matrix.pivots, &matrix.rows, plan.columns.len()),
                Vec::new(),
            )
        };
        let mut new_polys: Vec<Polynomial<F>> = new_rows
            .into_iter()
            .map(|row| decode(&row, &plan.columns, nvars, &order))
            .collect();
        new_polys.sort_by(|a, b| crate::groebner::compare_leading(a, b, &order));
        if let Some(t) = trace.as_deref_mut() {
            let leads = new_polys.iter().filter_map(|p| p.lm().cloned()).collect();
            t.rounds.push(TraceRound { plan, leads, live });
        }
        for poly in new_polys {
            state.update(poly);
        }
    }
    if let Some(t) = trace {
        t.active = state.active.clone();
    }
    let basis = state
        .basis
        .into_iter()
        .zip(state.active)
        .filter_map(|(poly, active)| active.then_some(poly))
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
                .execute_rows(&basis, round.live.iter().copied())?;
            let rows = Zp::echelonize(&matrix.pivots, &matrix.rows, round.plan.columns.len());
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
}

#[derive(Debug, Clone)]
struct Pair {
    i: usize,
    j: usize,
    lcm: Monomial,
    degree: u32,
    mask: u64,
}

struct State<F> {
    basis: Vec<Polynomial<F>>,
    active: Vec<bool>,
    // Active indices by (term count, index), so the first divisor found is the shortest.
    by_size: Vec<usize>,
    pairs: Vec<Pair>,
    masks: Vec<u64>,
}

struct Matrix<F> {
    pivots: Vec<SparseRow<F>>,
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

#[derive(Default)]
struct MonomialTable {
    ids: HashMap<Key, u32>,
    monomials: Vec<Monomial>,
    scratch: Vec<u32>,
}
impl MonomialTable {
    fn intern(&mut self, m: &Monomial) -> (u32, bool) {
        match self.ids.get(m.exps()) {
            Some(&id) => (id, false),
            None => (self.insert(m.clone()), true),
        }
    }

    fn intern_product(&mut self, a: &Monomial, b: &Monomial) -> (u32, bool) {
        self.scratch.clear();
        self.scratch
            .extend(a.exps().iter().zip(b.exps()).map(|(x, y)| x + y));
        match self.ids.get(self.scratch.as_slice()) {
            Some(&id) => (id, false),
            None => (self.insert(Monomial::new(self.scratch.as_slice())), true),
        }
    }

    #[allow(clippy::expect_used)] // More than 2^32 matrix columns cannot fit in practical memory.
    fn insert(&mut self, m: Monomial) -> u32 {
        let id = u32::try_from(self.monomials.len()).expect("F4 matrix exceeds u32 columns");
        self.monomials.push(m.clone());
        self.ids.insert(Key(m), id);
        id
    }
}

struct PlannedRow {
    basis: usize,
    // Support at planning time, kept only when recording a trace.
    terms: Vec<Monomial>,
    columns: Vec<u32>,
}
struct MatrixPlan {
    columns: Vec<Monomial>,
    pivots: Vec<PlannedRow>,
    rows: Vec<PlannedRow>,
}
impl MatrixPlan {
    fn execute<F: Field>(&self, basis: &[Polynomial<F>]) -> Option<Matrix<F>> {
        self.execute_rows(basis, 0..self.rows.len())
    }

    fn execute_rows<F: Field>(
        &self,
        basis: &[Polynomial<F>],
        rows: impl Iterator<Item = usize>,
    ) -> Option<Matrix<F>> {
        let encode = |row: &PlannedRow| {
            let p = basis.get(row.basis)?;
            let mut coefficients = Vec::with_capacity(p.terms.len());
            let mut columns = Vec::with_capacity(p.terms.len());
            if row.terms.is_empty() {
                // Planned from this basis without a trace: terms align with columns.
                return Some(SparseRow {
                    columns: row.columns.clone(),
                    coefficients: p.terms.iter().map(|(_, c)| c.clone()).collect(),
                });
            }
            let mut k = 0;
            for (m, c) in &p.terms {
                while row.terms.get(k).is_some_and(|n| n != m) {
                    k += 1;
                }
                row.terms.get(k)?;
                columns.push(row.columns[k]);
                coefficients.push(c.clone());
            }
            // Every known pivot must retain its learned leading column.
            if columns.first() != row.columns.first() {
                return None;
            }
            Some(SparseRow {
                columns,
                coefficients,
            })
        };
        Some(Matrix {
            pivots: self.pivots.iter().map(encode).collect::<Option<_>>()?,
            rows: rows.map(|i| encode(&self.rows[i])).collect::<Option<_>>()?,
        })
    }
}

impl<F: F4Field> State<F> {
    fn lm(&self, i: usize) -> &Monomial {
        self.basis[i]
            .leading_monomial()
            .unwrap_or_else(|| unreachable!("basis polynomials are nonzero"))
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
    fn update(&mut self, h: Polynomial<F>) {
        let h = if h.leading_coefficient().is_some_and(num_traits::One::is_one) {
            h
        } else {
            h.make_monic()
        };
        let Some(lm_h) = h.leading_monomial().cloned() else {
            return;
        };
        let t = self.basis.len();
        let mask_h = divisibility_mask(&lm_h);
        self.masks.push(mask_h);
        self.basis.push(h);
        self.active.push(true);

        let mut candidates: VecDeque<Pair> = (0..t)
            .filter(|&g| self.active[g])
            .map(|g| self.pair(g, t))
            .collect();
        let mut kept: Vec<Pair> = Vec::new();
        while let Some(p) = candidates.pop_front() {
            let coprime = self.lm(p.i).is_coprime(&lm_h);
            let dominated = |q: &Pair| q.mask & !p.mask == 0 && q.lcm.divides(&p.lcm);
            if coprime || !(candidates.iter().any(dominated) || kept.iter().any(dominated)) {
                kept.push(p);
            }
        }
        kept.retain(|p| !self.lm(p.i).is_coprime(&lm_h));

        let basis = &self.basis;
        let lm = |i: usize| {
            basis[i]
                .leading_monomial()
                .unwrap_or_else(|| unreachable!("basis polynomials are nonzero"))
        };
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
        let size = |g: usize| (basis[g].terms.len(), g);
        let at = self.by_size.partition_point(|&g| size(g) < size(t));
        self.by_size.insert(at, t);
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
        self.by_size
            .iter()
            .copied()
            .find(|&i| self.masks[i] & !mask == 0 && self.lm(i).divides(m))
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
        let mut pivots = Vec::new();
        let mut rows = Vec::new();
        let mut queue: Vec<u32> = Vec::new();
        let visit = |i: usize, mult: &Monomial, table: &mut MonomialTable, queue: &mut Vec<u32>| {
            let mut columns = Vec::with_capacity(self.basis[i].terms.len());
            for (m, _) in &self.basis[i].terms {
                let (id, fresh) = table.intern_product(m, mult);
                if fresh {
                    queue.push(id);
                }
                columns.push(id);
            }
            PlannedRow {
                basis: i,
                terms: if traced {
                    self.basis[i].support()
                } else {
                    Vec::new()
                },
                columns,
            }
        };
        let mut pivot_leads: HashSet<u32> = by_lcm.keys().copied().collect();
        let mut groups: Vec<_> = by_lcm.into_iter().collect();
        groups.sort_by(|(a, _), (b, _)| {
            order.compare(&table.monomials[*b as usize], &table.monomials[*a as usize])
        });
        for (_, mut group) in groups {
            group.sort_by_key(|(i, _)| (self.basis[*i].terms.len(), *i));
            for (k, (i, mult)) in group.into_iter().enumerate() {
                let row = visit(i, &mult, &mut table, &mut queue);
                if k == 0 {
                    pivots.push(row);
                } else {
                    rows.push(row);
                }
            }
        }
        while let Some(id) = queue.pop() {
            if pivot_leads.contains(&id) {
                continue;
            }
            let m = &table.monomials[id as usize];
            if let Some(i) = self.reducer(m)
                && let Some(mult) = m.quo(self.lm(i))
            {
                pivots.push(visit(i, &mult, &mut table, &mut queue));
                pivot_leads.insert(id);
            }
        }
        let mut ids: Vec<usize> = (0..table.monomials.len()).collect();
        ids.sort_by(|&a, &b| order.compare(&table.monomials[b], &table.monomials[a]));
        let index = column_index(&ids);
        let columns = ids
            .into_iter()
            .map(|id| table.monomials[id].clone())
            .collect();
        for row in pivots.iter_mut().chain(&mut rows) {
            for c in &mut row.columns {
                *c = index[*c as usize];
            }
        }
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
    let basis: Vec<_> = basis.iter().map(PolynomialExt::make_monic).collect();
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
    let mut tails: Vec<Vec<u32>> = basis
        .iter()
        .map(|g| {
            g.terms[1..]
                .iter()
                .map(|(m, _)| fresh(table.intern(m), &mut queue))
                .collect()
        })
        .collect();
    let mut pivots: Vec<(usize, Monomial, Vec<u32>)> = Vec::new();
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
        pivots.push((i, mult, columns));
    }
    let mut ids: Vec<usize> = (0..table.monomials.len()).collect();
    ids.sort_by(|&a, &b| order.compare(&table.monomials[b], &table.monomials[a]));
    let index = column_index(&ids);
    let columns: Vec<Monomial> = ids.iter().map(|&id| table.monomials[id].clone()).collect();
    let row = |g: &Polynomial<F>, cols: &mut [u32], skip: usize| {
        cols.iter_mut().for_each(|c| *c = index[*c as usize]);
        SparseRow {
            columns: cols.to_vec(),
            coefficients: g.terms[skip..].iter().map(|(_, c)| c.clone()).collect(),
        }
    };
    let pivot_rows: Vec<_> = pivots
        .iter_mut()
        .map(|(i, _, cols)| row(&basis[*i], cols, 0))
        .collect();
    let tail_rows: Vec<_> = basis
        .iter()
        .zip(&mut tails)
        .map(|(g, cols)| row(g, cols, 1))
        .collect();
    let reduced = F::reduce_rows(&pivot_rows, &tail_rows, columns.len());
    let mut out: Vec<Polynomial<F>> = basis
        .into_iter()
        .zip(reduced)
        .map(|(g, tail)| {
            let mut p = decode(&tail, &columns, nvars, &order);
            p.terms.insert(0, g.terms[0].clone());
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

/// Pivot rows indexed by leading column.
trait Pivots<F>: Sync {
    fn get(&self, col: u32) -> Option<&SparseRow<F>>;
}

struct Table<'a, F> {
    rows: &'a [SparseRow<F>],
    index: Vec<u32>,
}
impl<'a, F> Table<'a, F> {
    fn new(rows: &'a [SparseRow<F>], ncols: usize) -> Self {
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
    fn get(&self, col: u32) -> Option<&SparseRow<F>> {
        self.rows.get(self.index[col as usize] as usize)
    }
}

// Known pivots plus one slot per column for rows claimed during elimination.
struct Claimed<'a, F> {
    known: &'a Table<'a, F>,
    slots: Vec<OnceLock<(usize, SparseRow<F>)>>,
}
impl<F: Send + Sync> Pivots<F> for Claimed<'_, F> {
    fn get(&self, col: u32) -> Option<&SparseRow<F>> {
        self.known
            .get(col)
            .or_else(|| self.slots[col as usize].get().map(|(_, r)| r))
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
        normalize_residues(row, self.p);
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

// Reduce `rows` by the monic `pivots` and by each other, then back substitute. Workers
// claim free leading columns concurrently, and a row that loses a claim keeps reducing
// by the winner. The reduced echelon form is unique, so the result does not depend on
// scheduling, though which rows are reported live (independent) may.
fn echelon<F: Clone + Send + Sync, K: Kernel<F>>(
    kernel: &K,
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    let known = Table::new(pivots, ncols);
    let claimed = Claimed {
        known: &known,
        slots: (0..ncols).map(|_| OnceLock::new()).collect(),
    };
    let mut order: Vec<usize> = (0..rows.len())
        .filter(|&i| !rows[i].columns.is_empty())
        .collect();
    order.sort_by_key(|&i| (rows[i].columns[0], rows[i].columns.len()));
    map_rows(
        &order,
        || kernel.scratch(),
        |&i, s| {
            let mut r = kernel.reduce(&rows[i], &claimed, s);
            while let Some(&lead) = r.columns.first() {
                kernel.normalize(&mut r);
                match claimed.slots[lead as usize].set((i, r)) {
                    Ok(()) => return,
                    Err((_, lost)) => r = kernel.reduce(&lost, &claimed, s),
                }
            }
        },
    );
    let (live, new): (Vec<usize>, Vec<SparseRow<F>>) = claimed
        .slots
        .into_iter()
        .filter_map(OnceLock::into_inner)
        .unzip();
    // The column sweep clears every pivot column it meets, fill-in included, so each
    // tail reduces independently against the unreduced echelon form.
    let table = Table::new(&new, ncols);
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
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    echelon(&Generic(ncols), pivots, rows, ncols).0
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
            for (&pc, pv) in piv.columns.iter().zip(&piv.coefficients).skip(1) {
                let pc = pc as usize;
                buf[pc] = buf[pc].clone() - c.clone() * pv.clone();
            }
            hi = hi.max(span(&piv.columns).1);
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
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    echelonize_modular_traced(pivots, rows, ncols).0
}

fn echelonize_modular_traced<F: ModularField + Send + Sync>(
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    match modulus(pivots, rows) {
        None => (
            echelonize_generic(pivots, rows, ncols),
            (0..rows.len()).collect(),
        ),
        Some(p) if p <= u64::from(u32::MAX) => echelon_residues::<F, u32>(p, pivots, rows, ncols),
        Some(p) => echelon_residues::<F, u64>(p, pivots, rows, ncols),
    }
}

fn echelon_residues<F: ModularField + Send + Sync, C: Residue>(
    p: u64,
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    let (known, pending) = residues::<F, C>(p, pivots, rows);
    let (reduced, live) = echelon(&Dense { p, ncols }, &known, &pending, ncols);
    (lift(reduced, p), live)
}

fn reduce_rows_modular<F: ModularField + Send + Sync>(
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    match modulus(pivots, rows) {
        None => rows.to_vec(),
        Some(p) if p <= u64::from(u32::MAX) => reduce_residues::<F, u32>(p, pivots, rows, ncols),
        Some(p) => reduce_residues::<F, u64>(p, pivots, rows, ncols),
    }
}

fn reduce_residues<F: ModularField + Send + Sync, C: Residue>(
    p: u64,
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    let (known, pending) = residues::<F, C>(p, pivots, rows);
    let kernel = Dense { p, ncols };
    let table = Table::new(&known, ncols);
    let reduced = map_rows(
        &pending,
        || Kernel::<C>::scratch(&kernel),
        |r, s| kernel.reduce(r, &table, s),
    );
    lift(reduced, p)
}

// The modulus, or `None` when every coefficient is zero and so none can be recovered.
fn modulus<F: ModularField>(pivots: &[SparseRow<F>], rows: &[SparseRow<F>]) -> Option<u64> {
    pivots
        .iter()
        .chain(rows)
        .flat_map(|r| &r.coefficients)
        .map(ModularField::modulus)
        .find(|&m| m != 0)
}

type Residues<C> = (Vec<SparseRow<C>>, Vec<SparseRow<C>>);

// Residues of monic-normalized pivots and of rows.
fn residues<F: ModularField, C: Residue>(
    p: u64,
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
) -> Residues<C> {
    let convert = |r: &SparseRow<F>| SparseRow {
        columns: r.columns.clone(),
        coefficients: r
            .coefficients
            .iter()
            .map(|v| C::new(v.residue_mod(p)))
            .collect(),
    };
    let mut known: Vec<_> = pivots.iter().map(convert).collect();
    for r in &mut known {
        normalize_residues(r, p);
    }
    (known, rows.iter().map(convert).collect())
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

fn normalize_residues<C: Residue>(r: &mut SparseRow<C>, p: u64) {
    let Some(&lead) = r.coefficients.first() else {
        return;
    };
    if lead.into() == 1 {
        return;
    }
    let inv = inv_mod(lead.into(), p);
    for v in &mut r.coefficients {
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
    const LIMIT: u64 = 1 << 63;
    let small = p < (1 << 31);
    let ncols = buf.len();
    // `buf` is all zero on entry and is cleared again while compressing.
    // At most ncols triangular pivot eliminations can contribute to any cell.
    let deferred = small
        && (ncols as u128) * u128::from(p - 1).pow(2) + u128::from(p - 1) <= u128::from(u64::MAX);
    for (&c, &v) in row.columns.iter().zip(&row.coefficients) {
        buf[c as usize] = v.into();
    }
    let (lo, mut hi) = span(&row.columns);
    let mut j = lo;
    while j <= hi && j < ncols {
        if buf[j] == 0 {
            j += 1;
            continue;
        }
        let Some(piv) = pivots.get(j as u32) else {
            j += 1;
            continue;
        };
        let v = std::mem::take(&mut buf[j]) % p;
        if v != 0 {
            let c = p - v;
            hi = hi.max(span(&piv.columns).1);
            let tail = piv.columns[1..].iter().zip(&piv.coefficients[1..]);
            if deferred {
                for (&pc, &pv) in tail {
                    buf[pc as usize] += c * pv.into();
                }
            } else if small {
                for (&pc, &pv) in tail {
                    let acc = buf[pc as usize] + c * pv.into();
                    buf[pc as usize] = if acc >= LIMIT { acc % p } else { acc };
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
        columns: Vec::new(),
        coefficients: Vec::new(),
    };
    for (c, cell) in buf
        .iter_mut()
        .enumerate()
        .take(hi.saturating_add(1))
        .skip(lo)
    {
        let v = std::mem::take(cell) % p;
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
