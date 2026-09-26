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
use std::collections::{HashMap, HashSet, VecDeque};

/// A sparse matrix row with columns in ascending index order (descending monomial order).
#[derive(Debug, Clone)]
pub struct SparseRow<F> {
    pub columns: Vec<usize>,
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
        let plan = state.symbolic_preprocessing(&selected, &order);
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
    finish_basis(basis, canonicalize)
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
    /// A missing pivot or a changed row space invalidates the trace; the caller runs full F4.
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
            let matrix = round.plan.execute(&basis)?;
            let live: Vec<_> = round.live.iter().map(|&i| matrix.rows[i].clone()).collect();
            let rows = Zp::echelonize(&matrix.pivots, &live, round.plan.columns.len());
            // Zero rows never enter the new-row echelonization. Check their dependencies
            // against the completed pivots so unlucky learning primes cannot hide a relation.
            let mut used = vec![false; matrix.rows.len()];
            for &i in &round.live {
                used[i] = true;
            }
            let skipped: Vec<_> = matrix
                .rows
                .iter()
                .zip(used)
                .filter_map(|(r, used)| (!used).then_some(r.clone()))
                .collect();
            if !skipped.is_empty() {
                let mut pivots = matrix.pivots;
                pivots.extend(rows.iter().cloned());
                if !Zp::echelonize(&pivots, &skipped, round.plan.columns.len()).is_empty() {
                    return None;
                }
            }
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
        finish_basis(
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
}

struct State<F> {
    basis: Vec<Polynomial<F>>,
    active: Vec<bool>,
    pairs: Vec<Pair>,
    masks: Vec<u64>,
}

struct Matrix<F> {
    pivots: Vec<SparseRow<F>>,
    rows: Vec<SparseRow<F>>,
}

#[derive(Default)]
struct MonomialTable {
    ids: HashMap<Monomial, u32>,
    monomials: Vec<Monomial>,
}
impl MonomialTable {
    #[allow(clippy::expect_used)] // More than 2^32 matrix columns cannot fit in practical memory.
    fn intern(&mut self, m: Monomial) -> (u32, bool) {
        if let Some(&id) = self.ids.get(&m) {
            return (id, false);
        }
        let id = u32::try_from(self.monomials.len()).expect("F4 matrix exceeds u32 columns");
        self.monomials.push(m.clone());
        self.ids.insert(m, id);
        (id, true)
    }
}

struct PlannedRow {
    basis: usize,
    terms: Vec<Monomial>,
    columns: Vec<usize>,
}
struct MatrixPlan {
    columns: Vec<Monomial>,
    pivots: Vec<PlannedRow>,
    rows: Vec<PlannedRow>,
}
impl MatrixPlan {
    fn execute<F: Field>(&self, basis: &[Polynomial<F>]) -> Option<Matrix<F>> {
        let encode = |row: &PlannedRow| {
            let p = basis.get(row.basis)?;
            let mut coefficients = Vec::with_capacity(p.terms.len());
            let mut columns = Vec::with_capacity(p.terms.len());
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
            rows: self.rows.iter().map(encode).collect::<Option<_>>()?,
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
        let degree = lcm.degree();
        Pair { i, j, lcm, degree }
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
        self.masks.push(divisibility_mask(&lm_h));
        self.basis.push(h);
        self.active.push(true);

        let mut candidates: VecDeque<Pair> = (0..t)
            .filter(|&g| self.active[g])
            .map(|g| self.pair(g, t))
            .collect();
        let mut kept: Vec<Pair> = Vec::new();
        while let Some(p) = candidates.pop_front() {
            let coprime = self.lm(p.i).is_coprime(&lm_h);
            let dominated = |q: &Pair| q.lcm.divides(&p.lcm);
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
        self.pairs.retain(|p| {
            !lm_h.divides(&p.lcm) || lm(p.i).lcm(&lm_h) == p.lcm || lm(p.j).lcm(&lm_h) == p.lcm
        });
        self.pairs.extend(kept);
        for g in 0..t {
            if self.active[g] && lm_h.divides(lm(g)) {
                self.active[g] = false;
            }
        }
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
        (0..self.basis.len())
            .filter(|&i| self.active[i] && self.masks[i] & !mask == 0 && self.lm(i).divides(m))
            .min_by_key(|&i| self.basis[i].terms.len())
    }

    fn symbolic_preprocessing(&self, pairs: &[Pair], order: &MonomialOrder) -> MatrixPlan {
        let mut table = MonomialTable::default();
        let mut multiples: HashSet<(usize, Monomial)> = HashSet::new();
        let mut by_lcm: HashMap<u32, Vec<(usize, Monomial)>> = HashMap::new();
        for p in pairs {
            let (id, _) = table.intern(p.lcm.clone());
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
                let (id, fresh) = table.intern(m * mult);
                if fresh {
                    queue.push(id);
                }
                columns.push(id as usize);
            }
            PlannedRow {
                basis: i,
                terms: self.basis[i].support(),
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
        let mut index = vec![0; ids.len()];
        for (col, &id) in ids.iter().enumerate() {
            index[id] = col;
        }
        let columns = ids
            .into_iter()
            .map(|id| table.monomials[id].clone())
            .collect();
        for row in pivots.iter_mut().chain(&mut rows) {
            for c in &mut row.columns {
                *c = index[*c];
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
        .map(|(c, v)| term(v.clone(), columns[*c].clone()))
        .collect();
    Polynomial {
        terms,
        nvars,
        order: order.clone(),
    }
}

fn pivot_table<F>(pivots: &[SparseRow<F>], ncols: usize) -> Vec<usize> {
    let mut table = vec![usize::MAX; ncols];
    for (k, row) in pivots.iter().enumerate() {
        if let Some(&c) = row.columns.first() {
            table[c] = k;
        }
    }
    table
}

fn echelonize_generic<F: Field>(
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    let table = pivot_table(pivots, ncols);
    let reduce = |row: &SparseRow<F>, pivots: &[SparseRow<F>], table: &[usize]| {
        let mut buf: Vec<F> = vec![F::zero(); ncols];
        let (mut lo, mut hi) = (ncols, 0);
        for (c, v) in row.columns.iter().zip(&row.coefficients) {
            buf[*c] = v.clone();
            lo = lo.min(*c);
            hi = hi.max(*c);
        }
        let mut j = lo;
        while j <= hi && j < ncols {
            if !buf[j].is_zero() && table[j] != usize::MAX {
                let piv = &pivots[table[j]];
                let c = std::mem::replace(&mut buf[j], F::zero());
                for (pc, pv) in piv.columns.iter().zip(&piv.coefficients).skip(1) {
                    buf[*pc] = buf[*pc].clone() - c.clone() * pv.clone();
                }
                hi = hi.max(*piv.columns.last().unwrap_or(&0));
            }
            j += 1;
        }
        compress(buf, lo)
    };
    let reduced: Vec<SparseRow<F>> = rows.iter().map(|r| reduce(r, pivots, &table)).collect();
    echelon_new_rows(
        reduced,
        ncols,
        |r, p, t| reduce(r, p, t),
        |r| {
            let inv = r.coefficients[0]
                .inverse()
                .unwrap_or_else(|| unreachable!("nonzero coefficient"));
            for v in &mut r.coefficients {
                *v = v.clone() * inv.clone();
            }
        },
    )
}

fn compress<F: Field>(buf: Vec<F>, lo: usize) -> SparseRow<F> {
    let mut row = SparseRow {
        columns: Vec::new(),
        coefficients: Vec::new(),
    };
    for (c, v) in buf.into_iter().enumerate().skip(lo) {
        if !v.is_zero() {
            row.columns.push(c);
            row.coefficients.push(v);
        }
    }
    row
}

// Echelonize rows already reduced against the main pivots, then back substitute.
fn echelon_new_rows<F: Clone>(
    rows: Vec<SparseRow<F>>,
    ncols: usize,
    reduce: impl Fn(&SparseRow<F>, &[SparseRow<F>], &[usize]) -> SparseRow<F>,
    normalize: impl Fn(&mut SparseRow<F>),
) -> Vec<SparseRow<F>> {
    echelon_new_rows_traced(rows, ncols, reduce, normalize).0
}

fn echelon_new_rows_traced<F: Clone>(
    rows: Vec<SparseRow<F>>,
    ncols: usize,
    reduce: impl Fn(&SparseRow<F>, &[SparseRow<F>], &[usize]) -> SparseRow<F>,
    normalize: impl Fn(&mut SparseRow<F>),
) -> (Vec<SparseRow<F>>, Vec<usize>) {
    let mut rows: Vec<_> = rows
        .into_iter()
        .enumerate()
        .filter(|(_, r)| !r.columns.is_empty())
        .collect();
    rows.sort_by_key(|(_, r)| (r.columns[0], r.columns.len()));
    let mut live = Vec::new();
    let mut new: Vec<SparseRow<F>> = Vec::new();
    let mut table = vec![usize::MAX; ncols];
    for (source, row) in rows {
        let mut r = reduce(&row, &new, &table);
        if r.columns.is_empty() {
            continue;
        }
        live.push(source);
        normalize(&mut r);
        table[r.columns[0]] = new.len();
        new.push(r);
    }
    for k in (0..new.len()).rev() {
        let lead = new[k].columns[0];
        table[lead] = usize::MAX;
        let mut r = reduce(&new[k], &new, &table);
        table[lead] = k;
        if r.columns.first() == Some(&lead) {
            normalize(&mut r);
        }
        new[k] = r;
    }
    (new, live)
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
    let modulus = pivots
        .iter()
        .chain(rows)
        .flat_map(|r| &r.coefficients)
        .map(ModularField::modulus)
        .find(|&m| m != 0);
    let Some(p) = modulus else {
        return (
            echelonize_generic(pivots, rows, ncols),
            (0..rows.len()).collect(),
        );
    };
    let to_u64 = |r: &SparseRow<F>| SparseRow {
        columns: r.columns.clone(),
        coefficients: r.coefficients.iter().map(|v| v.residue_mod(p)).collect(),
    };
    let mut pivots_u64: Vec<SparseRow<u64>> = pivots.iter().map(to_u64).collect();
    for r in &mut pivots_u64 {
        normalize_u64(r, p);
    }
    let rows_u64: Vec<SparseRow<u64>> = rows.iter().map(to_u64).collect();
    let table = pivot_table(&pivots_u64, ncols);
    // Known pivots and new rows occupy separate blocks. Each worker reuses one
    // dense buffer while reducing a block against the immutable known pivots.
    let reduce_block = |block: &[SparseRow<u64>]| {
        let mut scratch = vec![0; ncols];
        block
            .iter()
            .map(|r| reduce_dense_u64(r, &pivots_u64, &table, &mut scratch, p))
            .collect::<Vec<_>>()
    };
    #[cfg(feature = "parallel")]
    let reduced: Vec<SparseRow<u64>> = if rows_u64.len() >= 64 {
        rows_u64
            .par_chunks(32)
            .flat_map_iter(reduce_block)
            .collect()
    } else {
        reduce_block(&rows_u64)
    };
    #[cfg(not(feature = "parallel"))]
    let reduced: Vec<SparseRow<u64>> = reduce_block(&rows_u64);
    let scratch = std::cell::RefCell::new(vec![0; ncols]);
    let reduce = |r: &SparseRow<u64>, piv: &[SparseRow<u64>], t: &[usize]| {
        reduce_dense_u64(r, piv, t, &mut scratch.borrow_mut(), p)
    };

    let (reduced, live) = echelon_new_rows_traced(reduced, ncols, reduce, |r| normalize_u64(r, p));
    let reduced = reduced
        .into_iter()
        .map(|r| SparseRow {
            columns: r.columns,
            coefficients: r
                .coefficients
                .into_iter()
                .map(|v| F::from_residue(v, p))
                .collect(),
        })
        .collect();
    (reduced, live)
}

fn normalize_u64(r: &mut SparseRow<u64>, p: u64) {
    let Some(&lead) = r.coefficients.first() else {
        return;
    };
    if lead == 1 {
        return;
    }
    let inv = inv_mod(lead, p);
    for v in &mut r.coefficients {
        *v = mul_mod(*v, inv, p);
    }
}

// Monagan-Pearce dense reduction: sums are reduced lazily while they fit in 64 bits.
fn reduce_dense_u64(
    row: &SparseRow<u64>,
    pivots: &[SparseRow<u64>],
    table: &[usize],
    buf: &mut [u64],
    p: u64,
) -> SparseRow<u64> {
    const LIMIT: u64 = 1 << 63;
    let small = p < (1 << 31);
    let ncols = buf.len();
    buf.fill(0);
    // At most ncols triangular pivot eliminations can contribute to any cell.
    let deferred = small
        && (ncols as u128) * u128::from(p - 1).pow(2) + u128::from(p - 1) <= u128::from(u64::MAX);
    let (mut lo, mut hi) = (ncols, 0);
    for (c, v) in row.columns.iter().zip(&row.coefficients) {
        buf[*c] = *v;
        lo = lo.min(*c);
        hi = hi.max(*c);
    }
    let mut j = lo;
    while j <= hi && j < ncols {
        let v = buf[j] % p;
        buf[j] = v;
        if v != 0 && table[j] != usize::MAX {
            buf[j] = 0;
            let piv = &pivots[table[j]];
            let c = p - v;
            hi = hi.max(*piv.columns.last().unwrap_or(&0));
            if deferred {
                for (pc, pv) in piv.columns.iter().zip(&piv.coefficients).skip(1) {
                    buf[*pc] += c * pv;
                }
            } else if small {
                for (pc, pv) in piv.columns.iter().zip(&piv.coefficients).skip(1) {
                    let acc = buf[*pc] + c * pv;
                    buf[*pc] = if acc >= LIMIT { acc % p } else { acc };
                }
            } else {
                for (pc, pv) in piv.columns.iter().zip(&piv.coefficients).skip(1) {
                    buf[*pc] = add_mod(buf[*pc], mul_mod(c, *pv, p), p);
                }
            }
        }
        j += 1;
    }
    let mut out = SparseRow {
        columns: Vec::new(),
        coefficients: Vec::new(),
    };
    for (c, v) in buf
        .iter()
        .copied()
        .enumerate()
        .take(hi.saturating_add(1))
        .skip(lo)
    {
        let v = v % p;
        if v != 0 {
            out.columns.push(c);
            out.coefficients.push(v);
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
    fn unlucky_learning_prime_cannot_hide_a_zero_row() {
        let system = "x + y; x - y";
        let (_, trace) = learn(input(2, system)).expect("learn");
        assert!(
            trace
                .rounds
                .iter()
                .any(|r| r.live.len() < r.plan.rows.len())
        );
        assert!(trace.replay(input(3, system)).is_none());
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
                    let columns: Vec<_> = (i..ncols).collect();
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
                    columns: (0..ncols).collect(),
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
