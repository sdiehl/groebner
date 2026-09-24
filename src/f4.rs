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

use crate::field::{Field, ModularField};
use crate::finite_field::{PrimeField, Zp};
use crate::groebner::{finish_basis, prepare_input, GroebnerError};
use crate::monomial::{Monomial, MonomialOrder};
use crate::polynomial::{Polynomial, Term};
use num_rational::BigRational;
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

impl F4Field for BigRational {}

impl<const P: u32> F4Field for PrimeField<P> {
    fn echelonize(
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_modular(pivots, rows, ncols)
    }
}

impl F4Field for Zp {
    fn echelonize(
        pivots: &[SparseRow<Self>],
        rows: &[SparseRow<Self>],
        ncols: usize,
    ) -> Vec<SparseRow<Self>> {
        echelonize_modular(pivots, rows, ncols)
    }
}

/// Compute a Groebner basis with F4 under the polynomials' own order.
pub fn groebner_basis_f4<F: F4Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    let input = prepare_input(polynomials)?;
    let order = input[0].order.clone();
    let nvars = input[0].nvars;
    let mut state = State {
        basis: Vec::new(),
        active: Vec::new(),
        pairs: Vec::new(),
    };
    for poly in input {
        state.update(poly);
    }
    while !state.pairs.is_empty() {
        let selected = state.select();
        let matrix = state.symbolic_preprocessing(&selected, &order);
        let new_rows = F::echelonize(&matrix.pivots, &matrix.rows, matrix.columns.len());
        let mut new_polys: Vec<Polynomial<F>> = new_rows
            .into_iter()
            .map(|row| decode(&row, &matrix.columns, nvars, &order))
            .collect();
        new_polys.sort_by(|a, b| crate::groebner::compare_leading(a, b, &order));
        for poly in new_polys {
            state.update(poly);
        }
    }
    let basis = state
        .basis
        .into_iter()
        .zip(state.active)
        .filter_map(|(poly, active)| active.then_some(poly))
        .collect();
    finish_basis(basis, canonicalize)
}

/// Compatibility wrapper around [`groebner_basis_f4`] for a const-modulus field.
#[deprecated(since = "0.3.0", note = "use groebner_basis_f4")]
pub fn groebner_basis_f4_mod<const P: u32>(
    polynomials: Vec<Polynomial<PrimeField<P>>>,
    modulus: u32,
    _order: MonomialOrder,
) -> Result<Vec<Polynomial<PrimeField<P>>>, GroebnerError> {
    if modulus != P {
        return Err(GroebnerError::InvalidPrimeField {
            expected: P,
            actual: modulus,
        });
    }
    groebner_basis_f4(polynomials, true)
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
}

struct Matrix<F> {
    columns: Vec<Monomial>,
    pivots: Vec<SparseRow<F>>,
    rows: Vec<SparseRow<F>>,
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
        let h = if h.leading_coefficient().is_some_and(Field::is_one) {
            h
        } else {
            h.make_monic()
        };
        let Some(lm_h) = h.leading_monomial().cloned() else {
            return;
        };
        let t = self.basis.len();
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
        (0..self.basis.len())
            .filter(|&i| self.active[i] && self.lm(i).divides(m))
            .min_by_key(|&i| self.basis[i].terms.len())
    }

    fn symbolic_preprocessing(&self, pairs: &[Pair], order: &MonomialOrder) -> Matrix<F> {
        let mut multiples: HashSet<(usize, Monomial)> = HashSet::new();
        let mut by_lcm: HashMap<Monomial, Vec<(usize, Monomial)>> = HashMap::new();
        for p in pairs {
            for i in [p.i, p.j] {
                let Some(mult) = p.lcm.divide(self.lm(i)) else {
                    continue;
                };
                if multiples.insert((i, mult.clone())) {
                    by_lcm.entry(p.lcm.clone()).or_default().push((i, mult));
                }
            }
        }

        let mut pivots: Vec<(usize, Monomial)> = Vec::new();
        let mut rows: Vec<(usize, Monomial)> = Vec::new();
        let mut seen: HashSet<Monomial> = HashSet::new();
        let mut queue: Vec<Monomial> = Vec::new();
        let visit =
            |i: usize, mult: &Monomial, seen: &mut HashSet<Monomial>, queue: &mut Vec<Monomial>| {
                for t in &self.basis[i].terms {
                    let m = t.monomial.multiply(mult);
                    if seen.insert(m.clone()) {
                        queue.push(m);
                    }
                }
            };
        for (lcm, mut group) in by_lcm {
            group.sort_by_key(|(i, _)| self.basis[*i].terms.len());
            seen.insert(lcm);
            for (k, (i, mult)) in group.into_iter().enumerate() {
                visit(i, &mult, &mut seen, &mut queue);
                if k == 0 {
                    pivots.push((i, mult));
                } else {
                    rows.push((i, mult));
                }
            }
        }
        let mut pivot_leads: HashSet<Monomial> = pivots
            .iter()
            .map(|(i, mult)| self.lm(*i).multiply(mult))
            .collect();
        while let Some(m) = queue.pop() {
            if pivot_leads.contains(&m) {
                continue;
            }
            if let Some(i) = self.reducer(&m) {
                let Some(mult) = m.divide(self.lm(i)) else {
                    continue;
                };
                visit(i, &mult, &mut seen, &mut queue);
                pivot_leads.insert(m);
                pivots.push((i, mult));
            }
        }

        let mut columns: Vec<Monomial> = seen.into_iter().collect();
        columns.sort_by(|a, b| b.compare(a, order));
        let index: HashMap<&Monomial, usize> =
            columns.iter().enumerate().map(|(k, m)| (m, k)).collect();
        let encode = |(i, mult): &(usize, Monomial)| {
            let poly = &self.basis[*i];
            let mut cols = Vec::with_capacity(poly.terms.len());
            let mut coeffs = Vec::with_capacity(poly.terms.len());
            for t in &poly.terms {
                cols.push(index[&t.monomial.multiply(mult)]);
                coeffs.push(t.coefficient.clone());
            }
            SparseRow {
                columns: cols,
                coefficients: coeffs,
            }
        };
        let mut pivots: Vec<SparseRow<F>> = pivots.iter().map(encode).collect();
        pivots.sort_by_key(|r| r.columns[0]);
        let rows = rows.iter().map(encode).collect();
        Matrix {
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
        .map(|(c, v)| Term::new(v.clone(), columns[*c].clone()))
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
                    buf[*pc] = buf[*pc].subtract(&c.multiply(pv));
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
                *v = v.multiply(&inv);
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
    mut rows: Vec<SparseRow<F>>,
    ncols: usize,
    reduce: impl Fn(&SparseRow<F>, &[SparseRow<F>], &[usize]) -> SparseRow<F>,
    normalize: impl Fn(&mut SparseRow<F>),
) -> Vec<SparseRow<F>> {
    rows.retain(|r| !r.columns.is_empty());
    rows.sort_by_key(|r| (r.columns[0], r.columns.len()));
    let mut new: Vec<SparseRow<F>> = Vec::new();
    let mut table = vec![usize::MAX; ncols];
    for row in rows {
        let mut r = reduce(&row, &new, &table);
        if r.columns.is_empty() {
            continue;
        }
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
    new
}

fn echelonize_modular<F: ModularField>(
    pivots: &[SparseRow<F>],
    rows: &[SparseRow<F>],
    ncols: usize,
) -> Vec<SparseRow<F>> {
    let modulus = pivots
        .iter()
        .chain(rows)
        .flat_map(|r| &r.coefficients)
        .map(ModularField::modulus)
        .find(|&m| m != 0);
    let Some(p) = modulus else {
        return echelonize_generic(pivots, rows, ncols);
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
    let reduce = |r: &SparseRow<u64>, piv: &[SparseRow<u64>], t: &[usize]| {
        reduce_dense_u64(r, piv, t, ncols, p)
    };

    #[cfg(feature = "parallel")]
    let reduced: Vec<SparseRow<u64>> = rows_u64
        .par_iter()
        .map(|r| reduce(r, &pivots_u64, &table))
        .collect();
    #[cfg(not(feature = "parallel"))]
    let reduced: Vec<SparseRow<u64>> = rows_u64
        .iter()
        .map(|r| reduce(r, &pivots_u64, &table))
        .collect();

    echelon_new_rows(reduced, ncols, reduce, |r| normalize_u64(r, p))
        .into_iter()
        .map(|r| SparseRow {
            columns: r.columns,
            coefficients: r
                .coefficients
                .into_iter()
                .map(|v| F::from_residue(v, p))
                .collect(),
        })
        .collect()
}

fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(p)) as u64
}

fn inv_mod(a: u64, p: u64) -> u64 {
    let (mut t, mut new_t) = (0i128, 1i128);
    let (mut r, mut new_r) = (i128::from(p), i128::from(a));
    while new_r != 0 {
        let q = r / new_r;
        (t, new_t) = (new_t, t - q * new_t);
        (r, new_r) = (new_r, r - q * new_r);
    }
    t.rem_euclid(i128::from(p)) as u64
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
    ncols: usize,
    p: u64,
) -> SparseRow<u64> {
    const LIMIT: u64 = 1 << 63;
    let small = p < (1 << 31);
    let mut buf = vec![0u64; ncols];
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
            if small {
                for (pc, pv) in piv.columns.iter().zip(&piv.coefficients).skip(1) {
                    let acc = buf[*pc] + c * pv;
                    buf[*pc] = if acc >= LIMIT { acc % p } else { acc };
                }
            } else {
                for (pc, pv) in piv.columns.iter().zip(&piv.coefficients).skip(1) {
                    let acc = buf[*pc] + mul_mod(c, *pv, p);
                    buf[*pc] = if acc >= p { acc - p } else { acc };
                }
            }
        }
        j += 1;
    }
    let mut out = SparseRow {
        columns: Vec::new(),
        coefficients: Vec::new(),
    };
    for (c, v) in buf.into_iter().enumerate().skip(lo) {
        let v = v % p;
        if v != 0 {
            out.columns.push(c);
            out.coefficients.push(v);
        }
    }
    out
}
