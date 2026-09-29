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

mod eliminate;
mod finish;
mod monomials;
mod symbolic;
mod trace;

use crate::finite_field::{PrimeField, Zp};
use crate::groebner::{GroebnerError, prepare_input};
use crate::monomial::{Monomial, MonomialOrder};
use crate::polynomial::{Polynomial, term};
use crate::{Field, par};
use eliminate::{
    Row, Table, echelonize_generic, echelonize_modular, echelonize_modular_traced, reduce_generic,
    reduce_rows_modular,
};
use finish::{FinishPlan, finish};
use monomials::pack;
use num_rational::BigRational;
use rustc_hash::FxHashMap as HashMap;
use symbolic::{Element, State};
use trace::TraceRound;
pub(crate) use trace::{F4Trace, PRIMES, learn};

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

fn compute<F: F4Field>(
    polynomials: Vec<Polynomial<F>>,
    canonicalize: bool,
    trace: Option<&mut F4Trace>,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    par::on_pool(|| compute_on(polynomials, canonicalize, trace))
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
    let mut state = State::new(nvars);
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
    if let Some(t) = trace.as_deref_mut() {
        t.active = state.active.clone();
    }
    // Terms share one monomial per distinct exponent vector, as decoded rows would.
    let mut cache = HashMap::default();
    let mut basis: Vec<_> = state
        .basis
        .into_iter()
        .zip(state.active)
        .filter_map(|(g, active)| active.then(|| g.polynomial(nvars, &order, &mut cache)))
        .collect();
    basis.retain(|p| !p.is_zero());
    match trace {
        Some(t) if canonicalize => {
            let plan = FinishPlan::new(&basis)?;
            let basis = plan.execute(basis);
            t.finish = Some(plan);
            Ok(basis)
        }
        _ => finish(basis, canonicalize),
    }
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
