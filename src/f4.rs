//! Sparse F4-style Groebner basis computation over prime fields.
//!
//! This module batches critical pairs by minimal degree, performs symbolic preprocessing, encodes
//! the resulting Macaulay rows into a sparse matrix, row-reduces over [`PrimeField`], and decodes
//! new basis elements. The implementation is intentionally separate from the Buchberger path so it
//! can be tested against [`crate::groebner_basis`] as a reference implementation.
//!
//! # Example
//! ```
//! use groebner::{groebner_basis_f4_mod, MonomialOrder, PolynomialRing, PrimeField};
//!
//! type F32003 = PrimeField<32003>;
//!
//! let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::Lex)?;
//! let f1 = ring.parse("x^2 - y")?;
//! let f2 = ring.parse("x*y - 1")?;
//! let basis = groebner_basis_f4_mod(vec![f1, f2], F32003::modulus(), ring.order())?;
//!
//! assert!(!basis.is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::field::Field;
use crate::finite_field::PrimeField;
use crate::groebner::GroebnerError;
use crate::monomial::{Monomial, MonomialOrder};
use crate::polynomial::{Polynomial, Term};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
struct F4Pair {
    i: usize,
    j: usize,
    lcm: Monomial,
    degree: u32,
}

impl F4Pair {
    fn new<const P: u32>(
        i: usize,
        j: usize,
        basis: &[Polynomial<PrimeField<P>>],
    ) -> Result<Self, GroebnerError> {
        Self::new_from_polynomials(i, j, &basis[i], &basis[j])
    }

    fn new_from_polynomials<const P: u32>(
        i: usize,
        j: usize,
        poly_i: &Polynomial<PrimeField<P>>,
        poly_j: &Polynomial<PrimeField<P>>,
    ) -> Result<Self, GroebnerError> {
        let lm_i = poly_i
            .leading_monomial()
            .ok_or(GroebnerError::NoLeadingMonomial(i))?;
        let lm_j = poly_j
            .leading_monomial()
            .ok_or(GroebnerError::NoLeadingMonomial(j))?;
        let lcm = lm_i.lcm(lm_j);
        let degree = lcm.degree();
        Ok(Self { i, j, lcm, degree })
    }
}

#[derive(Debug, Clone)]
struct SparseRow<F> {
    columns: Vec<usize>,
    coefficients: Vec<F>,
}

impl<F: Field> SparseRow<F> {
    fn new(columns: Vec<usize>, coefficients: Vec<F>) -> Self {
        let mut row = Self {
            columns,
            coefficients,
        };
        row.remove_zeros();
        row
    }

    fn is_zero(&self) -> bool {
        self.columns.is_empty()
    }

    fn leading_column(&self) -> Option<usize> {
        self.columns.first().copied()
    }

    fn leading_coefficient(&self) -> Option<&F> {
        self.coefficients.first()
    }

    fn normalize(&mut self) -> Result<(), GroebnerError> {
        let Some(leading) = self.leading_coefficient() else {
            return Ok(());
        };
        let inverse = leading.inverse().ok_or({
            GroebnerError::Polynomial(crate::polynomial::PolynomialError::DivisionByZero)
        })?;
        for coefficient in &mut self.coefficients {
            *coefficient = coefficient.multiply(&inverse);
        }
        self.remove_zeros();
        Ok(())
    }

    fn subtract_scaled(&mut self, other: &Self, scale: &F) {
        if scale.is_zero() || other.is_zero() {
            return;
        }

        let mut columns = Vec::with_capacity(self.columns.len() + other.columns.len());
        let mut coefficients =
            Vec::with_capacity(self.coefficients.len() + other.coefficients.len());
        let mut i = 0;
        let mut j = 0;

        while i < self.columns.len() && j < other.columns.len() {
            match self.columns[i].cmp(&other.columns[j]) {
                std::cmp::Ordering::Less => {
                    columns.push(self.columns[i]);
                    coefficients.push(self.coefficients[i].clone());
                    i += 1;
                }
                std::cmp::Ordering::Greater => {
                    let scaled = other.coefficients[j].multiply(scale).negate();
                    if !scaled.is_zero() {
                        columns.push(other.columns[j]);
                        coefficients.push(scaled);
                    }
                    j += 1;
                }
                std::cmp::Ordering::Equal => {
                    let scaled = other.coefficients[j].multiply(scale);
                    let coefficient = self.coefficients[i].subtract(&scaled);
                    if !coefficient.is_zero() {
                        columns.push(self.columns[i]);
                        coefficients.push(coefficient);
                    }
                    i += 1;
                    j += 1;
                }
            }
        }

        while i < self.columns.len() {
            columns.push(self.columns[i]);
            coefficients.push(self.coefficients[i].clone());
            i += 1;
        }

        while j < other.columns.len() {
            let scaled = other.coefficients[j].multiply(scale).negate();
            if !scaled.is_zero() {
                columns.push(other.columns[j]);
                coefficients.push(scaled);
            }
            j += 1;
        }

        self.columns = columns;
        self.coefficients = coefficients;
    }

    fn remove_zeros(&mut self) {
        let mut columns = Vec::with_capacity(self.columns.len());
        let mut coefficients = Vec::with_capacity(self.coefficients.len());
        for (column, coefficient) in self.columns.iter().copied().zip(&self.coefficients) {
            if !coefficient.is_zero() {
                columns.push(column);
                coefficients.push(coefficient.clone());
            }
        }
        self.columns = columns;
        self.coefficients = coefficients;
    }
}

/// Compute a Groebner basis using a sparse F4-style algorithm over `GF(P)`.
///
/// The runtime `p` argument is checked against the const generic modulus so callers can keep the
/// API shape close to command-line F4 tools while the coefficient type remains statically known.
pub fn groebner_basis_f4_mod<const P: u32>(
    polynomials: Vec<Polynomial<PrimeField<P>>>,
    p: u32,
    order: MonomialOrder,
) -> Result<Vec<Polynomial<PrimeField<P>>>, GroebnerError> {
    if p != P {
        return Err(GroebnerError::InvalidPrimeField {
            expected: P,
            actual: p,
        });
    }
    if polynomials.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }

    let nvars = polynomials[0].nvars;
    let mut basis: Vec<_> = polynomials
        .into_iter()
        .filter(|polynomial| !polynomial.is_zero())
        .map(|polynomial| polynomial.make_monic())
        .collect();
    if basis.is_empty() {
        return Err(GroebnerError::EmptyInput);
    }

    let mut pairs = initial_pairs(&basis)?;
    while !pairs.is_empty() {
        let selected = select_min_degree_pairs(&mut pairs);
        let rows = symbolic_preprocessing(&basis, &selected)?;
        if rows.is_empty() {
            continue;
        }

        let reduced_rows = reduce_sparse_matrix(&rows, order)?;
        let mut new_polynomials = decode_new_polynomials(reduced_rows, &basis)?;
        if new_polynomials.is_empty() {
            continue;
        }

        for polynomial in new_polynomials.drain(..) {
            let reduced = polynomial.reduce(&basis).map_err(GroebnerError::from)?;
            if reduced.is_zero() {
                continue;
            }
            let monic = reduced.make_monic();
            let new_index = basis.len();
            for (i, existing) in basis.iter().enumerate() {
                pairs.push(F4Pair::new_from_polynomials(
                    i, new_index, existing, &monic,
                )?);
            }
            basis.push(monic);
        }
    }

    canonicalize_basis(basis, nvars, order)
}

fn initial_pairs<const P: u32>(
    basis: &[Polynomial<PrimeField<P>>],
) -> Result<Vec<F4Pair>, GroebnerError> {
    let mut pairs = Vec::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            pairs.push(F4Pair::new(i, j, basis)?);
        }
    }
    Ok(pairs)
}

fn select_min_degree_pairs(pairs: &mut Vec<F4Pair>) -> Vec<F4Pair> {
    let Some(min_degree) = pairs.iter().map(|pair| pair.degree).min() else {
        return Vec::new();
    };

    let mut selected = Vec::new();
    let mut remaining = Vec::new();
    for pair in pairs.drain(..) {
        if pair.degree == min_degree {
            selected.push(pair);
        } else {
            remaining.push(pair);
        }
    }
    *pairs = remaining;
    selected
}

fn symbolic_preprocessing<const P: u32>(
    basis: &[Polynomial<PrimeField<P>>],
    pairs: &[F4Pair],
) -> Result<Vec<Polynomial<PrimeField<P>>>, GroebnerError> {
    let mut rows = Vec::new();
    let mut row_keys = HashSet::new();
    let mut pending_monomials = Vec::new();
    let mut seen_monomials = HashSet::new();

    for pair in pairs {
        if pair.i >= basis.len() || pair.j >= basis.len() {
            continue;
        }
        let Some(lm_i) = basis[pair.i].leading_monomial() else {
            return Err(GroebnerError::NoLeadingMonomial(pair.i));
        };
        let Some(lm_j) = basis[pair.j].leading_monomial() else {
            return Err(GroebnerError::NoLeadingMonomial(pair.j));
        };
        if pair.lcm == lm_i.multiply(lm_j) {
            continue;
        }

        let multiplier_i = pair
            .lcm
            .divide(lm_i)
            .ok_or(crate::polynomial::PolynomialError::DivisionFailed)?;
        let multiplier_j = pair
            .lcm
            .divide(lm_j)
            .ok_or(crate::polynomial::PolynomialError::DivisionFailed)?;
        add_symbolic_row(
            basis,
            pair.i,
            &multiplier_i,
            &mut rows,
            &mut row_keys,
            &mut pending_monomials,
            &mut seen_monomials,
        );
        add_symbolic_row(
            basis,
            pair.j,
            &multiplier_j,
            &mut rows,
            &mut row_keys,
            &mut pending_monomials,
            &mut seen_monomials,
        );
    }

    while let Some(monomial) = pending_monomials.pop() {
        let Some((basis_index, leading)) = find_reducer(basis, &monomial) else {
            continue;
        };
        let multiplier = monomial
            .divide(leading)
            .ok_or(crate::polynomial::PolynomialError::DivisionFailed)?;
        add_symbolic_row(
            basis,
            basis_index,
            &multiplier,
            &mut rows,
            &mut row_keys,
            &mut pending_monomials,
            &mut seen_monomials,
        );
    }

    Ok(rows)
}

fn add_symbolic_row<const P: u32>(
    basis: &[Polynomial<PrimeField<P>>],
    basis_index: usize,
    multiplier: &Monomial,
    rows: &mut Vec<Polynomial<PrimeField<P>>>,
    row_keys: &mut HashSet<(usize, Monomial)>,
    pending_monomials: &mut Vec<Monomial>,
    seen_monomials: &mut HashSet<Monomial>,
) {
    if !row_keys.insert((basis_index, multiplier.clone())) {
        return;
    }

    let row = basis[basis_index].multiply_monomial(multiplier);
    for term in &row.terms {
        if seen_monomials.insert(term.monomial.clone()) {
            pending_monomials.push(term.monomial.clone());
        }
    }
    rows.push(row);
}

fn find_reducer<'a, const P: u32>(
    basis: &'a [Polynomial<PrimeField<P>>],
    monomial: &Monomial,
) -> Option<(usize, &'a Monomial)> {
    basis.iter().enumerate().find_map(|(index, polynomial)| {
        polynomial
            .leading_monomial()
            .filter(|leading| leading.divides(monomial))
            .map(|leading| (index, leading))
    })
}

fn reduce_sparse_matrix<const P: u32>(
    rows: &[Polynomial<PrimeField<P>>],
    order: MonomialOrder,
) -> Result<Vec<Polynomial<PrimeField<P>>>, GroebnerError> {
    let Some(first) = rows.first() else {
        return Ok(Vec::new());
    };
    let nvars = first.nvars;
    let columns = collect_columns(rows, order);
    let column_indices: HashMap<_, _> = columns
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, monomial)| (monomial, index))
        .collect();

    let mut encoded: Vec<_> = rows
        .iter()
        .map(|row| encode_row(row, &column_indices))
        .collect();
    encoded.sort_by_key(|row| row.leading_column().unwrap_or(usize::MAX));

    let mut pivots: HashMap<usize, SparseRow<PrimeField<P>>> = HashMap::new();
    let mut pivot_order = Vec::new();
    for mut row in encoded {
        reduce_row(&mut row, &pivots)?;
        if row.is_zero() {
            continue;
        }
        row.normalize()?;
        if let Some(column) = row.leading_column() {
            pivots.insert(column, row);
            pivot_order.push(column);
        }
    }

    pivot_order.sort_unstable();
    let mut reduced = Vec::new();
    for column in pivot_order {
        if let Some(row) = pivots.remove(&column) {
            let polynomial = decode_row(row, &columns, nvars, order);
            if !polynomial.is_zero() {
                reduced.push(polynomial);
            }
        }
    }
    Ok(reduced)
}

fn collect_columns<const P: u32>(
    rows: &[Polynomial<PrimeField<P>>],
    order: MonomialOrder,
) -> Vec<Monomial> {
    let mut seen = HashSet::new();
    let mut columns = Vec::new();
    for row in rows {
        for term in &row.terms {
            if seen.insert(term.monomial.clone()) {
                columns.push(term.monomial.clone());
            }
        }
    }
    columns.sort_by(|a, b| b.compare(a, order));
    columns
}

fn encode_row<const P: u32>(
    row: &Polynomial<PrimeField<P>>,
    column_indices: &HashMap<Monomial, usize>,
) -> SparseRow<PrimeField<P>> {
    let mut entries: Vec<_> = row
        .terms
        .iter()
        .filter_map(|term| {
            column_indices
                .get(&term.monomial)
                .copied()
                .map(|column| (column, term.coefficient))
        })
        .collect();
    entries.sort_by_key(|(column, _)| *column);

    let (columns, coefficients): (Vec<_>, Vec<_>) = entries.into_iter().unzip();
    SparseRow::new(columns, coefficients)
}

fn reduce_row<const P: u32>(
    row: &mut SparseRow<PrimeField<P>>,
    pivots: &HashMap<usize, SparseRow<PrimeField<P>>>,
) -> Result<(), GroebnerError> {
    while let Some(column) = row.leading_column() {
        let Some(pivot) = pivots.get(&column) else {
            break;
        };
        let coefficient = row
            .leading_coefficient()
            .ok_or({
                GroebnerError::Polynomial(crate::polynomial::PolynomialError::NoLeadingCoefficient)
            })?
            .to_owned();
        row.subtract_scaled(pivot, &coefficient);
    }
    Ok(())
}

fn decode_row<const P: u32>(
    row: SparseRow<PrimeField<P>>,
    columns: &[Monomial],
    nvars: usize,
    order: MonomialOrder,
) -> Polynomial<PrimeField<P>> {
    let terms = row
        .columns
        .into_iter()
        .zip(row.coefficients)
        .filter_map(|(column, coefficient)| {
            columns
                .get(column)
                .cloned()
                .map(|monomial| Term::new(coefficient, monomial))
        })
        .collect();
    Polynomial::new(terms, nvars, order)
}

fn decode_new_polynomials<const P: u32>(
    rows: Vec<Polynomial<PrimeField<P>>>,
    basis: &[Polynomial<PrimeField<P>>],
) -> Result<Vec<Polynomial<PrimeField<P>>>, GroebnerError> {
    let mut existing_leading = HashSet::new();
    for polynomial in basis {
        if let Some(leading) = polynomial.leading_monomial() {
            existing_leading.insert(leading.clone());
        }
    }

    let mut decoded = Vec::new();
    let mut new_leading = HashSet::new();
    for row in rows {
        let monic = row.make_monic();
        let Some(leading) = monic.leading_monomial() else {
            continue;
        };
        if existing_leading.contains(leading) || !new_leading.insert(leading.clone()) {
            continue;
        }
        decoded.push(monic);
    }
    Ok(decoded)
}

fn canonicalize_basis<const P: u32>(
    mut basis: Vec<Polynomial<PrimeField<P>>>,
    nvars: usize,
    order: MonomialOrder,
) -> Result<Vec<Polynomial<PrimeField<P>>>, GroebnerError> {
    basis.retain(|polynomial| !polynomial.is_zero());
    for polynomial in &mut basis {
        *polynomial = polynomial.make_monic();
    }

    let mut i = 0;
    while i < basis.len() {
        let Some(lm_i) = basis[i].leading_monomial().cloned() else {
            basis.remove(i);
            continue;
        };
        let redundant = basis.iter().enumerate().any(|(j, polynomial)| {
            i != j
                && polynomial
                    .leading_monomial()
                    .is_some_and(|lm_j| lm_j.divides(&lm_i))
        });
        if redundant {
            basis.remove(i);
        } else {
            i += 1;
        }
    }

    let mut reduced = Vec::new();
    for i in 0..basis.len() {
        let others: Vec<_> = basis
            .iter()
            .enumerate()
            .filter_map(|(j, polynomial)| {
                if i == j {
                    None
                } else {
                    Some(polynomial.clone())
                }
            })
            .collect();
        let polynomial = basis[i]
            .reduce(&others)
            .map_err(GroebnerError::from)?
            .make_monic();
        if !polynomial.is_zero() {
            reduced.push(polynomial);
        }
    }

    reduced.sort_by(|a, b| match (a.leading_monomial(), b.leading_monomial()) {
        (Some(ma), Some(mb)) => mb.compare(ma, order),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });

    if reduced.is_empty() {
        Ok(vec![Polynomial::zero(nvars, order)])
    } else {
        Ok(reduced)
    }
}
