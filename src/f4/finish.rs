//! Minimization and interreduction of a finished F4 basis with one Macaulay matrix,
//! planned once and replayable at other primes.

use super::eliminate::{Kernel, Row, Table};
use super::monomials::MonomialTable;
use super::trace::{LaneKernel, Lanes, PRIMES};
use super::{F4Field, SparseRow, column_index, decode, reducers};
use crate::finite_field::Zp;
use crate::groebner::{GroebnerError, finish_basis};
use crate::monomial::{Monomial, divisibility_mask};
use crate::polynomial::Polynomial;
use crate::{Field, ModularField, par};

/// Minimize and, when `canonicalize` is set, interreduce with one Macaulay matrix: the tails
/// of the minimal basis are the rows, and symbolic preprocessing supplies the reducers.
pub(super) fn finish<F: F4Field>(
    mut basis: Vec<Polynomial<F>>,
    canonicalize: bool,
) -> Result<Vec<Polynomial<F>>, GroebnerError> {
    if !canonicalize {
        return finish_basis(basis, false);
    }
    basis.retain(|p| !p.is_zero());
    Ok(FinishPlan::new(&basis)?.execute(basis))
}

// The coefficient-independent half of the final interreduction: the minimal elements
// of a basis, and the reducers and columns of their tails. It depends only on the
// supports, so a trace replays it at other primes.
pub(super) struct FinishPlan {
    len: usize,
    keep: Vec<usize>,
    leads: Vec<Monomial>,
    columns: Vec<Monomial>,
    pivots: Vec<(usize, Vec<u32>)>,
    tails: Vec<Vec<u32>>,
    // Positions in `keep`, by descending leading monomial.
    sorted: Vec<usize>,
}

impl FinishPlan {
    pub(super) fn new<F: Field>(basis: &[Polynomial<F>]) -> Result<Self, GroebnerError> {
        let first = basis.first().ok_or(GroebnerError::EmptyInput)?;
        let order = first.order.clone();
        let minimal = crate::groebner::minimal_mask(basis);
        let keep: Vec<usize> = (0..basis.len()).filter(|&i| minimal[i]).collect();
        let kept: Vec<&Polynomial<F>> = keep.iter().map(|&i| &basis[i]).collect();
        let masks: Vec<u64> = kept
            .iter()
            .map(|g| divisibility_mask(&g.terms[0].0))
            .collect();
        let reducer = |m: &Monomial| {
            let mask = divisibility_mask(m);
            (0..kept.len())
                .filter(|&i| masks[i] & !mask == 0 && kept[i].terms[0].0.divides(m))
                .min_by_key(|&i| kept[i].terms.len())
        };
        let mut table = MonomialTable::default();
        let mut queue = Vec::new();
        let fresh = |(id, fresh): (u32, bool), queue: &mut Vec<u32>| {
            if fresh {
                queue.push(id);
            }
            id
        };
        let tails: Vec<Vec<u32>> = kept
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
            let Some(mult) = m.quo(&kept[i].terms[0].0) else {
                continue;
            };
            let mut columns = vec![id];
            for (t, _) in &kept[i].terms[1..] {
                columns.push(fresh(table.intern_product(t, &mult), &mut queue));
            }
            pivots.push((i, columns));
        }
        let ids = table.columns(&order);
        let index = column_index(&ids);
        let reindex = |mut cols: Vec<u32>| {
            cols.iter_mut().for_each(|c| *c = index[*c as usize]);
            cols
        };
        let leads: Vec<Monomial> = kept.iter().map(|g| g.terms[0].0.clone()).collect();
        let mut sorted: Vec<usize> = (0..keep.len()).collect();
        sorted.sort_by(|&a, &b| order.compare(&leads[b], &leads[a]));
        Ok(Self {
            len: basis.len(),
            keep,
            leads,
            columns: ids.iter().map(|&id| table.monomials[id].clone()).collect(),
            pivots: pivots.into_iter().map(|(i, c)| (i, reindex(c))).collect(),
            tails: tails.into_iter().map(reindex).collect(),
            sorted,
        })
    }

    pub(super) fn execute<F: F4Field>(&self, basis: Vec<Polynomial<F>>) -> Vec<Polynomial<F>> {
        let (nvars, order) = (basis[0].nvars, basis[0].order.clone());
        let mut basis: Vec<_> = basis.into_iter().map(Some).collect();
        let kept: Vec<Polynomial<F>> = self
            .keep
            .iter()
            .filter_map(|&i| basis[i].take())
            .map(|g| {
                if g.lc().is_some_and(num_traits::One::is_one) {
                    g
                } else {
                    g.monic()
                }
            })
            .collect();
        let pivot_rows = reducers(
            kept.len(),
            |i| kept[i].terms.iter().map(|(_, c)| c.clone()).collect(),
            self.pivots.iter().cloned(),
        );
        let tail_rows: Vec<_> = kept
            .iter()
            .zip(&self.tails)
            .map(|(g, cols)| SparseRow {
                columns: cols.clone(),
                coefficients: g.terms[1..].iter().map(|(_, c)| c.clone()).collect(),
            })
            .collect();
        // The rows hold everything else, so only the leading terms are kept.
        let mut leads: Vec<_> = kept
            .into_iter()
            .map(|mut g| Some(g.terms.swap_remove(0)))
            .collect();
        let reduced = F::reduce_rows(&pivot_rows, &tail_rows, self.columns.len());
        drop((pivot_rows, tail_rows));
        self.sorted
            .iter()
            .filter_map(|&k| {
                let mut p = decode(&reduced[k], &self.columns, nvars, &order);
                p.terms.insert(0, leads[k].take()?);
                Some(p)
            })
            .collect()
    }

    // `execute` at the lanes' primes, for the active elements of a lane replay, or
    // `None` when their supports differ from the plan's or a residue splits them.
    pub(super) fn execute_lanes(
        &self,
        basis: &[&Polynomial<Zp>],
        coefficients: &[&Vec<Lanes>],
        primes: [u64; PRIMES],
    ) -> Option<Vec<Vec<Polynomial<Zp>>>> {
        if basis.len() != self.len {
            return None;
        }
        let same = self
            .keep
            .iter()
            .zip(&self.leads)
            .zip(&self.tails)
            .all(|((&i, lead), tail)| {
                let terms = &basis[i].terms;
                terms.len() == tail.len() + 1
                    && terms[0].0 == *lead
                    && terms[1..]
                        .iter()
                        .zip(tail)
                        .all(|((m, _), &c)| *m == self.columns[c as usize])
            });
        if !same {
            return None;
        }
        let ncols = self.columns.len();
        let kernel = LaneKernel::new(primes, ncols);
        let mut monic = Vec::with_capacity(self.keep.len());
        for &i in &self.keep {
            let mut row = SparseRow {
                columns: Vec::new(),
                coefficients: coefficients[i].clone(),
            };
            if row.coefficients[0].contains(&0) {
                return None;
            }
            kernel.normalize(&mut row);
            monic.push(row.coefficients);
        }
        let table = Table::new(
            self.pivots.iter().map(|(k, columns)| Row {
                columns,
                coefficients: &monic[*k],
            }),
            ncols,
        );
        let tails: Vec<_> = monic
            .iter()
            .zip(&self.tails)
            .map(|(cs, columns)| SparseRow {
                columns: columns.clone(),
                coefficients: cs[1..].to_vec(),
            })
            .collect();
        let reduced = par::map_blocks(
            &tails,
            kernel.block(tails.len()),
            || Kernel::<Lanes>::scratch(&kernel),
            |r, s| kernel.reduce(r, &table, s),
        );
        if kernel.failed.into_inner() {
            return None;
        }
        let (nvars, order) = (basis[0].nvars, basis[0].order.clone());
        let lane = |l: usize| {
            self.sorted
                .iter()
                .map(|&k| {
                    let lead = (self.leads[k].clone(), Zp::from_residue(1, primes[l]));
                    let tail = reduced[k]
                        .columns
                        .iter()
                        .zip(&reduced[k].coefficients)
                        .filter(|(_, v)| v[l] != 0)
                        .map(|(&c, v)| {
                            (
                                self.columns[c as usize].clone(),
                                Zp::from_residue(v[l].into(), primes[l]),
                            )
                        });
                    Polynomial {
                        terms: std::iter::once(lead).chain(tail).collect(),
                        nvars,
                        order: order.clone(),
                    }
                })
                .collect()
        };
        Some((0..PRIMES).map(lane).collect())
    }
}
