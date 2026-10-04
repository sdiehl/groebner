//! The F4 basis state: Gebauer-Moller pair selection and symbolic preprocessing into
//! matrix plans over interned monomial columns.

use super::monomials::{
    Fresh, LANES, MonomialTable, lane_divides, lane_max, lane_sum, pack, unpack,
};
use super::{F4Field, Reducers, SparseRow, column_index, reducers};
use crate::monomial::{Monomial, MonomialOrder, divisibility_mask};
use crate::polynomial::Polynomial;
use crate::{Field, par};
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use std::borrow::Cow;

#[derive(Debug, Clone)]
pub(super) struct Pair {
    i: usize,
    j: usize,
    lcm: Monomial,
    degree: u32,
    mask: u64,
}

// A basis element: its coefficients in term order, the packed monomials of its terms when
// all fit, and otherwise, or when a trace records the support, the monomials themselves.
pub(super) struct Element<F> {
    pub(super) lm: Monomial,
    coefficients: Vec<F>,
    keys: Option<Box<[u128]>>,
    monomials: Vec<Monomial>,
}

impl<F: Field> Element<F> {
    pub(super) fn new(p: &Polynomial<F>, traced: bool) -> Self {
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

    pub(super) fn decode(
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

    pub(super) fn polynomial(
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
                    (m.clone(), c)
                })
                .collect(),
            _ => self.monomials.into_iter().zip(self.coefficients).collect(),
        };
        Polynomial {
            terms,
            nvars,
            order: order.clone(),
        }
    }
}

pub(super) struct State<F> {
    pub(super) basis: Vec<Element<F>>,
    pub(super) active: Vec<bool>,
    // Active indices by (term count, index), so the first divisor found is the shortest.
    by_size: Vec<usize>,
    pub(super) pairs: Vec<Pair>,
    masks: Vec<u64>,
    // Packed leading monomials, when they fit.
    keys: Vec<Option<u128>>,
    nvars: usize,
}

pub(super) struct Matrix<F> {
    pub(super) pivots: Reducers<F>,
    pub(super) rows: Vec<SparseRow<F>>,
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
    let free = par::map_big(&sorted, |(lp, mp, p, coprime)| {
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

pub(super) struct PlannedRow {
    basis: usize,
    // Support at planning time, kept only when recording a trace.
    terms: Vec<Monomial>,
    // Taken by the matrix unless the plan is traced.
    columns: Vec<u32>,
}
pub(super) struct MatrixPlan {
    pub(super) columns: Vec<Monomial>,
    pivots: Vec<PlannedRow>,
    pub(super) rows: Vec<PlannedRow>,
}
impl MatrixPlan {
    // Some rows of a traced plan, over a replayed basis whose coefficients, in term
    // order, are `coefficients(i)`.
    pub(super) fn execute_rows<F: Field, C>(
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
    pub(super) fn new(nvars: usize) -> Self {
        Self {
            basis: Vec::new(),
            active: Vec::new(),
            by_size: Vec::new(),
            pairs: Vec::new(),
            masks: Vec::new(),
            keys: Vec::new(),
            nvars,
        }
    }

    fn lm(&self, i: usize) -> &Monomial {
        &self.basis[i].lm
    }

    // The matrix of a plan made from this basis. Unless the plan is traced, each row's
    // columns are freed as they are copied, so only one row is held twice, and the copies
    // lie together, which elimination reads faster than the plan's scattered vectors.
    pub(super) fn matrix(&self, plan: &mut MatrixPlan, traced: bool) -> Matrix<F> {
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
    pub(super) fn update(&mut self, mut h: Element<F>) {
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
                let cands = par::map_big(&others, |&g| {
                    let k = self.keys[g].unwrap_or_default();
                    let lcm = lane_max(k, h);
                    (lcm, lane_sum(lcm), lcm == k + h, mask(g))
                });
                survivors(&cands, |&a, &b| lane_divides(a, b))
            }
            _ => {
                let cands = par::map_big(&others, |&g| {
                    let lcm = self.lm(g).lcm(lm_h);
                    let degree = lcm.degree();
                    (lcm, degree, self.lm(g).is_coprime(lm_h), mask(g))
                });
                survivors(&cands, Monomial::divides)
            }
        };
        kept.into_iter().map(|k| self.pair(others[k], t)).collect()
    }

    pub(super) fn select(&mut self) -> Vec<Pair> {
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

    pub(super) fn symbolic_preprocessing(
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
            let found = par::map_weighted(
                &level,
                |(i, _, _)| self.basis[*i].len(),
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
            let found = par::map_rows(
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
        par::for_rows_mut(&mut pivots, reindex);
        par::for_rows_mut(&mut rows, reindex);
        pivots.sort_by_key(|r| r.columns[0]);
        MatrixPlan {
            columns,
            pivots,
            rows,
        }
    }
}
