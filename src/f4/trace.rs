//! Traces of F4 matrix plans learned at one prime and replayed at others, one prime
//! at a time or four side by side in lanes of one sweep.

use super::eliminate::{Kernel, Pivots, Table, echelon, fits, span};
use super::finish::{FinishPlan, finish};
use super::symbolic::MatrixPlan;
use super::{F4Field, SparseRow, compute, decode};
use crate::finite_field::Zp;
use crate::groebner::{GroebnerError, prepare_input};
use crate::monomial::Monomial;
use crate::polynomial::Polynomial;
use crate::{ModularField, par};
use polycore::modp::{inv as inv_mod, mul as mul_mod};
use std::sync::atomic::{AtomicBool, Ordering};

/// Coefficient-independent F4 plans, learned at a prime and reusable at other primes.
#[derive(Default)]
pub(crate) struct F4Trace {
    pub(super) input_leads: Vec<Monomial>,
    pub(super) rounds: Vec<TraceRound>,
    pub(super) active: Vec<bool>,
    pub(super) finish: Option<FinishPlan>,
}
pub(super) struct TraceRound {
    pub(super) plan: MatrixPlan,
    pub(super) leads: Vec<Monomial>,
    pub(super) live: Vec<usize>,
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
        let active: Vec<usize> = (0..basis.len()).filter(|&i| self.active[i]).collect();
        let finished = self.finish.as_ref().and_then(|plan| {
            let elements: Vec<_> = active.iter().map(|&i| &basis[i]).collect();
            let lanes: Vec<_> = active.iter().map(|&i| &coefficients[i]).collect();
            plan.execute_lanes(&elements, &lanes, primes)
        });
        if finished.is_some() {
            return finished;
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
                                (m.clone(), Zp::from_residue(v[l].into(), primes[l]))
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
pub(super) type Lanes = [u32; PRIMES];

// Dense reduction of lane residues. A leading coefficient that vanishes at only some
// primes splits the images' structure, which is recorded in `failed`.
pub(super) struct LaneKernel {
    p: [u64; PRIMES],
    ncols: usize,
    deferred: bool,
    pub(super) failed: AtomicBool,
}
impl LaneKernel {
    pub(super) fn new(p: [u64; PRIMES], ncols: usize) -> Self {
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
    // A replay round has few rows, each costly, so a block per worker would leave most
    // of the pool idle.
    fn block(&self, n: usize) -> usize {
        n.div_ceil(4 * par::threads()).clamp(1, 32)
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
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: AVX2 support was just detected.
        #[allow(unsafe_code)]
        return unsafe { accumulate_lanes_avx2(buf, tail, c) };
    }
    for (&pc, pv) in tail {
        let cell = &mut buf[pc as usize];
        for l in 0..PRIMES {
            cell[l] += u64::from(c[l]) * u64::from(pv[l]);
        }
    }
}

// One widening multiply and add covers the four lanes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[allow(unsafe_code)]
unsafe fn accumulate_lanes_avx2<'a>(
    buf: &mut [[u64; PRIMES]],
    tail: impl Iterator<Item = (&'a u32, &'a Lanes)>,
    c: Lanes,
) {
    use std::arch::x86_64::{
        __m128i, __m256i, _mm_loadu_si128, _mm256_add_epi64, _mm256_cvtepu32_epi64,
        _mm256_loadu_si256, _mm256_mul_epu32, _mm256_storeu_si256,
    };
    const { assert!(PRIMES == 4) };
    // SAFETY: each pointer addresses a whole array of four `u32` or four `u64` lanes.
    unsafe {
        let m = _mm256_cvtepu32_epi64(_mm_loadu_si128(c.as_ptr().cast::<__m128i>()));
        for (&pc, pv) in tail {
            let cell = buf[pc as usize].as_mut_ptr().cast::<__m256i>();
            let v = _mm256_cvtepu32_epi64(_mm_loadu_si128(pv.as_ptr().cast::<__m128i>()));
            let acc = _mm256_add_epi64(_mm256_loadu_si256(cell), _mm256_mul_epu32(v, m));
            _mm256_storeu_si256(cell, acc);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::f4::groebner_basis_f4_direct;
    use crate::{MonomialOrder, PolynomialRing};

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
}
