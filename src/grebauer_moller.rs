//! The selection  strategy of Gebauer and Möller
//!
//! From the paper "On an installation of Buchberger's algorithm."

use crate::field::Field;
use crate::groebner::CriticalPair;
use crate::polynomial::Polynomial;
use std::collections::HashSet;

/// Filter pairs using the Gebauer–Möller criteria (B, M, F)
pub fn filter_gm_pairs<F: Field>(
    basis: &[Polynomial<F>],
    pairs: Vec<CriticalPair>,
) -> Vec<CriticalPair> {
    // B-criterion: Remove pairs (i, j) if lcm(i, j) is divisible by lcm(i, k) for some k != j
    let mut to_remove = vec![false; pairs.len()];
    for idx in 0..pairs.len() {
        let pair = &pairs[idx];
        for k in 0..basis.len() {
            if k == pair.j || k == pair.i {
                continue;
            }
            if let (Some(lm_i), Some(lm_k)) = (
                basis[pair.i].leading_monomial(),
                basis[k].leading_monomial(),
            ) {
                let lcm_ik = lm_i.lcm(lm_k);
                if lcm_ik != pair.lcm && pair.lcm.divides(&lcm_ik) {
                    to_remove[idx] = true;
                    break;
                }
            }
        }
    }
    // M-criterion: Remove pairs (i, j) if there exists another pair (i, k) with lcm(i, k) dividing lcm(i, j)
    for idx in 0..pairs.len() {
        if to_remove[idx] {
            continue;
        }
        let pair = &pairs[idx];
        for kidx in 0..pairs.len() {
            if idx == kidx || to_remove[kidx] {
                continue;
            }
            let other = &pairs[kidx];
            if pair.i == other.i && other.lcm.divides(&pair.lcm) && other.lcm != pair.lcm {
                to_remove[idx] = true;
                break;
            }
        }
    }
    // F-criterion: Remove duplicate pairs with the same lcm
    let mut seen = HashSet::new();
    for idx in 0..pairs.len() {
        if to_remove[idx] {
            continue;
        }
        let lcm = &pairs[idx].lcm;
        if !seen.insert(lcm.clone()) {
            to_remove[idx] = true;
        }
    }
    // Retain only pairs not marked for removal
    pairs
        .into_iter()
        .enumerate()
        .filter(|(idx, _)| !to_remove[*idx])
        .map(|(_, p)| p)
        .collect()
}
