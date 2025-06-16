use std::cmp::{Ordering, Reverse};
use std::collections::{BinaryHeap, HashSet};

// --- Monomials, Terms, Polynomials (adapted for NUM_VARS=2) ---
const NUM_VARS: usize = 2; // For x, y

#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Monomial {
    exponents: [u32; NUM_VARS],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MonomialOrder {
    Lex,
    GrLex,
    GRevLex,
}

impl Monomial {
    pub fn new(exponents: [u32; NUM_VARS]) -> Self {
        Monomial { exponents }
    }
    fn one() -> Self {
        Monomial {
            exponents: [0; NUM_VARS],
        }
    }

    fn compare_with_order(&self, other: &Self, order: MonomialOrder) -> Ordering {
        match order {
            MonomialOrder::Lex => {
                for i in 0..NUM_VARS {
                    match self.exponents[i].cmp(&other.exponents[i]) {
                        Ordering::Less => return Ordering::Less,
                        Ordering::Greater => return Ordering::Greater,
                        Ordering::Equal => continue,
                    }
                }
                Ordering::Equal
            }
            MonomialOrder::GrLex => {
                for i in 0..NUM_VARS {
                    match self.exponents[i].cmp(&other.exponents[i]) {
                        Ordering::Less => return Ordering::Less,
                        Ordering::Greater => return Ordering::Greater,
                        Ordering::Equal => continue,
                    }
                }
                Ordering::Equal
            }
            MonomialOrder::GRevLex => {
                for i in 0..NUM_VARS {
                    match self.exponents[i].cmp(&other.exponents[i]) {
                        Ordering::Less => return Ordering::Greater,
                        Ordering::Greater => return Ordering::Less,
                        Ordering::Equal => continue,
                    }
                }
                Ordering::Equal
            }
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        let mut new_exponents = [0; NUM_VARS];
        for (i, (self_exp, other_exp)) in self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .enumerate()
            .take(NUM_VARS)
        {
            new_exponents[i] = self_exp + other_exp;
        }
        Monomial {
            exponents: new_exponents,
        }
    }

    fn is_multiple_of(&self, other: &Self) -> bool {
        // self is multiple of other
        for i in 0..NUM_VARS {
            if self.exponents[i] < other.exponents[i] {
                return false;
            }
        }
        true
    }

    fn divide(&self, other: &Self) -> Option<Self> {
        // self / other
        if !self.is_multiple_of(other) {
            return None;
        }
        let mut new_exponents = [0; NUM_VARS];
        for (i, (self_exp, other_exp)) in self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .enumerate()
            .take(NUM_VARS)
        {
            new_exponents[i] = self_exp - other_exp;
        }
        Some(Monomial {
            exponents: new_exponents,
        })
    }

    fn lcm(&self, other: &Self) -> Self {
        let mut new_exponents = [0; NUM_VARS];
        for (i, (self_exp, other_exp)) in self
            .exponents
            .iter()
            .zip(other.exponents.iter())
            .enumerate()
            .take(NUM_VARS)
        {
            new_exponents[i] = (*self_exp).max(*other_exp);
        }
        Monomial {
            exponents: new_exponents,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Term {
    coeff: f64,
    mono: Monomial,
}
impl Term {
    fn new(coeff: f64, mono: Monomial) -> Self {
        Term { coeff, mono }
    }
}

#[derive(Debug, Clone)] // Removed PartialEq for custom implementation
pub struct Polynomial {
    terms: Vec<Term>,
    order: MonomialOrder,
}

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        if self.order != other.order || self.terms.len() != other.terms.len() {
            return false;
        }
        for (t1, t2) in self.terms.iter().zip(other.terms.iter()) {
            if t1.mono != t2.mono || (t1.coeff - t2.coeff).abs() > 1e-9 {
                return false;
            }
        }
        true
    }
}

impl Polynomial {
    fn new(mut terms: Vec<Term>, order: MonomialOrder) -> Self {
        terms.retain(|t| t.coeff.abs() > 1e-9);
        terms.sort_by(|a, b| b.mono.compare_with_order(&a.mono, order));
        Polynomial { terms, order }
    }
    fn from_term(term: Term, order: MonomialOrder) -> Self {
        if term.coeff.abs() < 1e-9 {
            Polynomial::zero(order)
        } else {
            Polynomial {
                terms: vec![term],
                order,
            }
        }
    }
    fn zero(order: MonomialOrder) -> Self {
        Polynomial {
            terms: Vec::new(),
            order,
        }
    }
    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
    fn leading_term(&self) -> Option<&Term> {
        self.terms.first()
    }
    fn leading_monomial(&self) -> Option<&Monomial> {
        self.terms.first().map(|t| &t.mono)
    }
    fn leading_coefficient(&self) -> Option<f64> {
        self.terms.first().map(|t| t.coeff)
    }

    fn add(&self, other: &Self) -> Self {
        let mut result_terms = Vec::new();
        let mut i = 0;
        let mut j = 0;
        while i < self.terms.len() && j < other.terms.len() {
            match self.terms[i]
                .mono
                .compare_with_order(&other.terms[j].mono, self.order)
            {
                Ordering::Greater => {
                    result_terms.push(self.terms[i].clone());
                    i += 1;
                }
                Ordering::Less => {
                    result_terms.push(other.terms[j].clone());
                    j += 1;
                }
                Ordering::Equal => {
                    let nc = self.terms[i].coeff + other.terms[j].coeff;
                    if nc.abs() > 1e-9 {
                        result_terms.push(Term::new(nc, self.terms[i].mono.clone()));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        result_terms.extend_from_slice(&self.terms[i..]);
        result_terms.extend_from_slice(&other.terms[j..]);
        Polynomial::new(result_terms, self.order)
    }
    fn negate(&self) -> Self {
        Polynomial::new(
            self.terms
                .iter()
                .map(|t| Term::new(-t.coeff, t.mono.clone()))
                .collect(),
            self.order,
        )
    }
    fn subtract(&self, other: &Self) -> Self {
        self.add(&other.negate())
    }
    fn multiply_by_term(&self, term_m: &Term) -> Self {
        if term_m.coeff.abs() < 1e-9 || self.is_zero() {
            return Polynomial::zero(self.order);
        }
        Polynomial::new(
            self.terms
                .iter()
                .map(|t| Term::new(t.coeff * term_m.coeff, t.mono.multiply(&term_m.mono)))
                .collect(),
            self.order,
        )
    }
    fn s_polynomial(p1: &Polynomial, p2: &Polynomial) -> Option<Polynomial> {
        if p1.is_zero() || p2.is_zero() {
            return None;
        }
        let lt1 = p1.leading_term()?;
        let lt2 = p2.leading_term()?;
        let lcm_lm = lt1.mono.lcm(&lt2.mono);
        let t1_mono = lcm_lm.divide(&lt1.mono)?;
        let t1 = Term::new(1.0 / lt1.coeff, t1_mono);
        let t2_mono = lcm_lm.divide(&lt2.mono)?;
        let t2 = Term::new(1.0 / lt2.coeff, t2_mono);
        Some(p1.multiply_by_term(&t1).subtract(&p2.multiply_by_term(&t2)))
    }
    fn reduce_by_many(&self, divisors: &[Polynomial]) -> Polynomial {
        let mut p = self.clone();
        let mut r = Polynomial::zero(self.order);
        'reduction_loop: loop {
            if p.is_zero() {
                break;
            }
            let ltp_clone = p.leading_term().unwrap().clone();
            for div_poly in divisors {
                if div_poly.is_zero() {
                    continue;
                }
                let ltd = div_poly.leading_term().unwrap();
                if ltp_clone.mono.is_multiple_of(&ltd.mono) {
                    if let Some(quot_m) = ltp_clone.mono.divide(&ltd.mono) {
                        let quot_c = ltp_clone.coeff / ltd.coeff;
                        let quot_t = Term::new(quot_c, quot_m);
                        p = p.subtract(&div_poly.multiply_by_term(&quot_t));
                        continue 'reduction_loop;
                    }
                }
            }
            // If we get here, no reduction was possible
            r = r.add(&Polynomial::from_term(ltp_clone.clone(), self.order));
            p = p.subtract(&Polynomial::from_term(ltp_clone, self.order));
        }
        Polynomial::new(r.terms, r.order)
    }
    fn normalize(&self) -> Option<Polynomial> {
        if self.is_zero() {
            return None;
        }
        let lc = self.leading_coefficient().unwrap();
        if lc.abs() < 1e-9 {
            return None;
        }
        Some(self.multiply_by_term(&Term::new(1.0 / lc, Monomial::one())))
    }
}

// --- Buchberger Algorithm Data Structures ---
#[derive(Debug, Clone, Eq)]
struct CriticalPairItem {
    idx1: usize,
    idx2: usize,
    lcm_lm: Monomial,
    order: MonomialOrder, // For Ord/PartialOrd to use in comparison
}

impl PartialEq for CriticalPairItem {
    fn eq(&self, other: &Self) -> bool {
        // Sort order doesn't depend on idx order
        (self.idx1.min(self.idx2) == other.idx1.min(other.idx2))
            && (self.idx1.max(self.idx2) == other.idx1.max(other.idx2))
            && self.lcm_lm == other.lcm_lm // Technically, indices are enough if LCM is just for sorting
    }
}

impl Ord for CriticalPairItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.lcm_lm
            .compare_with_order(&other.lcm_lm, self.order)
            .then_with(|| {
                (self.idx1.min(self.idx2), self.idx1.max(self.idx2)) // Tie-breaking
                    .cmp(&(other.idx1.min(other.idx2), other.idx1.max(other.idx2)))
            })
    }
}

impl PartialOrd for CriticalPairItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// --- Buchberger's Algorithm ---
pub fn buchberger(initial_polys: Vec<Polynomial>, order: MonomialOrder) -> Vec<Polynomial> {
    if initial_polys.is_empty() {
        return Vec::new();
    }

    // 1. Initial reduction of input polynomials
    let mut f_processed = Vec::new();
    if !initial_polys.is_empty() {
        let mut current_f = initial_polys
            .iter()
            .filter_map(|p| p.normalize())
            .collect::<Vec<_>>();
        loop {
            let mut next_f = Vec::new();
            for i in 0..current_f.len() {
                let p = &current_f[i];
                // Reduce by polynomials from previous list `current_f` up to i-1
                // A simpler inter-reduction: reduce p by what's already in `next_f`.
                let temp_divs = current_f
                    .iter()
                    .enumerate()
                    .filter(|(k, _)| *k < i)
                    .map(|(_, poly)| poly.clone())
                    .collect::<Vec<_>>();
                let r = p.reduce_by_many(&temp_divs);

                if !r.is_zero() {
                    if let Some(rn) = r.normalize() {
                        next_f.push(rn);
                    }
                }
            }
            if next_f == current_f {
                // Stable state
                f_processed = next_f;
                break;
            }
            current_f = next_f;
            if current_f.is_empty() {
                break;
            } // All reduced to zero
        }
    }
    if f_processed.is_empty()
        && !initial_polys.is_empty()
        && initial_polys.iter().all(|p| p.is_zero())
    {
        return Vec::new();
    } else if f_processed.is_empty()
        && !initial_polys.iter().any(|p| p.is_zero())
        && !initial_polys.is_empty()
    {
        // This case can happen if initial reduction yields nothing but original polys were non-zero
        // Fallback: just use normalized initial polys if initial reduction to empty is suspicious
        f_processed = initial_polys.iter().filter_map(|p| p.normalize()).collect();
    }

    let mut g: Vec<Polynomial> = f_processed;
    if g.is_empty() {
        return Vec::new();
    }

    // 2. Initialize critical pairs (using min-heap for normal strategy)
    // `Reverse` makes BinaryHeap a min-heap
    let mut crit_pairs_heap: BinaryHeap<Reverse<CriticalPairItem>> = BinaryHeap::new();
    let mut existing_pairs_set: HashSet<(usize, usize)> = HashSet::new(); // To avoid duplicate pairs in heap

    for i in 0..g.len() {
        for j in (i + 1)..g.len() {
            if let (Some(lm_i), Some(lm_j)) = (g[i].leading_monomial(), g[j].leading_monomial()) {
                let key = (i.min(j), i.max(j));
                if existing_pairs_set.insert(key) {
                    crit_pairs_heap.push(Reverse(CriticalPairItem {
                        idx1: i,
                        idx2: j,
                        lcm_lm: lm_i.lcm(lm_j),
                        order,
                    }));
                }
            }
        }
    }

    // 3. Main loop
    let mut iter_count = 0;
    const MAX_ITERS: usize = 100; // Safety break

    while let Some(Reverse(pair_item)) = crit_pairs_heap.pop() {
        iter_count += 1;
        if iter_count > MAX_ITERS {
            println!("WARN: Max iterations reached in Buchberger");
            break;
        }

        let p1 = &g[pair_item.idx1];
        let p2 = &g[pair_item.idx2];

        if let Some(s_poly) = Polynomial::s_polynomial(p1, p2) {
            if s_poly.is_zero() {
                continue;
            }

            let h_reduced = s_poly.reduce_by_many(&g);

            if !h_reduced.is_zero() {
                if let Some(h_normalized) = h_reduced.normalize() {
                    let new_poly_idx = g.len();
                    for (i, poly) in g.iter().enumerate() {
                        if let (Some(lm_i), Some(lm_h)) =
                            (poly.leading_monomial(), h_normalized.leading_monomial())
                        {
                            let key = (i.min(new_poly_idx), i.max(new_poly_idx));
                            if existing_pairs_set.insert(key) {
                                crit_pairs_heap.push(Reverse(CriticalPairItem {
                                    idx1: i,
                                    idx2: new_poly_idx,
                                    lcm_lm: lm_i.lcm(lm_h),
                                    order,
                                }));
                            }
                        }
                    }
                    g.push(h_normalized);
                }
            }
        }
    }

    // 4. Reduce the resulting Groebner basis G
    // First, make it minimal (no LM divides another)
    let mut minimal_g = Vec::new();
    g.sort_by(|a, b| {
        a.leading_monomial()
            .unwrap_or(&Monomial::one())
            .compare_with_order(b.leading_monomial().unwrap_or(&Monomial::one()), order)
    });

    for p_i in g.iter() {
        if p_i.is_zero() {
            continue;
        }
        let lm_pi = p_i.leading_monomial().unwrap();
        let is_minimal_lm = !minimal_g.iter().any(|p_j: &Polynomial| {
            lm_pi.is_multiple_of(p_j.leading_monomial().unwrap())
                && lm_pi != p_j.leading_monomial().unwrap()
        });
        if is_minimal_lm {
            // Also ensure no existing poly in minimal_g has LM divisible by lm_pi
            minimal_g.retain(|p_j: &Polynomial| {
                !(p_j.leading_monomial().unwrap().is_multiple_of(lm_pi)
                    && p_j.leading_monomial().unwrap() != lm_pi)
            });
            minimal_g.push(p_i.clone());
            // Re-sort and unique LMs after potential removals and addition
            minimal_g.sort_by(|a, b| {
                a.leading_monomial()
                    .unwrap()
                    .compare_with_order(b.leading_monomial().unwrap(), order)
            });
            minimal_g.dedup_by_key(|p| p.leading_monomial().cloned());
        }
    }

    // Second, make it reduced (each poly is reduced w.r.t. others)
    let mut reduced_g = Vec::new();
    for (i, p_i) in minimal_g.iter().enumerate() {
        let mut other_divs = Vec::new();
        for (j, poly) in minimal_g.iter().enumerate() {
            if i == j {
                continue;
            }
            other_divs.push(poly.clone());
        }
        let h = p_i.reduce_by_many(&other_divs);
        if !h.is_zero() {
            // Normalization should already be handled, but good to ensure.
            if let Some(hn) = h.normalize() {
                reduced_g.push(hn);
            }
        }
    }

    // Final sort by leading monomial
    reduced_g.sort_by(|a, b| {
        a.leading_monomial()
            .unwrap_or(&Monomial::one())
            .compare_with_order(b.leading_monomial().unwrap_or(&Monomial::one()), order)
    });
    reduced_g
}

// --- Test ---
#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create polynomials for tests more easily
    // polys_spec: Vec of Vec of (coeff, [exp_x, exp_y])
    fn create_poly_from_spec(
        spec: Vec<(f64, [u32; NUM_VARS])>,
        order: MonomialOrder,
    ) -> Polynomial {
        let terms = spec
            .into_iter()
            .map(|(c, exps)| Term::new(c, Monomial::new(exps)))
            .collect();
        Polynomial::new(terms, order)
    }

    // Helper to compare polynomial lists (ignoring order of polys, but terms within polys are ordered)
    fn compare_poly_sets(set1: &[Polynomial], set2: &[Polynomial]) -> bool {
        if set1.len() != set2.len() {
            return false;
        }
        let mut set1_sorted = set1.to_vec();
        let mut set2_sorted = set2.to_vec();

        let order = set1_sorted.first().map_or(MonomialOrder::Lex, |p| p.order);

        set1_sorted.sort_by(|a, b| {
            a.leading_monomial()
                .unwrap()
                .compare_with_order(b.leading_monomial().unwrap(), order)
        });
        set2_sorted.sort_by(|a, b| {
            a.leading_monomial()
                .unwrap()
                .compare_with_order(b.leading_monomial().unwrap(), order)
        });

        set1_sorted == set2_sorted
    }

    #[test]
    fn test_buchberger_example_cox_little_oshea_chap2_sec7_ex2() {
        // Example: f1 = x^2 - y, f2 = xy - 1 with Lex order x > y
        // Expected Reduced Groebner Basis: {x - y^2, y^3 - 1} (or equivalent with LCs)
        // (after normalization: {x - y^2, y^3 - 1})

        let order = MonomialOrder::Lex;

        let f1 = create_poly_from_spec(vec![(1.0, [2, 0]), (-1.0, [0, 1])], order); // x^2 - y
        let f2 = create_poly_from_spec(vec![(1.0, [1, 1]), (-1.0, [0, 0])], order); // xy - 1

        let initial_polys = vec![f1, f2];
        let gb = buchberger(initial_polys, order);

        println!("Computed Groebner Basis:");
        for p in &gb {
            println!("{:?}", p.terms);
        }

        // Expected basis elements (normalized)
        let expected_g1 = create_poly_from_spec(vec![(1.0, [1, 0]), (-1.0, [0, 2])], order); // x - y^2
        let expected_g2 = create_poly_from_spec(vec![(1.0, [0, 3]), (-1.0, [0, 0])], order); // y^3 - 1
        let expected_basis = vec![expected_g1, expected_g2];

        // Note: Floating point comparisons can be tricky.
        // The Polynomial PartialEq handles term-wise comparison with tolerance via new().
        // compare_poly_sets sorts by LM then compares.
        assert!(
            compare_poly_sets(&gb, &expected_basis),
            "Groebner basis does not match expected"
        );
    }

    #[test]
    fn test_buchberger_already_groebner() {
        let order = MonomialOrder::Lex;
        let f1 = create_poly_from_spec(vec![(1.0, [1, 0])], order); // x
        let f2 = create_poly_from_spec(vec![(1.0, [0, 1])], order); // y
        // This is already a Groebner basis {x, y}
        let initial_polys = vec![f1.clone(), f2.clone()];
        let gb = buchberger(initial_polys, order);

        assert!(
            compare_poly_sets(&gb, &[f1, f2]),
            "Failed for already Groebner basis"
        );
    }

    #[test]
    fn test_buchberger_single_poly() {
        let order = MonomialOrder::Lex;
        let f1 = create_poly_from_spec(vec![(2.0, [1, 1]), (-4.0, [0, 0])], order); // 2xy - 4
        let initial_polys = vec![f1];
        let gb = buchberger(initial_polys, order);

        let expected_f1_normalized =
            create_poly_from_spec(vec![(1.0, [1, 1]), (-2.0, [0, 0])], order); // xy - 2
        assert!(
            compare_poly_sets(&gb, &[expected_f1_normalized]),
            "Failed for single polynomial"
        );
    }

    #[test]
    fn test_buchberger_zero_poly() {
        let order = MonomialOrder::Lex;
        let f1 = Polynomial::zero(order);
        let initial_polys = vec![f1];
        let gb = buchberger(initial_polys, order);
        assert!(gb.is_empty(), "Failed for zero polynomial input");

        let f2 = create_poly_from_spec(vec![(1.0, [1, 0])], order); // x
        let initial_polys_2 = vec![f2.clone(), Polynomial::zero(order)];
        let gb2 = buchberger(initial_polys_2, order);
        assert!(
            compare_poly_sets(&gb2, &[f2]),
            "Failed for input with a zero polynomial"
        );
    }
}
