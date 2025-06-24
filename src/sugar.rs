// Sugar strategy for Groebner basis computation
// The "sugar" of a polynomial is typically the total degree of the original input polynomial that led to the current one.
// from the paper “One sugar cube, please” or selection strategies in the Buchberger algorithm"

use crate::polynomial::Polynomial;

/// Struct to associate a polynomial with its sugar value.
#[derive(Debug, Clone)]
pub struct SugaredPolynomial<P> {
    pub poly: Polynomial<P>,
    pub sugar: usize,
}

impl<P: Clone> SugaredPolynomial<P> {
    /// Create a new SugaredPolynomial, computing the sugar as the total degree.
    pub fn new(poly: Polynomial<P>) -> Self {
        let sugar = total_degree(&poly);
        Self { poly, sugar }
    }
}

// Only implement Eq/PartialEq/Ord/PartialOrd if Polynomial<P>: PartialEq + Eq
impl<P: PartialEq + Eq> PartialEq for SugaredPolynomial<P> {
    fn eq(&self, other: &Self) -> bool {
        self.poly == other.poly && self.sugar == other.sugar
    }
}
impl<P: PartialEq + Eq> Eq for SugaredPolynomial<P> {}
impl<P: PartialEq + Eq> Ord for SugaredPolynomial<P> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.sugar.cmp(&other.sugar)
    }
}
impl<P: PartialEq + Eq> PartialOrd for SugaredPolynomial<P> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Compute the total degree of a polynomial (max of sum of exponents in any term's monomial)
fn total_degree<P>(poly: &Polynomial<P>) -> usize {
    poly.terms
        .iter()
        .map(|t| {
            t.monomial
                .exponents
                .iter()
                .map(|&e| e as usize)
                .sum::<usize>()
        })
        .max()
        .unwrap_or(0)
}

// Utility function to select the next S-polynomial to reduce based on sugar value.
pub fn select_next_by_sugar<P: Clone>(
    queue: &mut Vec<SugaredPolynomial<P>>,
) -> Option<SugaredPolynomial<P>> {
    if queue.is_empty() {
        None
    } else {
        // Find the polynomial with the lowest sugar value
        let min_idx = queue
            .iter()
            .enumerate()
            .min_by_key(|(_, sp)| sp.sugar)
            .map(|(i, _)| i)?;
        Some(queue.remove(min_idx))
    }
}
