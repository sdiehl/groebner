//! Groebner-specific operations on polycore polynomials.
use crate::monomial::divisibility_mask;
use crate::{Field, Monomial};
pub use polycore::{Poly as Polynomial, Term};

/// Division against a basis, reducing by the shortest applicable divisor at each step.
pub trait PolynomialExt<F: Field> {
    /// Full normal form of `self` modulo `basis`: no remaining term is divisible by any leading monomial.
    fn normal_form(&self, basis: &[Polynomial<F>]) -> Polynomial<F>;
    /// Multivariate division: quotients `q` and remainder `r` with `self = sum q[i] * divisors[i] + r`.
    fn divide_with_remainder(
        &self,
        divisors: &[Polynomial<F>],
    ) -> (Vec<Polynomial<F>>, Polynomial<F>);
}

impl<F: Field> PolynomialExt<F> for Polynomial<F> {
    fn normal_form(&self, basis: &[Polynomial<F>]) -> Polynomial<F> {
        divide_by(self, basis.iter().enumerate(), |_, _, _| {})
    }

    fn divide_with_remainder(
        &self,
        divisors: &[Polynomial<F>],
    ) -> (Vec<Polynomial<F>>, Polynomial<F>) {
        let mut quotients = vec![Vec::new(); divisors.len()];
        let remainder = divide_by(self, divisors.iter().enumerate(), |i, c, m| {
            quotients[i].push((m.clone(), c.clone()));
        });
        let quotients = quotients
            .into_iter()
            .map(|terms| Polynomial {
                terms,
                nvars: self.nvars,
                order: self.order.clone(),
            })
            .collect();
        (quotients, remainder)
    }
}

/// Division loop reporting each quotient term `c * m` against divisor `i`, in descending order.
/// Unlike `Poly::divide`, which takes the first divisor that applies, this takes the shortest.
pub(crate) fn divide_by<'a, F: Field + 'a>(
    polynomial: &Polynomial<F>,
    divisors: impl Iterator<Item = (usize, &'a Polynomial<F>)>,
    mut record: impl FnMut(usize, &F, &Monomial),
) -> Polynomial<F> {
    let leads: Vec<_> = divisors
        .filter_map(|(i, g)| g.lt().map(|lt| (i, lt, g, divisibility_mask(&lt.0))))
        .collect();
    let mut work = polynomial.clone();
    let mut head = 0;
    let mut rest = Vec::new();
    while head < work.terms.len() {
        let lt = &work.terms[head];
        let mask = divisibility_mask(&lt.0);
        let reducer = leads
            .iter()
            .filter(|(_, glt, _, divisor_mask)| divisor_mask & !mask == 0 && glt.0.divides(&lt.0))
            .min_by_key(|(_, _, g, _)| g.terms.len())
            .and_then(|(i, glt, g, _)| lt.0.quo(&glt.0).map(|m| (*i, glt, g, m)));
        if let Some((i, glt, g, m)) = reducer {
            let c = lt.1.clone() / glt.1.clone();
            record(i, &c, &m);
            work.terms.drain(..head);
            work = work.sub_mul(&c, &m, g);
            head = 0;
        } else {
            rest.push(lt.clone());
            head += 1;
        }
    }
    Polynomial {
        terms: rest,
        nvars: polynomial.nvars,
        order: polynomial.order.clone(),
    }
}
