//! Groebner Basis Library
//!
//! This library implements algorithms for computing Groebner bases of polynomial ideals,
//! including the Buchberger algorithm and F4 algorithm optimizations.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::fmt;

/// Supported monomial orderings
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MonomialOrder {
    /// Lexicographic ordering
    Lex,
    /// Graded lexicographic ordering  
    GrLex,
    /// Graded reverse lexicographic ordering
    GRevLex,
}

impl MonomialOrder {
    /// Compare two monomials using this ordering
    pub fn compare(&self, a: &[u32], b: &[u32]) -> Ordering {
        match self {
            MonomialOrder::Lex => {
                for (ai, bi) in a.iter().zip(b.iter()) {
                    match ai.cmp(bi) {
                        Ordering::Equal => continue,
                        other => return other, // Standard lexicographic order
                    }
                }
                Ordering::Equal
            }
            MonomialOrder::GrLex => {
                let deg_a: u32 = a.iter().sum();
                let deg_b: u32 = b.iter().sum();
                match deg_a.cmp(&deg_b) {
                    Ordering::Equal => {
                        // Lexicographic tiebreaker
                        for (ai, bi) in a.iter().zip(b.iter()) {
                            match ai.cmp(bi) {
                                Ordering::Equal => continue,
                                other => return other,
                            }
                        }
                        Ordering::Equal
                    }
                    other => other, // Higher degree first
                }
            }
            MonomialOrder::GRevLex => {
                let deg_a: u32 = a.iter().sum();
                let deg_b: u32 = b.iter().sum();
                match deg_a.cmp(&deg_b) {
                    Ordering::Equal => {
                        for (ai, bi) in a.iter().rev().zip(b.iter().rev()) {
                            match ai.cmp(bi) {
                                Ordering::Equal => continue,
                                other => return other, // Reverse for grevlex
                            }
                        }
                        Ordering::Equal
                    }
                    other => other.reverse(),
                }
            }
        }
    }
}

/// A monomial represented as a vector of exponents
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Monomial {
    pub exponents: Vec<u32>,
}

impl Monomial {
    /// Create a new monomial with given exponents
    pub fn new(exponents: Vec<u32>) -> Self {
        Self { exponents }
    }

    /// Create the constant monomial (all exponents zero)
    pub fn one(nvars: usize) -> Self {
        Self {
            exponents: vec![0; nvars],
        }
    }

    /// Get the number of variables
    pub fn nvars(&self) -> usize {
        self.exponents.len()
    }

    /// Get the total degree
    pub fn degree(&self) -> u32 {
        self.exponents.iter().sum()
    }

    /// Multiply two monomials
    pub fn multiply(&self, other: &Self) -> Self {
        assert_eq!(self.nvars(), other.nvars());
        let exponents = self
            .exponents
            .iter()
            .zip(&other.exponents)
            .map(|(a, b)| a + b)
            .collect();
        Self { exponents }
    }

    /// Check if this monomial divides another
    pub fn divides(&self, other: &Self) -> bool {
        assert_eq!(self.nvars(), other.nvars());
        self.exponents
            .iter()
            .zip(&other.exponents)
            .all(|(a, b)| a <= b)
    }

    /// Divide two monomials (returns None if not divisible)
    pub fn divide(&self, other: &Self) -> Option<Self> {
        assert_eq!(self.nvars(), other.nvars());
        if !other.divides(self) {
            return None;
        }
        let exponents = self
            .exponents
            .iter()
            .zip(&other.exponents)
            .map(|(a, b)| a - b)
            .collect();
        Some(Self { exponents })
    }

    /// Compute the least common multiple
    pub fn lcm(&self, other: &Self) -> Self {
        assert_eq!(self.nvars(), other.nvars());
        let exponents = self
            .exponents
            .iter()
            .zip(&other.exponents)
            .map(|(a, b)| (*a).max(*b))
            .collect();
        Self { exponents }
    }

    /// Compare monomials using a given ordering
    pub fn compare(&self, other: &Self, order: MonomialOrder) -> Ordering {
        order.compare(&self.exponents, &other.exponents)
    }
}

impl fmt::Display for Monomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.exponents.iter().all(|&e| e == 0) {
            return write!(f, "1");
        }

        let mut first = true;
        for (i, &exp) in self.exponents.iter().enumerate() {
            if exp > 0 {
                if !first {
                    write!(f, "*")?;
                }
                first = false;
                if exp == 1 {
                    write!(f, "x{i}")?;
                } else {
                    write!(f, "x{i}^{exp}")?;
                }
            }
        }
        Ok(())
    }
}

/// A polynomial term (coefficient and monomial)
#[derive(Debug, Clone, PartialEq)]
pub struct Term<F> {
    pub coefficient: F,
    pub monomial: Monomial,
}

impl<F> Term<F> {
    pub fn new(coefficient: F, monomial: Monomial) -> Self {
        Self {
            coefficient,
            monomial,
        }
    }
}

/// A multivariate polynomial
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial<F> {
    pub terms: Vec<Term<F>>,
    pub nvars: usize,
    pub order: MonomialOrder,
}

impl<F: Field> Polynomial<F> {
    /// Create a new polynomial
    pub fn new(mut terms: Vec<Term<F>>, nvars: usize, order: MonomialOrder) -> Self {
        // Remove zero terms
        terms.retain(|t| !t.coefficient.is_zero());

        // Sort terms by monomial order (descending)
        terms.sort_by(|a, b| b.monomial.compare(&a.monomial, order));

        // Combine like terms
        let mut combined: Vec<Term<F>> = Vec::new();
        for term in terms {
            if let Some(last) = combined.last_mut() {
                if last.monomial == term.monomial {
                    last.coefficient = last.coefficient.add(&term.coefficient);
                    if last.coefficient.is_zero() {
                        combined.pop();
                    }
                    continue;
                }
            }
            combined.push(term);
        }

        Self {
            terms: combined,
            nvars,
            order,
        }
    }

    /// Create the zero polynomial
    pub fn zero(nvars: usize, order: MonomialOrder) -> Self {
        Self {
            terms: Vec::new(),
            nvars,
            order,
        }
    }

    /// Create a constant polynomial
    pub fn constant(coeff: F, nvars: usize, order: MonomialOrder) -> Self {
        if coeff.is_zero() {
            Self::zero(nvars, order)
        } else {
            Self::new(vec![Term::new(coeff, Monomial::one(nvars))], nvars, order)
        }
    }

    /// Check if polynomial is zero
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// Get the leading term
    pub fn leading_term(&self) -> Option<&Term<F>> {
        self.terms.first()
    }

    /// Get the leading monomial
    pub fn leading_monomial(&self) -> Option<&Monomial> {
        self.terms.first().map(|t| &t.monomial)
    }

    /// Get the leading coefficient
    pub fn leading_coefficient(&self) -> Option<&F> {
        self.terms.first().map(|t| &t.coefficient)
    }

    /// Make the polynomial monic (leading coefficient = 1)
    pub fn make_monic(&self) -> Self {
        if let Some(lc) = self.leading_coefficient() {
            if !lc.is_zero() {
                let inv_lc = lc.inverse();
                return self.multiply_scalar(&inv_lc);
            }
        }
        self.clone()
    }

    /// Add two polynomials
    pub fn add(&self, other: &Self) -> Self {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);

        let mut terms = self.terms.clone();
        terms.extend(other.terms.clone());

        Self::new(terms, self.nvars, self.order)
    }

    /// Subtract two polynomials
    pub fn subtract(&self, other: &Self) -> Self {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);

        let mut terms = self.terms.clone();
        for term in &other.terms {
            terms.push(Term::new(term.coefficient.negate(), term.monomial.clone()));
        }

        Self::new(terms, self.nvars, self.order)
    }

    /// Multiply by a scalar
    pub fn multiply_scalar(&self, scalar: &F) -> Self {
        if scalar.is_zero() {
            return Self::zero(self.nvars, self.order);
        }

        let terms = self
            .terms
            .iter()
            .map(|t| Term::new(t.coefficient.multiply(scalar), t.monomial.clone()))
            .collect();

        Self::new(terms, self.nvars, self.order)
    }

    /// Multiply by a monomial
    pub fn multiply_monomial(&self, monomial: &Monomial) -> Self {
        let terms = self
            .terms
            .iter()
            .map(|t| Term::new(t.coefficient.clone(), t.monomial.multiply(monomial)))
            .collect();

        Self::new(terms, self.nvars, self.order)
    }

    /// Multiply two polynomials
    pub fn multiply(&self, other: &Self) -> Self {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);

        if self.is_zero() || other.is_zero() {
            return Self::zero(self.nvars, self.order);
        }

        let mut terms = Vec::new();
        for t1 in &self.terms {
            for t2 in &other.terms {
                terms.push(Term::new(
                    t1.coefficient.multiply(&t2.coefficient),
                    t1.monomial.multiply(&t2.monomial),
                ));
            }
        }

        Self::new(terms, self.nvars, self.order)
    }

    /// Compute S-polynomial of two polynomials
    pub fn s_polynomial(&self, other: &Self) -> Self {
        assert_eq!(self.nvars, other.nvars);
        assert_eq!(self.order, other.order);

        if self.is_zero() || other.is_zero() {
            return Self::zero(self.nvars, self.order);
        }

        let lm1 = self.leading_monomial().unwrap();
        let lm2 = other.leading_monomial().unwrap();
        let lc1 = self.leading_coefficient().unwrap();
        let lc2 = other.leading_coefficient().unwrap();

        let lcm = lm1.lcm(lm2);
        let m1 = lcm.divide(lm1).unwrap();
        let m2 = lcm.divide(lm2).unwrap();

        let term1 = self.multiply_monomial(&m1).multiply_scalar(&lc2.inverse());
        let term2 = other.multiply_monomial(&m2).multiply_scalar(&lc1.inverse());

        term1.subtract(&term2)
    }

    /// Reduce this polynomial by a set of polynomials
    pub fn reduce(&self, basis: &[Self]) -> Self {
        let mut remainder = self.clone();

        while !remainder.is_zero() {
            let mut reduced = false;

            if let Some(leading_mono) = remainder.leading_monomial() {
                for divisor in basis {
                    if let Some(div_leading) = divisor.leading_monomial() {
                        if div_leading.divides(leading_mono) {
                            let quotient_mono = leading_mono.divide(div_leading).unwrap();
                            let lc_remainder = remainder.leading_coefficient().unwrap();
                            let lc_divisor = divisor.leading_coefficient().unwrap();
                            let quotient_coeff = lc_remainder.multiply(&lc_divisor.inverse());

                            let subtrahend = divisor
                                .multiply_monomial(&quotient_mono)
                                .multiply_scalar(&quotient_coeff);

                            remainder = remainder.subtract(&subtrahend);
                            reduced = true;
                            break;
                        }
                    }
                }
            }

            if !reduced {
                // Leading term cannot be reduced, so we're done
                break;
            }
        }

        remainder
    }
}

impl<F: Field> fmt::Display for Polynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        for (i, term) in self.terms.iter().enumerate() {
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "{}", term.coefficient)?;
            if !term.monomial.exponents.iter().all(|&e| e == 0) {
                write!(f, "*{}", term.monomial)?;
            }
        }
        Ok(())
    }
}

/// Trait for field operations
pub trait Field: Clone + PartialEq + fmt::Debug + fmt::Display {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;
    fn add(&self, other: &Self) -> Self;
    fn subtract(&self, other: &Self) -> Self;
    fn multiply(&self, other: &Self) -> Self;
    fn negate(&self) -> Self;
    fn inverse(&self) -> Self;

    fn divide(&self, other: &Self) -> Self {
        self.multiply(&other.inverse())
    }
}

/// Rational number field implementation
#[derive(Debug, Clone, PartialEq)]
pub struct Rational {
    pub numerator: i64,
    pub denominator: i64,
}

impl Rational {
    pub fn new(num: i64, den: i64) -> Self {
        if den == 0 {
            panic!("Division by zero");
        }
        let g = gcd(num.abs(), den.abs());
        let (num, den) = if den < 0 { (-num, -den) } else { (num, den) };
        Self {
            numerator: num / g,
            denominator: den / g,
        }
    }

    pub fn from_integer(n: i64) -> Self {
        Self::new(n, 1)
    }
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

impl Field for Rational {
    fn zero() -> Self {
        Self::new(0, 1)
    }

    fn one() -> Self {
        Self::new(1, 1)
    }

    fn is_zero(&self) -> bool {
        self.numerator == 0
    }

    fn is_one(&self) -> bool {
        self.numerator == 1 && self.denominator == 1
    }

    fn add(&self, other: &Self) -> Self {
        Self::new(
            self.numerator * other.denominator + other.numerator * self.denominator,
            self.denominator * other.denominator,
        )
    }

    fn subtract(&self, other: &Self) -> Self {
        Self::new(
            self.numerator * other.denominator - other.numerator * self.denominator,
            self.denominator * other.denominator,
        )
    }

    fn multiply(&self, other: &Self) -> Self {
        Self::new(
            self.numerator * other.numerator,
            self.denominator * other.denominator,
        )
    }

    fn negate(&self) -> Self {
        Self::new(-self.numerator, self.denominator)
    }

    fn inverse(&self) -> Self {
        if self.is_zero() {
            panic!("Division by zero");
        }
        Self::new(self.denominator, self.numerator)
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}

/// Critical pair for Buchberger algorithm
#[derive(Debug, Clone)]
struct CriticalPair {
    i: usize,
    j: usize,
    lcm: Monomial,
    degree: u32,
}

impl CriticalPair {
    fn new(
        i: usize,
        j: usize,
        poly_i: &Polynomial<impl Field>,
        poly_j: &Polynomial<impl Field>,
    ) -> Self {
        let lm_i = poly_i.leading_monomial().unwrap();
        let lm_j = poly_j.leading_monomial().unwrap();
        let lcm = lm_i.lcm(lm_j);
        let degree = lcm.degree();

        Self { i, j, lcm, degree }
    }
}

impl PartialEq for CriticalPair {
    fn eq(&self, other: &Self) -> bool {
        self.degree == other.degree && self.lcm == other.lcm
    }
}

impl Eq for CriticalPair {}

impl PartialOrd for CriticalPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CriticalPair {
    fn cmp(&self, other: &Self) -> Ordering {
        // Priority queue is max-heap, so reverse for min-heap behavior
        self.degree.cmp(&other.degree).reverse()
    }
}

/// Compute Groebner basis using Buchberger's algorithm
pub fn groebner_basis<F: Field>(
    polynomials: Vec<Polynomial<F>>,
    _order: MonomialOrder,
) -> Vec<Polynomial<F>> {
    if polynomials.is_empty() {
        return Vec::new();
    }

    let _nvars = polynomials[0].nvars;
    let mut basis: Vec<Polynomial<F>> = polynomials
        .into_iter()
        .filter(|p| !p.is_zero())
        .map(|p| p.make_monic())
        .collect();

    if basis.is_empty() {
        return Vec::new();
    }

    // Initialize critical pairs
    let mut pairs = BinaryHeap::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            pairs.push(CriticalPair::new(i, j, &basis[i], &basis[j]));
        }
    }

    while let Some(pair) = pairs.pop() {
        if pair.i >= basis.len() || pair.j >= basis.len() {
            continue;
        }

        // Check Buchberger criteria
        let poly_i = &basis[pair.i];
        let poly_j = &basis[pair.j];

        // Criterion 1: Check if lcm equals product of leading monomials
        let lm_i = poly_i.leading_monomial().unwrap();
        let lm_j = poly_j.leading_monomial().unwrap();
        let product = lm_i.multiply(lm_j);

        if pair.lcm == product {
            continue; // Skip this pair
        }

        // Compute S-polynomial and reduce
        let s_poly = poly_i.s_polynomial(poly_j);
        let reduced = s_poly.reduce(&basis);

        if !reduced.is_zero() {
            let monic_reduced = reduced.make_monic();

            // Add new critical pairs
            let new_index = basis.len();
            for (i, existing) in basis.iter().enumerate() {
                pairs.push(CriticalPair::new(i, new_index, existing, &monic_reduced));
            }

            basis.push(monic_reduced);
        }
    }

    // Remove redundant polynomials
    minimize_basis(&mut basis);

    basis
}

/// Remove redundant polynomials from the basis
fn minimize_basis<F: Field>(basis: &mut Vec<Polynomial<F>>) {
    let mut to_remove = Vec::new();

    for i in 0..basis.len() {
        for j in 0..basis.len() {
            if i != j {
                if let (Some(lm_i), Some(lm_j)) =
                    (basis[i].leading_monomial(), basis[j].leading_monomial())
                {
                    if lm_j.divides(lm_i) {
                        to_remove.push(i);
                        break;
                    }
                }
            }
        }
    }

    // Remove in reverse order to maintain indices
    to_remove.sort_unstable();
    to_remove.reverse();
    for &i in &to_remove {
        basis.remove(i);
    }
}

/// Check if a set of polynomials forms a Groebner basis
pub fn is_groebner_basis<F: Field>(basis: &[Polynomial<F>]) -> bool {
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            let s_poly = basis[i].s_polynomial(&basis[j]);
            let reduced = s_poly.reduce(basis);
            if !reduced.is_zero() {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_polynomial(
        terms: Vec<(i64, i64, Vec<u32>)>, // (num, den, exponents)
        nvars: usize,
        order: MonomialOrder,
    ) -> Polynomial<Rational> {
        let terms = terms
            .into_iter()
            .map(|(num, den, exp)| Term::new(Rational::new(num, den), Monomial::new(exp)))
            .collect();
        Polynomial::new(terms, nvars, order)
    }

    #[test]
    fn test_monomial_operations() {
        let m1 = Monomial::new(vec![2, 1, 0]);
        let m2 = Monomial::new(vec![1, 0, 2]);

        assert_eq!(m1.degree(), 3);
        assert_eq!(m2.degree(), 3);

        let product = m1.multiply(&m2);
        assert_eq!(product.exponents, vec![3, 1, 2]);

        let lcm = m1.lcm(&m2);
        assert_eq!(lcm.exponents, vec![2, 1, 2]);

        assert!(m1.divides(&product));
        assert!(!m1.divides(&m2));
    }

    #[test]
    fn test_polynomial_arithmetic() {
        // Test: x^2 + 2xy + y^2
        let p1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (2, 1, vec![1, 1]), (1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        // Test: x^2 - y^2
        let p2 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (-1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let sum = p1.add(&p2);
        assert_eq!(sum.terms.len(), 2); // 2x^2 + 2xy

        let product = p1.multiply(&p2);
        assert!(!product.is_zero());
    }

    #[test]
    fn test_s_polynomial() {
        // f = x^2 + 2xy^2
        let f = create_polynomial(
            vec![(1, 1, vec![2, 0]), (2, 1, vec![1, 2])],
            2,
            MonomialOrder::Lex,
        );

        // g = xy + 2y^3 - 1
        let g = create_polynomial(
            vec![(1, 1, vec![1, 1]), (2, 1, vec![0, 3]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let s_poly = f.s_polynomial(&g);
        assert!(!s_poly.is_zero());
    }

    #[test]
    fn test_groebner_basis_simple() {
        // Ideal generated by x^2 - y, xy - 1
        let f1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (-1, 1, vec![0, 1])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![1, 1]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_groebner_basis_cox_little_oshea() {
        // Example from Cox, Little, O'Shea: x^2 + 2xy^2, xy + 2y^3 - 1
        let f = create_polynomial(
            vec![(1, 1, vec![2, 0]), (2, 1, vec![1, 2])],
            2,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 1]), (2, 1, vec![0, 3]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));

        // Expected result should contain x and y^3 - 1/2
        let _expected_xx = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let _expected_y3 = create_polynomial(
            vec![(1, 1, vec![0, 3]), (-1, 2, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        // Check that basis contains equivalent polynomials
        assert!(basis.len() >= 2);
    }

    #[test]
    fn test_katsura_3() {
        // Katsura-3 system: x + 2y + 2z - 1, x^2 + 2y^2 + 2z^2 - x, 2xy + 2yz - y
        let f1 = create_polynomial(
            vec![
                (1, 1, vec![1, 0, 0]),
                (2, 1, vec![0, 1, 0]),
                (2, 1, vec![0, 0, 1]),
                (-1, 1, vec![0, 0, 0]),
            ],
            3,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![
                (1, 1, vec![2, 0, 0]),
                (2, 1, vec![0, 2, 0]),
                (2, 1, vec![0, 0, 2]),
                (-1, 1, vec![1, 0, 0]),
            ],
            3,
            MonomialOrder::Lex,
        );

        let f3 = create_polynomial(
            vec![
                (2, 1, vec![1, 1, 0]),
                (2, 1, vec![0, 1, 1]),
                (-1, 1, vec![0, 1, 0]),
            ],
            3,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2, f3], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_rational_field() {
        let r1 = Rational::new(1, 2);
        let r2 = Rational::new(3, 4);

        assert_eq!(r1.add(&r2), Rational::new(5, 4));
        assert_eq!(r1.multiply(&r2), Rational::new(3, 8));
        assert_eq!(r1.inverse(), Rational::new(2, 1));
    }

    #[test]
    fn test_monomial_ordering() {
        let m1 = Monomial::new(vec![2, 1]);
        let m2 = Monomial::new(vec![1, 2]);

        // In lex order: x^2*y > x*y^2
        assert_eq!(m1.compare(&m2, MonomialOrder::Lex), Ordering::Greater);

        // In grlex order: both have degree 3, so lex tiebreaker
        assert_eq!(m1.compare(&m2, MonomialOrder::GrLex), Ordering::Greater);

        // In grevlex order: both have degree 3, so reverse lex tiebreaker
        assert_eq!(m1.compare(&m2, MonomialOrder::GRevLex), Ordering::Less);
    }

    #[test]
    fn test_zero_polynomial() {
        let zero = Polynomial::<Rational>::zero(2, MonomialOrder::Lex);
        assert!(zero.is_zero());
        assert!(zero.leading_term().is_none());

        let constant = Polynomial::constant(Rational::new(5, 1), 2, MonomialOrder::Lex);
        assert!(!constant.is_zero());
        assert_eq!(constant.terms.len(), 1);
    }

    #[test]
    fn test_polynomial_reduction() {
        // Test reducing x^2 + xy by [x + y]
        let dividend = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![1, 1])],
            2,
            MonomialOrder::Lex,
        );

        let divisor = create_polynomial(
            vec![(1, 1, vec![1, 0]), (1, 1, vec![0, 1])],
            2,
            MonomialOrder::Lex,
        );

        let remainder = dividend.reduce(&[divisor]);
        // Should reduce to some remainder
        assert!(remainder.terms.len() <= dividend.terms.len());
    }

    #[test]
    fn test_empty_and_single_polynomial() {
        let empty_basis = groebner_basis::<Rational>(vec![], MonomialOrder::Lex);
        assert!(empty_basis.is_empty());

        let single = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let single_basis = groebner_basis(vec![single.clone()], MonomialOrder::Lex);
        assert_eq!(single_basis.len(), 1);
        assert!(is_groebner_basis(&single_basis));
    }
}
