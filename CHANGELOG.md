# Changelog

## 0.4.0 (unreleased)

- Rebase fields, monomials, orders, and sparse polynomials on `polycore`.
- Use parallel multi-modular F4 automatically over `BigRational`, with leading-term
  voting, support alignment, bounded CRT retries, and a fresh full modular check.
- Add `groebner_basis_f4_rational(polys, certify)` for optional exact checks and
  `f4::groebner_basis_f4_direct` for the direct coefficient-field algorithm.
- Learn F4 matrix plans at the first prime and replay them at later primes. Changed
  pivots, supports, or zero-row dependencies fall back to a full run.
- Intern matrix monomials, filter divisibility with degree masks, and reduce row
  blocks in parallel with reusable buffers and bounded deferred reduction.
- Parse parentheses, polynomial powers, and division by constant expressions through
  `polycore::Ring`, retaining grouped numbers, implicit products, and formatting.
- Correct corpus fingerprints to omit terms that vanish modulo the reference prime.

### Migration from 0.3

`Zp`, `PrimeField`, `MonomialOrder`, and `Polynomial<F>` remain exported names for
`Fp`, `Gf`, `Order`, and `Poly<F>`. These are the polycore types, without conversion.

- `PrimeField` now takes a `u64` const modulus. Its `value()` returns `u64`.
- `Field` uses `+`, `-`, `*`, `/`, and unary `-`. Import `num_traits::{Zero, One}`
  for concrete-type identities; `Field::inverse` still returns `Option`.
- Terms are `(Monomial, coefficient)` tuples. Replace `Term::new(c, m)` with `(m, c)`;
  read `.0` and `.1` in place of `.monomial` and `.coefficient`.
- Use `m.exps()` instead of the public `m.exponents` field. `Order::Weighted` is
  now a tuple variant; prefer `Order::weighted`.
- Import `PolynomialExt` for `s_polynomial`, `make_monic`, and compatibility helpers.
  Polycore's `reduce` and `divide` return values directly. The old fallible
  shortest-reducer versions are `normal_form` and `divide_with_remainder`.
- Import `MonomialExt` for the old monomial method names, or use polycore's `var`,
  `quo`, multiplication operator, and `order.compare(&a, &b)`.
- Remove calls to the deprecated `groebner_basis_f4_mod`; use `groebner_basis_f4`.
- Rational F4 returns a reduced basis even when `canonicalize` is false. The default
  reconstruction check is probabilistic; pass `certify = true` for exact Buchberger
  and input-reduction checks.


## 0.3.1 (2026-09-26)

- Add `RationalFunction` for coefficients in Q(a).
- Add `PolynomialRing::with_parameter` for parsing and printing parameters.
- Add `specialize` to substitute a value for the parameter.
- Add defaulted `ParseCoefficient` hooks for parameters and printing.
- Add `Polynomial::divide` returning quotients and remainder.
- Add `Ideal::lift` and `LiftBasis` for membership certificates.
- Add `Ideal::generators` for the original generating polynomials.
- Add `verify_lift` to check certificates by ring arithmetic.
- Add cyclic8-9, eco12, katsura8-10, and noon9 benchmarks.
- Drop the direct `num-bigint` dependency.
- Move to edition 2024 with minimum Rust 1.98.

## 0.3.0 (2026-09-25)

- Add generic `groebner_basis_f4` over any `F4Field`, including `BigRational`.
- Add dense modular F4 path for `PrimeField` and `Zp`.
- Add `Zp` for runtime primes up to 64 bits.
- Add `Ideal` with membership, normal forms, and elimination.
- Add `Ideal` dimension, standard monomials, and multiplication matrices.
- Add `Ideal::radical_contains` and `Ideal::change_order`.
- Add FGLM change of order for zero-dimensional ideals.
- Add weighted, block, and elimination monomial orders.
- Add `PolynomialRing::with_modulus` and `ParseCoefficient`.
- Breaking: `groebner_basis` and variants take order from the polynomials.
- Deprecate `groebner_basis_f4_mod` in favor of `groebner_basis_f4`.

## 0.2.1

- Add `format_latex` for LaTeX polynomial output.
- Add `groebner_basis_incremental`.
- Add `groebner_basis_parallel` and `is_groebner_basis_parallel` via Rayon.
- Add default `parallel` feature.
- Validate test suite against Mathematica-generated corpus.
- Bump criterion to 0.8.2.

## 0.2.0 (2026-07-07)

- Add sparse F4 algorithm via `groebner_basis_f4_mod`.
- Add `PrimeField<P>` for modular arithmetic over machine primes.
- Make `PolynomialRing` generic over the coefficient field.
- Add criterion benchmarks with cyclic7 and katsura7 systems.
- Add `buchberger` and `f4` examples.

## 0.1.2 (2026-07-07)

- Add `PolynomialRing` with string parsing and formatting (#2).
- Variable list defines the lexicographic variable order.
- Expand stress tests.

## 0.1.1 (2025-06-25)

- Add sugar selection strategy.
- Add Gebauer-Moller pair criteria via `filter_gm_pairs`.
- Add `groebner_basis_with_strategy` and `SelectionStrategy`.
- Pass selection strategy by reference.
- Add stress tests.
- Improve documentation and clippy lints.

## 0.1.0 (2025-06-23)

- Initial release.
- Buchberger algorithm over rational coefficients.
- Lex, GrLex, and GRevLex monomial orders.
- `GroebnerError` replaces panics and unwraps.
- Textbook test suite documented in `SUITE.md`.
