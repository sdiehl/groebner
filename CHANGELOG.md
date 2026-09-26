# Changelog

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
