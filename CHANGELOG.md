# Changelog

## Unreleased

- Intern packed F4 monomials in an open addressing table, and renumber and sort
  preprocessing columns in parallel. reimer8 runs about a tenth faster.
- Store F4 basis elements as coefficients and packed term keys, so a preprocessing
  product is one add of the term key and the packed multiplier. Monomials are kept
  only for terms that do not pack or when a trace records the support.
- Move F4 matrix columns out of the plan instead of copying them, and keep only the
  leading terms during final interreduction. Peak memory falls by a tenth to over a
  third, from 8.6 GB to 5.3 GB on henrion8.
- Prune F4 critical pairs with packed lcm keys, testing each new pair against
  lower degree lcms in parallel. The kept pairs match the sequential criterion.
- Search reducers for new preprocessing monomials in parallel, testing divisibility
  on packed keys.
- Run the F4 main loop on a pool thread, so each parallel step starts by work
  stealing rather than waking the pool from outside.
- Echelonize F4 matrices over prime fields by random linear combinations of row
  blocks, reducing about one row per new pivot. A block stops after enough
  consecutive zero combinations to miss a pivot with probability below 2^-40.
  Multi-modular runs keep the deterministic traced elimination.
- Intern F4 matrix monomials as packed exponent keys, in parallel per preprocessing
  level with sharded tables for new monomials.
- Reduce leftover F4 entries during the elimination sweep, skipping the output pass
  for rows that reduce to zero.
- Scatter F4 reducer rows in fixed chunks of eight and store residues as `u16`
  when the prime fits in 16 bits.
- Share coefficients among F4 reducer rows that are multiples of one basis element,
  cutting peak memory by about a fifth.
- Interreduce the final F4 basis with one Macaulay matrix instead of polynomial
  division, for both direct and per-prime modular runs.
- Echelonize new F4 rows in parallel: workers claim free pivot columns
  concurrently, and back substitution reduces every row independently.
- Intern F4 matrix monomials by exponent slice with FxHash, allocating a product
  only when it is a new column.
- Skip zero cells and clear only the touched span in dense modular row reduction.
- Store modular matrix rows as `u32` residues whenever the prime fits in 32 bits.
- Look up F4 reducers in active elements sorted by length, and filter
  Gebauer-Moller pair criteria by divisibility mask without allocating.
- Breaking: `SparseRow::columns` is now `Vec<u32>`.
- Breaking: `F4Field` methods take known pivots as `Reducers`.

## 0.4.1 (2026-09-27)

- Verify exact cofactor certificates for every reconstructed basis polynomial.
- Build `RationalFunction` on `polycore::RatFunc`.
- Add `Frac` coefficients in several parameters.
- Gate `Frac` behind the optional `parameters` feature.
- Add `PolynomialRing::with_parameters` and `ParseCoefficient::parameter`.
- Lower MSRV to Rust 1.88.

## 0.4.0 (2026-09-27)

- Rebase fields, monomials, orders, and sparse polynomials on `polycore`.
- Use parallel multi-modular F4 automatically over `BigRational`, with leading-term
  voting, support alignment, persistent CRT with separate recovery, parallel
  coefficient reconstruction, and a fresh full modular check. Validate candidates
  individually instead of computing a batch.
- Add `groebner_basis_f4_rational(polys, certify)` for optional exact checks and
  `f4::groebner_basis_f4_direct` for the direct coefficient-field algorithm.
- Learn F4 matrix plans at the first prime and replay only independent rows at later
  primes. Changed pivots or supports fall back to a full run; repeated failures
  replace the trace. A full modular check independently validates reconstruction.
- Intern matrix monomials, filter divisibility with degree masks, and reduce row
  blocks in parallel with reusable buffers and bounded deferred reduction.
- Parse parentheses, polynomial powers, and division by constant expressions through
  `polycore::Ring`, retaining grouped numbers, implicit products, and formatting.
- Correct corpus fingerprints to omit terms that vanish modulo the reference prime.
- Check a fixed 256-system corpus tier on every push and pull request, with and
  without Rayon. Missing fixtures, mismatches, crashes, and timeouts fail CI.

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
