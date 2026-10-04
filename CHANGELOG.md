# Changelog

## 0.5.1 (2026-10-05)

- Take Lehmer reconstruction and mixed radix CRT from `polycore` 0.1.3.
- Widen rational replay rounds automatically when the first image is quick.
- Reduce rows modulo primes above 2^31 with Shoup multiplication.
- Reduce rows modulo primes above 2^62 with `polycore` `MulBy`.
- Accumulate Monte Carlo combinations with NEON widening multiply-accumulate.
- Accumulate lane replay rows and Monte Carlo combinations with AVX2 on x86_64.
- Add katsura14, eco15, and cyclic10 benchmark inputs.
- Add opt-in `mimalloc` feature installing mimalloc as the global allocator.
- Build tests at `opt-level = 3` with overflow checks.

## 0.5.0 (2026-09-29)

- Hash packed monomials from the top product bits to avoid probe clustering.
- Sweep four random combinations per pass in the Monte Carlo echelon.
- Balance symbolic preprocessing levels by row length across threads.
- Add RationalOptions to fix the rational replay batch width.
- Abort failed rational reconstructions early and gcd small operands in place.
- Split lane replay rows finely so small rounds use every thread.
- Replay the final interreduction four primes at a time.
- Reconstruct rational coefficients with Lehmer gcd and half gcd.
- Validate rational candidates against a held-out replayed image.
- Reduce rational candidates modulo check primes in parallel.
- Accumulate lane replay rows with NEON widening multiply-accumulate.
- Reuse planned replay columns for basis elements with unchanged support.
- Accumulate CRT images as machine word mixed radix digits.
- Validate rational F4 untraced, overlapping the first learned image.
- Replay rational F4 traces at four primes per pass.
- Reduce large-prime sweeps with branchless conditional subtraction.
- Speed up rational F4 reconstruction by a third or more.
- Intern packed F4 monomials in an open addressing table.
- Renumber and sort preprocessing columns in parallel.
- Store F4 basis elements as coefficients and packed term keys.
- Move F4 matrix columns out of the plan instead of copying.
- Keep only leading terms during final interreduction, cutting peak memory.
- Prune F4 critical pairs in parallel with packed lcm keys.
- Search reducers for preprocessing monomials in parallel on packed keys.
- Run the F4 main loop on a pool thread.
- Echelonize F4 matrices by random linear combinations of row blocks.
- Keep deterministic traced elimination for multi-modular runs.
- Intern F4 matrix monomials as packed keys, sharded per level.
- Reduce leftover F4 entries during the elimination sweep.
- Scatter F4 reducer rows in fixed chunks of eight.
- Store residues as `u16` when the prime fits in 16 bits.
- Share coefficients among F4 reducer rows from one basis element.
- Interreduce the final F4 basis with one Macaulay matrix.
- Echelonize new F4 rows in parallel with concurrent pivot claims.
- Intern F4 matrix monomials by exponent slice with FxHash.
- Skip zero cells in dense modular row reduction.
- Store modular matrix rows as `u32` residues when possible.
- Look up F4 reducers in active elements sorted by length.
- Filter Gebauer-Moller pair criteria by divisibility mask without allocating.
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
