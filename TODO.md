# TODO for v0.2.1

## 1. Mathematical Structures & Generalizations

- [x] Support for Finite Fields (e.g., GF(p))
  - [x] Add field implementations for finite fields to enable computations over fields of prime order.
- [x] Support for Multivariate Polynomial Rings over Arbitrary Fields
  - Generalize code to allow user-defined field types and provide built-in options.
- [ ] Support for Polynomial Rings with Parameters (Coefficient Rings)
  - [x] One parameter: coefficients in Q(a) via `RationalFunction`.
  - [ ] Several parameters (needs multivariate polynomial GCD).
- [ ] Support for Modules over Polynomial Rings
  - Extend to handle submodules of free modules over polynomial rings (syzygies, free resolutions).

## 2. Algorithmic Optimizations & New Algorithms

- [x] Add criterion benchmarks
- [ ] Implement F4/F5 Algorithms
  - [x] Add sparse F4 implementation over prime fields.
  - [ ] Add faster Groebner basis algorithms (F5) for large systems.
- [x] Parallelization
  - [x] Parallelize batched S-polynomial reduction and F4 row encoding with Rayon.
- [x] Sparse Polynomial Representations
  - [x] Optimize storage and arithmetic for sparse polynomials.
- [x] Improved Reduction Strategies to minimize unnecessary S-polynomial computations
  - [x] Sugar strategy
  - [x] Gebauer–Möller criteria
- [x] Incremental Groebner Basis Computation
  - [x] Allow incremental updates to the basis when new generators are added.

## 3. Usability & API Improvements

- [x] Replace panics with Result-based error handling.
- [x] Allow users to define custom variable order for monomial orders.
- [x] Pretty Printing and LaTeX Output
- [x] Variable Naming and Symbolic Input
  - [x] Allow users to specify variable names and parse polynomials from strings.
