# TODO for v0.2.0

## 1. Mathematical Structures & Generalizations
- [ ] Support for Finite Fields (e.g., GF(p))
    - Add field implementations for finite fields to enable computations over fields of prime order.
    - Use https://github.com/zkcrypto/ff for finite fields.
- [ ] Support for Multivariate Polynomial Rings over Arbitrary Fields
    - Generalize code to allow user-defined field types and provide built-in options.
- [ ] Support for Polynomial Rings with Parameters (Coefficient Rings)
    - Allow coefficients to be polynomials themselves (nested polynomial rings).
- [ ] Support for Modules over Polynomial Rings
    - Extend to handle submodules of free modules over polynomial rings (syzygies, free resolutions).

## 2. Algorithmic Optimizations & New Algorithms
- [ ] Add criterion benchmarks
- [ ] Implement F4/F5 Algorithms
    - Add faster Groebner basis algorithms (F4, F5) for large systems.
- [ ] Parallelization
    - Parallelize S-polynomial reduction and other steps for performance.
- [ ] Sparse Polynomial Representations
    - Optimize storage and arithmetic for sparse polynomials.
- [ ] Improved Reduction Strategies
    - Implement strategies like sugar strategy, Gebauer–Möller criteria to minimize unnecessary S-polynomial computations.
- [ ] Incremental Groebner Basis Computation
    - Allow incremental updates to the basis when new generators are added.

## 3. Usability & API Improvements
- [x] Replace panics with Result-based error handling.
- [ ] Custom Monomial Orders
    - Allow users to define custom monomial orders.
- [ ] Pretty Printing and LaTeX Output
    - Add pretty-printing for polynomials and LaTeX export.
- [ ] Variable Naming and Symbolic Input
    - Allow users to specify variable names and parse polynomials from strings.