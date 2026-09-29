# Gröbner Basis

An optimized Rust implementation of the F4 and Buchberger algorithms for computing Gröbner bases. It achieves SOTA performance on standard benchmarks over prime fields and the rationals, using parallel sparse linear algebra, SIMD row reduction, and multi-modular rational reconstruction.

Examples:

- [Buchberger example](examples/buchberger.rs)
- [F4 example](examples/f4.rs)

## Usage

```bash
cargo add groebner
```

Parse polynomials with named variables under a monomial order and compute a reduced Gröbner basis with F4 (over `PrimeField<P>`, `Zp` or `BigRational`):

```rust
use groebner::{groebner_basis_f4, MonomialOrder, PolynomialRing, PrimeField};

let ring = PolynomialRing::<PrimeField<32003>>::new(["x", "y"], MonomialOrder::GRevLex)?;
let polys = ring.parse_many("x^2 - y; x*y - 1")?;
let basis = groebner_basis_f4(polys, true)?;
for p in &basis {
    println!("{}", ring.format(p)?);
}
# Ok::<(), Box<dyn std::error::Error>>(())
```

See the [API documentation](https://docs.rs/groebner) for Buchberger, runtime primes, ideals, lift
certificates, monomial orders, and parametric coefficients.

## Benchmarks

Reduced basis with `groebner_basis_f4` over GF(32003) in GRevLex, on Apple M5 silicon (10 cores).

| System    | Basis size |   Time |
| --------- | ---------: | -----: |
| chandra12 |       2048 | 0.42 s |
| chandra13 |       4096 | 1.52 s |
| cyclic8   |        372 | 0.17 s |
| cyclic9   |       1344 | 5.39 s |
| eco12     |        743 | 0.26 s |
| eco13     |       1465 | 0.87 s |
| eco14     |       2852 | 4.76 s |
| henrion7  |        415 | 0.32 s |
| henrion8  |       2344 | 26.4 s |
| katsura11 |       1050 | 0.48 s |
| katsura12 |       2091 | 2.42 s |
| katsura13 |       4140 | 13.3 s |
| noon9     |       3682 | 1.62 s |
| noon10    |      10273 | 10.1 s |
| reimer7   |        227 | 0.13 s |
| reimer8   |        612 | 1.68 s |

With a runtime modulus through `Zp`, same settings, by size of the prime.

| System    |  32003 | 2^31 - 1 | 2^62 - 57 |
| --------- | -----: | -------: | --------: |
| cyclic8   | 0.18 s |   0.21 s |    0.25 s |
| katsura11 | 0.51 s |   0.65 s |    0.81 s |
| noon9     | 1.69 s |   1.81 s |    1.98 s |
| reimer8   | 1.68 s |   1.84 s |    2.21 s |

Over the rationals with `BigRational` coefficients, same settings.

| System    | Basis size |   Time |
| --------- | ---------: | -----: |
| cyclic8   |        372 | 1.94 s |
| eco12     |        743 | 1.57 s |
| katsura9  |        272 | 0.23 s |
| katsura10 |        537 | 2.03 s |
| noon9     |       3682 | 3.92 s |
| reimer7   |        227 | 0.71 s |

Inputs are in [`benches/data`](benches/data).

## Test Suite

```bash
cargo test
cargo bench
```

The [test suite](SUITE.md) is the full list of known Gröbner bases for a variety of large multivariate polynomial systems from several textbooks and some trusted Mathematica generated corpus. Both the Rust algo implementations have to correctly produce the same textbook outputs and Mathematica for all inputs, up to re-ordering.

## References

The classic papers on this topic:

1. Cox, D., Little, J., O'Shea, D. "Ideals, Varieties, and Algorithms"
1. Buchberger, B. "Gröbner Bases: An Algorithmic Method in Polynomial Ideal Theory"
1. Giovini, A., Mora, T., Niesi, G., Robbiano, L., & Traverso, C. (1991, June). “One sugar cube, please” or selection strategies in the Buchberger algorithm. In Proceedings of the 1991 international symposium on Symbolic and algebraic computation (pp. 49-54).
1. Gebauer, R., & Möller, H. M. (1988). On an installation of Buchberger's algorithm. Journal of Symbolic computation, 6(2-3), 275-286.
1. Roune, B. H., & Stillman, M. (2012, July). Practical Gröbner basis computation. In Proceedings of the 37th International Symposium on Symbolic and Algebraic Computation (pp. 203-210).
1. Faugère, J. C. (1999). A new efficient algorithm for computing Gröbner bases (F4). Journal of Pure and Applied Algebra, 139(1-3), 61-88.
1. Faugère, J. C., Gianni, P., Lazard, D., & Mora, T. (1993). Efficient computation of zero-dimensional Gröbner bases by change of ordering. Journal of Symbolic Computation, 16(4), 329-344.
1. Monagan, M., & Pearce, R. (2015). A compact parallel implementation of F4. In Proceedings of PASCO 2015 (pp. 95-100).

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
