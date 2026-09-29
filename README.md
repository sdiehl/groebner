# Groebner Basis

An optimized implementation of the F4 and Buchberger algorithm for computing Groebner bases in Rust. It achieves SOTA performance on the standard benchmarks over prime fields, using parallel sparse linear algebra, and computes rational bases by multi-modular reconstruction.

Examples:

- [Buchberger example](examples/buchberger.rs)
- [F4 example](examples/f4.rs)

## Usage

```bash
cargo add groebner
```

Parse polynomials with named variables under a monomial order and compute a reduced Groebner basis with F4 (over `PrimeField<P>`, `Zp` or `BigRational`):

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
| chandra12 |       2048 | 0.45 s |
| chandra13 |       4096 | 1.82 s |
| cyclic8   |        372 | 0.17 s |
| cyclic9   |       1344 | 8.47 s |
| eco12     |        743 | 0.29 s |
| eco13     |       1465 | 1.13 s |
| eco14     |       2852 | 7.01 s |
| henrion7  |        415 | 0.32 s |
| henrion8  |       2344 | 54.7 s |
| katsura11 |       1050 | 0.62 s |
| katsura12 |       2091 | 3.50 s |
| katsura13 |       4140 | 23.5 s |
| noon9     |       3682 | 1.80 s |
| noon10    |      10273 | 12.6 s |
| reimer7   |        227 | 0.14 s |
| reimer8   |        612 | 2.65 s |

Inputs are in [`benches/data`](benches/data).

## Test Suite

```bash
cargo test
cargo bench
```

The [test suite](SUITE.md) is the full list of known Groebner bases for a variety of large multivariate polynomial systems from several textbooks and some trusted Mathematica generated corpus. Both the Rust algo implementations have to correctly produce the same textbook outputs and Mathematica for all inputs, up to re-ordering.

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
