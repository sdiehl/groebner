# Groebner Basis

This is an implementation of the F4 and Buchberger algorithm for computing Groebner bases in Rust.

See [examples/basic.rs](examples/basic.rs) for a usage example.

## Usage

To use in your project:

```bash
cargo add groebner
```

## Test Suite

To run the example, use:

```bash
cargo run --example basic
```

```bash
cargo test
```

The [test suite](SUITE.md) is the full list of known Groebner bases for a variety of large multivariate polynomial systems.

## References

1. Cox, D., Little, J., O'Shea, D. "Ideals, Varieties, and Algorithms"
1. Buchberger, B. "Gröbner Bases: An Algorithmic Method in Polynomial Ideal Theory"
1. Giovini, A., Mora, T., Niesi, G., Robbiano, L., & Traverso, C. (1991, June). “One sugar cube, please” or selection strategies in the Buchberger algorithm. In Proceedings of the 1991 international symposium on Symbolic and algebraic computation (pp. 49-54).
1. Gebauer, R., & Möller, H. M. (1988). On an installation of Buchberger's algorithm. Journal of Symbolic computation, 6(2-3), 275-286.
1. Roune, B. H., & Stillman, M. (2012, July). Practical Gröbner basis computation. In Proceedings of the 37th International Symposium on Symbolic and Algebraic Computation (pp. 203-210).

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
