# Groebner Basis

This is an implementation of the F4 and Buchberger algorithm for computing Groebner bases in Rust.

Examples:

- [Buchberger example](examples/buchberger.rs)
- [F4 example](examples/f4.rs)

## Usage

To use in your project:

```bash
cargo add groebner
```

Create polynomials from strings with an explicit variable order:

```rust
use groebner::{groebner_basis, MonomialOrder, PolynomialRing};
use num_rational::BigRational;

let ring = PolynomialRing::<BigRational>::new(["z3", "z1", "z2"], MonomialOrder::Lex)?;
let f1 = ring.parse("z1^2 - z2")?;
let f2 = ring.parse("z1*z2 - 1")?;
let basis = groebner_basis(vec![f1, f2], ring.order(), true)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

For modular computations over a machine prime:

```rust
use groebner::{groebner_basis, MonomialOrder, PolynomialRing, PrimeField};

type F32003 = PrimeField<32003>;

let ring = PolynomialRing::<F32003>::new(["x", "y"], MonomialOrder::GrLex)?;
let f1 = ring.parse("x^2 - y")?;
let f2 = ring.parse("x*y - 1")?;
let basis = groebner_basis(vec![f1, f2], ring.order(), true)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Test Suite

To run the examples, use:

```bash
cargo run --example buchberger
cargo run --example f4
```

To build the API documentation locally:

```bash
cargo doc --no-deps
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
