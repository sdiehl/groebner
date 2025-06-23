# Groebner Basis

This is an implementation of the F4 and Buchberger algorithm for computing Groebner bases in Rust.

See [examples/basic.rs](examples/basic.rs) for a usage example.

## Usage

```bash
cargo run --release
```

To run the example, use:
§
```bash
cargo run --example basic
```

## Test Suite

```bash
cargo test
```

The [test suite](SUITE.md) is the full list of known Groebner bases for a variety of large multivariate polynomial systems.

## References

1. Cox, D., Little, J., O'Shea, D. "Ideals, Varieties, and Algorithms"
1. Buchberger, B. "Gröbner Bases: An Algorithmic Method in Polynomial Ideal Theory"

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.