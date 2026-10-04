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

Enable the `mimalloc` feature for maximal performance. It installs mimalloc as the global allocator,
so leave it off if your binary already sets one.

```toml
groebner = { version = "0.5", features = ["mimalloc"] }
```

## Benchmarks

Reduced basis with `groebner_basis_f4` over GF(32003) in GRevLex, on Apple M5 silicon (10 cores).

| System    | Basis size |   Time |
| --------- | ---------: | -----: |
| chandra12 |       2048 | 0.41 s |
| chandra13 |       4096 | 1.42 s |
| cyclic8   |        372 | 0.18 s |
| cyclic9   |       1344 | 4.98 s |
| eco12     |        743 | 0.25 s |
| eco13     |       1465 | 0.84 s |
| eco14     |       2852 | 4.35 s |
| henrion7  |        415 | 0.35 s |
| henrion8  |       2344 | 24.1 s |
| katsura11 |       1050 | 0.47 s |
| katsura12 |       2091 | 2.12 s |
| katsura13 |       4140 | 11.8 s |
| noon9     |       3682 | 1.54 s |
| noon10    |      10273 | 8.64 s |
| reimer7   |        227 | 0.14 s |
| reimer8   |        612 | 1.54 s |

With a runtime modulus through `Zp`, same settings, by size of the prime.

| System    |  32003 | 2^31 - 1 | 2^62 - 57 | 2^63 - 25 |
| --------- | -----: | -------: | --------: | --------: |
| cyclic8   | 0.18 s |   0.21 s |    0.25 s |    0.30 s |
| katsura11 | 0.48 s |   0.60 s |    0.74 s |    0.97 s |
| noon9     | 1.56 s |   1.67 s |    1.80 s |    1.99 s |
| reimer8   | 1.59 s |   1.75 s |    2.00 s |    2.52 s |

Over the rationals with `BigRational` coefficients, same settings.

| System    | Basis size |   Time |
| --------- | ---------: | -----: |
| cyclic8   |        372 | 1.68 s |
| eco12     |        743 | 1.37 s |
| katsura9  |        272 | 0.23 s |
| katsura10 |        537 | 1.68 s |
| noon9     |       3682 | 3.22 s |
| reimer7   |        227 | 0.67 s |

Inputs are in [`benches/data`](benches/data).

## Test Suite

```bash
cargo test
cargo bench --features mimalloc
```

The [test suite](SUITE.md) is the full list of known Gröbner bases for a variety of large multivariate polynomial systems from several textbooks and some trusted Mathematica generated corpus. Both the Rust algo implementations have to correctly produce the same textbook outputs and Mathematica for all inputs, up to re-ordering.

## License

MIT Licensed. Copyright 2024-2026 Stephen Diehl. See [LICENSE](LICENSE) for details.

