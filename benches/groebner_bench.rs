#![allow(clippy::expect_used)]

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
#[cfg(feature = "parallel")]
use groebner::groebner_basis_parallel;
use groebner::{
    MonomialOrder, Polynomial, PolynomialRing, PrimeField, groebner_basis, groebner_basis_f4,
};
use num_rational::BigRational;
use std::hint::black_box;

type F32003 = PrimeField<32003>;

const CYCLIC7: &str = include_str!("data/cyclic7.txt");
const KATSURA7: &str = include_str!("data/katsura7.txt");

fn vars(nvars: usize) -> Vec<String> {
    (0..nvars).map(|index| format!("x{index}")).collect()
}

fn parse_prime_system(input: &str, nvars: usize, limit: usize) -> Vec<Polynomial<F32003>> {
    let variables = vars(nvars);
    let ring = PolynomialRing::<F32003>::new(variables, MonomialOrder::GrLex)
        .expect("benchmark ring should be valid");
    input
        .lines()
        .take(limit)
        .map(|line| ring.parse(line).expect("benchmark polynomial should parse"))
        .collect()
}

fn parse_rational_system(input: &str, nvars: usize, limit: usize) -> Vec<Polynomial<BigRational>> {
    let variables = vars(nvars);
    let ring = PolynomialRing::<BigRational>::new(variables, MonomialOrder::GrLex)
        .expect("benchmark ring should be valid");
    input
        .lines()
        .take(limit)
        .map(|line| ring.parse(line).expect("benchmark polynomial should parse"))
        .collect()
}

fn bench_axf4_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("axf4_parse");
    group.bench_function("cyclic7_gf32003", |b| {
        b.iter(|| black_box(parse_prime_system(CYCLIC7, 7, 7)));
    });
    group.bench_function("katsura7_gf32003", |b| {
        b.iter(|| black_box(parse_prime_system(KATSURA7, 8, 8)));
    });
    group.finish();
}

fn bench_groebner_small_axf4_families(c: &mut Criterion) {
    let mut group = c.benchmark_group("groebner_axf4_subsystems");
    group.sample_size(10);

    for (name, input, nvars, limit) in [
        ("cyclic7_first2", CYCLIC7, 7, 2),
        ("katsura7_first3", KATSURA7, 8, 3),
    ] {
        group.bench_with_input(
            BenchmarkId::new("buchberger_gf32003", name),
            &(input, nvars, limit),
            |b, (input, nvars, limit)| {
                b.iter_batched(
                    || parse_prime_system(input, *nvars, *limit),
                    |polys| {
                        black_box(
                            groebner_basis(polys, true).expect("GF(p) benchmark should compute"),
                        );
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );

        group.bench_with_input(
            BenchmarkId::new("f4_gf32003", name),
            &(input, nvars, limit),
            |b, (input, nvars, limit)| {
                b.iter_batched(
                    || parse_prime_system(input, *nvars, *limit),
                    |polys| {
                        black_box(
                            groebner_basis_f4(polys, true)
                                .expect("F4 GF(p) benchmark should compute"),
                        );
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );

        #[cfg(feature = "parallel")]
        group.bench_with_input(
            BenchmarkId::new("parallel_buchberger_gf32003", name),
            &(input, nvars, limit),
            |b, (input, nvars, limit)| {
                b.iter_batched(
                    || parse_prime_system(input, *nvars, *limit),
                    |polys| {
                        black_box(
                            groebner_basis_parallel(polys, true)
                                .expect("parallel GF(p) benchmark should compute"),
                        );
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );

        group.bench_with_input(
            BenchmarkId::new("rational", name),
            &(input, nvars, limit),
            |b, (input, nvars, limit)| {
                b.iter_batched(
                    || parse_rational_system(input, *nvars, *limit),
                    |polys| {
                        black_box(
                            groebner_basis(polys, true).expect("rational benchmark should compute"),
                        );
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }

    group.finish();
}

fn bench_f4_full_systems(c: &mut Criterion) {
    let mut group = c.benchmark_group("f4_axf4_full");
    group.sample_size(10);
    for (name, input, nvars) in [("cyclic7", CYCLIC7, 7), ("katsura7", KATSURA7, 8)] {
        group.bench_with_input(
            BenchmarkId::new("f4_gf32003", name),
            &(input, nvars),
            |b, (input, nvars)| {
                b.iter_batched(
                    || parse_prime_system(input, *nvars, usize::MAX),
                    |polys| {
                        black_box(
                            groebner_basis_f4(polys, true).expect("F4 benchmark should compute"),
                        );
                    },
                    criterion::BatchSize::SmallInput,
                );
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_axf4_parsing,
    bench_groebner_small_axf4_families,
    bench_f4_full_systems
);
criterion_main!(benches);
