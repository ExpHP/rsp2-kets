#[macro_use]
extern crate criterion;
extern crate rsp2_kets;
extern crate rand;
use criterion::Criterion;
use rsp2_kets::{Ket, Basis};

fn generate_ket(dim: usize) -> Ket {
    let real = (0..dim).map(|_| 0.5 - ::rand::random::<f64>()).collect();
    let imag = (0..dim).map(|_| 0.5 - ::rand::random::<f64>()).collect();
    Ket::new(real, imag)
}

fn generate_basis(num_kets: usize, dim: usize) -> Basis {
    let data = (0..num_kets * dim * 2).map(|_| 0.5 - ::rand::random::<f64>()).collect();
    Basis::new(data, dim)
}

fn bench_orthonormalize_30(c: &mut Criterion) {
    c.bench_function_over_inputs(
        "orthonormalize/30",
        |b, &&n| {
            let input = generate_basis(n, 30);
            b.iter(|| input.orthonormalize() );
        },
        &[16, 64, 256, 1024],
    );
}


fn bench_overlap(c: &mut Criterion) {
    c.bench_function_over_inputs(
        "overlap",
        |b, &&n| {
            let bra = generate_ket(n);
            let ket = generate_ket(n);
            b.iter(|| bra.overlap(&ket) );
        },
        &[16, 64, 256, 1024, 4096, 4096*2, 4096*4, 4096*8, 4096*16],
    );
}

criterion_group!{
    name = benches_overlap;
    config = Criterion::default();
    targets = bench_overlap,
}

criterion_group!{
    name = benches_orthonormalize;
    config = Criterion::default().sample_size(4);
    targets = bench_orthonormalize_30,
}

criterion_main!{
    benches_overlap,
    benches_orthonormalize,
}
