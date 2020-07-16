use criterion::{criterion_group, criterion_main, Criterion};
#[path = "../test_cases.rs"]
mod test_cases;
use test_cases::*;

pub fn criterion_benchmark(criterion: &mut Criterion) {
    let test_cases: TestCases = toml::from_str(include_str!("../test_cases.toml")).unwrap();
    criterion.bench_function("propagate all", |b| {
        b.iter(|| {
            let mut predictions = Vec::new();
            for test_case in test_cases.list.iter() {
                let constants = sgp4::Constants::from_elements_afspc_compatibility_mode(
                    &sgp4::Elements::from_tle(
                        None,
                        test_case.line1.as_bytes(),
                        test_case.line2.as_bytes(),
                    )
                    .unwrap(),
                )
                .unwrap();
                for state in &test_case.states {
                    if let State::Ok { time, .. } = state {
                        predictions
                            .push(constants.propagate_afspc_compatibility_mode(*time).unwrap());
                    }
                }
            }
            predictions
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
