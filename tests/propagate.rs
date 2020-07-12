#[path = "../test_cases.rs"]
mod test_cases;
use test_cases::*;

#[test]
fn propagate() -> sgp4::Result<()> {
    let test_cases: TestCases = toml::from_str(include_str!("../test_cases.toml")).unwrap();
    for test_case in test_cases.list.iter() {
        let constants =
            sgp4::Constants::from_elements_afspc_compatibility_mode(&sgp4::Elements::from_tle(
                None,
                test_case.line1.as_bytes(),
                test_case.line2.as_bytes(),
            )?)?;
        for state in &test_case.states {
            match state {
                State::Ok {
                    time,
                    position,
                    velocity,
                    ..
                } => {
                    let prediction = constants.propagate_afspc_compatibility_mode(*time)?;
                    for index in 0..3 {
                        assert!((position[index] - prediction.position[index]).abs() < 1.0e-6);
                        assert!((velocity[index] - prediction.velocity[index]).abs() < 1.0e-9);
                    }
                }
                State::Err { time, error } => {
                    let prediction = constants.propagate_afspc_compatibility_mode(*time);
                    if let Err(prediction_error) = prediction {
                        assert_eq!(&format!("{}", prediction_error), error);
                    } else {
                        panic!("propagation should have returned an errror");
                    }
                }
            }
        }
    }
    Ok(())
}
