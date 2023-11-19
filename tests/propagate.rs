#[path = "../test_cases.rs"]
mod test_cases;
use test_cases::*;

#[test]
fn propagate() -> anyhow::Result<()> {
    let test_cases: TestCases = toml::from_str(include_str!("../test_cases.toml")).unwrap();
    for test_case in test_cases.list.iter() {
        #[cfg(feature = "alloc")]
        let element =
            sgp4::Elements::from_tle(None, test_case.line1.as_bytes(), test_case.line2.as_bytes())
                .map_err(|error| anyhow::anyhow!("{error}"))?;

        #[cfg(not(feature = "alloc"))]
        let element =
            sgp4::Elements::from_tle(test_case.line1.as_bytes(), test_case.line2.as_bytes())
                .map_err(|error| anyhow::anyhow!("{error}"))?;

        let constants = sgp4::Constants::from_elements_afspc_compatibility_mode(&element)
            .map_err(|error| anyhow::anyhow!("{error}"))?;
        for state in &test_case.states {
            match state {
                State::Ok {
                    time,
                    position,
                    velocity,
                    ..
                } => {
                    let prediction = constants
                        .propagate_afspc_compatibility_mode(sgp4::MinutesSinceEpoch(*time))
                        .map_err(|error| anyhow::anyhow!("{error}"))?;
                    for index in 0..3 {
                        assert!((position[index] - prediction.position[index]).abs() < 1.0e-6);
                        assert!((velocity[index] - prediction.velocity[index]).abs() < 1.0e-9);
                    }
                }
                State::Err { time, error } => {
                    let prediction = constants
                        .propagate_afspc_compatibility_mode(sgp4::MinutesSinceEpoch(*time));
                    if let Err(prediction_error) = prediction {
                        match prediction_error {
                            sgp4::Error::OutOfRangePerturbedEccentricity { eccentricity: _, t } => {
                                if error != "diverging perturbed eccentricity" {
                                    panic!(
                                        "bad error type (expected {}, got {:?})",
                                        error, prediction_error
                                    );
                                } else {
                                    assert_eq!(*time, t);
                                }
                            }
                            sgp4::Error::NegativeSemiLatusRectum { t } => {
                                if error != "negative semi-latus rectum" {
                                    panic!(
                                        "bad error type (expected {}, got {:?})",
                                        error, prediction_error
                                    );
                                } else {
                                    assert_eq!(*time, t);
                                }
                            }
                            _ => panic!(
                                "bad error type (expected {}, got {:?})",
                                error, prediction_error
                            ),
                        }
                    } else {
                        panic!("propagation should have returned an error");
                    }
                }
            }
        }
    }
    Ok(())
}
