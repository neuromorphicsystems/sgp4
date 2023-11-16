/// Represents a propagation error caused by orbital elements divergence
#[derive(Debug, Clone)]
pub enum Error {
    OutOfRangeEccentricity {
        /// Eccentricity value (unitless)
        eccentricity: f64,

        /// Minutes since epoch
        t: f64,
    },

    OutOfRangePerturbedEccentricity {
        /// Eccentricity value (unitless)
        eccentricity: f64,

        /// Minutes since epoch
        t: f64,
    },

    NegativeSemiLatusRectum {
        /// Minutes since epoch
        t: f64,
    },
}

impl core::fmt::Display for Error {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::OutOfRangeEccentricity { eccentricity, t } => formatter.write_fmt(
                core::format_args!(
                    "The propagated eccentricity ({}) is outside the range [0, 1[ {} minutes after epoch (before adding third-body perturbations)",
                    eccentricity,
                    t,
                )
            ),
            Error::OutOfRangePerturbedEccentricity { eccentricity, t } => formatter.write_fmt(
                core::format_args!(
                    "The propagated eccentricity ({}) is outside the range [0, 1[ {} minutes after epoch (after adding third-body perturbations)",
                    eccentricity,
                    t,
                )
            ),
            Error::NegativeSemiLatusRectum { t } => formatter.write_fmt(
                core::format_args!("The propagated semi-latus rectum is negative {} minutes after epoch", t)
            )
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}
