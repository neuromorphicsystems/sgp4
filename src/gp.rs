/// Represents a propagation error caused by orbital elements divergence
#[derive(thiserror::Error, Debug, Clone)]
pub enum Error {
    #[error("The propagated eccentricity ({eccentricity}) is outside the range [0, 1[ {t} minutes after epoch (before adding third-body perturbations)")]
    OutOfRangeEccentricity {
        /// Eccentricity value (unitless)
        eccentricity: f64,

        /// Minutes since epoch
        t: f64,
    },

    #[error("The propagated eccentricity ({eccentricity}) is outside the range [0, 1[ {t} minutes after epoch (after adding third-body perturbations)")]
    OutOfRangePerturbedEccentricity {
        /// Eccentricity value (unitless)
        eccentricity: f64,

        /// Minutes since epoch
        t: f64,
    },

    #[error("The propagated semi-latus rectum is negative {t} minutes after epoch")]
    NegativeSemiLatusRectum {
        /// Minutes since epoch
        t: f64,
    },
}
