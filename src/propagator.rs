use crate::model;
use crate::third_body;
use serde::{Deserialize, Serialize};

/// Predicted satellite position and velocity after SGP4 propagation
///
/// The position and velocity are given in the True Equator, Mean Equinox (TEME) of epoch reference frame.
#[derive(Serialize, Deserialize, PartialEq, Debug,)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub struct Prediction {
    /// The three position components (x, y, z) in km
    pub position: [f64; 3],

    /// The three velocity components (x, y, z) in km.s⁻¹
    pub velocity: [f64; 3],
}

/// The Brouwer orbital elements
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub struct Orbit {
    /// Angle between the equator and the orbit plane in rad
    pub inclination: f64,

    /// Angle between vernal equinox and the point where the orbit crosses the equatorial plane in rad
    pub right_ascension: f64,

    /// Shape of the orbit
    pub eccentricity: f64,

    /// Angle between the ascending node and the orbit's point of closest approach to the earth in rad
    pub argument_of_perigee: f64,

    /// Angle of the satellite location measured from perigee in rad
    pub mean_anomaly: f64,

    /// Mean number of orbits per day in rad.min⁻¹
    pub mean_motion: f64,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub(crate) enum Elliptic {
    No {},
    Yes { k11: f64, k12: f64, k13: f64 },
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub(crate) enum HighAltitude {
    No {},
    Yes {
        c5: f64,
        d2: f64,
        d3: f64,
        d4: f64,
        eta: f64,
        k7: f64,
        k8: f64,
        k9: f64,
        k10: f64,
        elliptic: Elliptic,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub(crate) enum Resonance {
    OneDay {
        dr1: f64,
        dr2: f64,
        dr3: f64,
    },
    HalfDay {
        d2201: f64,
        d2211: f64,
        d3210: f64,
        d3222: f64,
        d4410: f64,
        d4422: f64,
        d5220: f64,
        d5232: f64,
        d5421: f64,
        d5433: f64,
        k14: f64,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub(crate) enum Resonant {
    No {
        a0: f64,
    },
    Yes {
        lambda_0: f64,
        lambda_dot_0: f64,
        sidereal_time_0: f64,
        resonance: Resonance,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub(crate) enum Method {
    NearEarth {
        a0: f64,
        k2: f64,
        k3: f64,
        k4: f64,
        k5: f64,
        k6: f64,
        high_altitude: HighAltitude,
    },
    DeepSpace {
        eccentricity_dot: f64,
        inclination_dot: f64,
        solar_perturbations: third_body::Perturbations,
        lunar_perturbations: third_body::Perturbations,
        resonant: Resonant,
    },
}

/// Propagator variables calculated from epoch quantities and used during propagation
///
/// Constants can be initialized from general perturbation elements.
/// They are not mutated during propagation, which means they can
/// be used by different threads in parallel
/// (for example to generate predictions at different times).
#[derive(Serialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub struct Constants<'a> {
    pub(crate) geopotential: &'a model::Geopotential,
    pub(crate) right_ascension_dot: f64,
    pub(crate) argument_of_perigee_dot: f64,
    pub(crate) mean_anomaly_dot: f64,
    pub(crate) c1: f64,
    pub(crate) c4: f64,
    pub(crate) k0: f64,
    pub(crate) k1: f64,
    pub(crate) method: Method,
    pub(crate) orbit_0: Orbit,
}
