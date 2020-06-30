use crate::model;
use crate::third_body;
use crate::tle;

#[derive(Debug, Clone)]
pub struct Error {
    message: String,
}

impl Error {
    pub fn new(message: &str) -> Error {
        Error {
            message: message.to_owned(),
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "{}", self.message)
    }
}

impl std::error::Error for Error {}

impl From<tle::Error> for Error {
    fn from(error: tle::Error) -> Self {
        Error::new(&error.to_string())
    }
}

pub type Result<T> = std::result::Result<T, Error>;

pub struct Prediction {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
}

pub enum Elliptic {
    No {},
    Yes { k12: f64, k13: f64 },
}

pub enum HighAltitude {
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
        k11: f64,
        elliptic: Elliptic,
    },
}

pub enum Resonance {
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

pub enum Resonant {
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

pub enum Method {
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

pub struct Orbit {
    pub inclination: f64,
    pub right_ascension: f64,
    pub eccentricity: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
    pub mean_motion: f64,
}

pub struct Constants<'a> {
    pub geopotential: &'a model::Geopotential,
    pub right_ascension_dot: f64,
    pub argument_of_perigee_dot: f64,
    pub mean_anomaly_dot: f64,
    pub c1: f64,
    pub c4: f64,
    pub k0: f64,
    pub k1: f64,
    pub method: Method,
    pub orbit_0: Orbit,
}
