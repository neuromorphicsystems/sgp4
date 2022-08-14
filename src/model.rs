use serde::{Deserialize, Serialize};

/// Model of the Earth radius and gravitational field
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub struct Geopotential {
    /// Equatorial radius of the earth in km
    // aₑ
    pub ae: f64,

    /// square root of earth's gravitational parameter in earth radii³ min⁻²
    // kₑ
    pub ke: f64,

    /// un-normalised second zonal harmonic
    // J₂
    pub j2: f64,

    /// un-normalised third zonal harmonic
    // J₃
    pub j3: f64,

    /// un-normalised fourth zonal harmonic
    // J₄
    pub j4: f64,
}

/// The geopotential model recommended by the IAU
///
/// This model is recommended to propagate orbits.
pub const WGS84: Geopotential = Geopotential {
    ae: 6378.137,
    ke: 0.07436685316871385,
    j2: 0.00108262998905,
    j3: -0.00000253215306,
    j4: -0.00000161098761,
};

/// The geopotential model used in the AFSPC implementation
///
/// This model should be used if compatibility with the AFSPC implementation is needed.
pub const WGS72: Geopotential = Geopotential {
    ae: 6378.135,
    ke: 0.07436691613317342,
    j2: 0.001082616,
    j3: -0.00000253881,
    j4: -0.00000165597,
};

/// Converts an epoch to sidereal time using the IAU expression
///
/// This is the recommended method to calculate the sidereal time.
///
/// # Arguments
///
/// * `epoch` - Years since UTC 1 January 2000 12h00
pub fn iau_epoch_to_sidereal_time(epoch: f64) -> f64 {
    // c₂₀₀₀ = y₂₀₀₀ / 100
    let c2000 = epoch / 100.0;

    // θ₀ = ¹/₂₄₀ (π / 180) (- 6.2 × 10⁻⁶ c₂₀₀₀³ + 0.093104 c₂₀₀₀²
    //      + (876600 × 3600 + 8640184.812866) c₂₀₀₀ + 67310.54841) mod 2π
    ((-6.2e-6 * c2000.powi(3)
        + 0.093104 * c2000.powi(2)
        + (876600.0 * 3600.0 + 8640184.812866) * c2000
        + 67310.54841)
        * (std::f64::consts::PI / 180.0)
        / 240.0)
        .rem_euclid(2.0 * std::f64::consts::PI)
}

/// Converts an epoch to sidereal time using the AFSPC expression
///
/// This function should be used if compatibility with the AFSPC implementation is needed.
///
/// # Arguments
///
/// * `epoch` - Years since UTC 1 January 2000 12h00
pub fn afspc_epoch_to_sidereal_time(epoch: f64) -> f64 {
    // d₁₉₇₀ = 365.25 (y₂₀₀₀ + 30) + 1
    let d1970 = (epoch + 30.0) * 365.25 + 1.0;

    // θ₀ = 1.7321343856509374 + 1.72027916940703639 × 10⁻² ⌊t₁₉₇₀ + 10⁻⁸⌋
    //      + (1.72027916940703639 × 10⁻² + 2π) (t₁₉₇₀ - ⌊t₁₉₇₀ + 10⁻⁸⌋)
    //      + 5.07551419432269442 × 10⁻¹⁵ t₁₉₇₀² mod 2π
    #[allow(clippy::excessive_precision)]
    (1.7321343856509374
        + 1.72027916940703639e-2 * (d1970 + 1.0e-8).floor()
        + (1.72027916940703639e-2 + 2.0 * std::f64::consts::PI)
            * (d1970 - (d1970 + 1.0e-8).floor())
        + d1970.powi(2) * 5.07551419432269442e-15)
        .rem_euclid(2.0 * std::f64::consts::PI)
}
