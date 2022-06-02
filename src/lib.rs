//! This crate implements the SGP4 algorithm for satellite propagation.
//!
//! It also provides methods to parse Two-Line Element Set (TLE) and Orbit Mean-Elements Message (OMM) data.
//!
//! A UTF-8 transcript of the mathematical expressions that compose SGP4
//! can be found in the repository [https://github.com/neuromorphicsystems/sgp4](https://github.com/neuromorphicsystems/sgp4).
//!
//! # Example
//!
//! The following standalone program downloads the lastest Galileo OMMs from Celestrak
//! and predicts the satellites' positions and velocities after 12 h and 24 h.
//!
//! ```
//! fn main() -> sgp4::Result<()> {
//!     match ureq::get("https://celestrak.com/NORAD/elements/gp.php")
//!         .query("GROUP", "galileo")
//!         .query("FORMAT", "json")
//!         .call()
//!     {
//!         Ok(response) => {
//!             let elements_group: Vec<sgp4::Elements> = response.into_json()?;
//!             for elements in &elements_group {
//!                 println!("{}", elements.object_name.as_ref().unwrap());
//!                 let constants = sgp4::Constants::from_elements(elements)?;
//!                 for hours in &[12, 24] {
//!                     println!("    t = {} min", hours * 60);
//!                     let prediction = constants.propagate((hours * 60) as f64)?;
//!                     println!("        r = {:?} km", prediction.position);
//!                     println!("        ṙ = {:?} km.s⁻¹", prediction.velocity);
//!                 }
//!             }
//!             Ok(())
//!         }
//!         Err(ureq::Error::Status(code, response)) => Err(sgp4::Error::new(format!(
//!             "network error {}: {}",
//!             code,
//!             response.into_string()?
//!         ))),
//!         Err(error) => Err(sgp4::Error::new(error.to_string())),
//!     }
//! }
//! ```
//! More examples can be found in the repository [https://github.com/neuromorphicsystems/sgp4/tree/master/examples](https://github.com/neuromorphicsystems/sgp4/tree/master/examples).
//!

mod deep_space;
mod gp;
mod model;
mod near_earth;
mod propagator;
mod third_body;

pub use deep_space::ResonanceState;
pub use gp::parse_2les;
pub use gp::parse_3les;
pub use gp::Classification;
pub use gp::Elements;
pub use gp::Error;
pub use gp::Result;
pub use model::afspc_epoch_to_sidereal_time;
pub use model::iau_epoch_to_sidereal_time;
pub use model::Geopotential;
pub use model::WGS72;
pub use model::WGS84;
pub use propagator::Constants;
pub use propagator::Orbit;
pub use propagator::Prediction;

impl Orbit {
    /// Creates a new Brouwer orbit representation from Kozai elements
    ///
    /// If the Kozai orbital elements are obtained from a TLE or OMM,
    /// the convenience function [sgp4::Constants::from_elements](struct.Constants.html#method.from_elements)
    /// can be used instead of manually mapping the `Elements` fields to the `Constants::new` parameters.
    ///
    /// # Arguments
    ///
    /// * `geopotential` - The model of Earth gravity to use in the conversion
    /// * `inclination` - Angle between the equator and the orbit plane in rad
    /// * `right_ascension` - Angle between vernal equinox and the point where the orbit crosses the equatorial plane in rad
    /// * `eccentricity` - The shape of the orbit
    /// * `argument_of_perigee` - Angle between the ascending node and the orbit's point of closest approach to the earth in rad
    /// * `mean_anomaly` - Angle of the satellite location measured from perigee in rad
    /// * `kozai_mean_motion` - Mean orbital angular velocity in rad.min⁻¹ (Kozai convention)
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let elements = sgp4::Elements::from_tle(
    ///     Some("ISS (ZARYA)".to_owned()),
    ///     "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992".as_bytes(),
    ///     "2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008".as_bytes(),
    /// )?;
    /// let orbit_0 = sgp4::Orbit::from_kozai_elements(
    ///     &sgp4::WGS84,
    ///     elements.inclination * (std::f64::consts::PI / 180.0),
    ///     elements.right_ascension * (std::f64::consts::PI / 180.0),
    ///     elements.eccentricity,
    ///     elements.argument_of_perigee * (std::f64::consts::PI / 180.0),
    ///     elements.mean_anomaly * (std::f64::consts::PI / 180.0),
    ///     elements.mean_motion * (std::f64::consts::PI / 720.0),
    /// )?;
    /// #     Ok(())
    /// # }
    /// ```
    pub fn from_kozai_elements(
        geopotential: &Geopotential,
        inclination: f64,
        right_ascension: f64,
        eccentricity: f64,
        argument_of_perigee: f64,
        mean_anomaly: f64,
        kozai_mean_motion: f64,
    ) -> Result<Self> {
        if kozai_mean_motion <= 0.0 {
            Err(Error::new(
                "the Kozai mean motion must be positive".to_owned(),
            ))
        } else {
            let mean_motion = {
                // a₁ = (kₑ / n₀)²ᐟ³
                let a1 = (geopotential.ke / kozai_mean_motion).powf(2.0 / 3.0);

                //      3      3 cos²I₀
                // p₀ = - J₂ -----------
                //      4    (1 − e₀²)³ᐟ²
                let p0 = 0.75 * geopotential.j2 * (3.0 * inclination.cos().powi(2) - 1.0)
                    / (1.0 - eccentricity.powi(2)).powf(3.0 / 2.0);

                // 𝛿₁ = p₀ / a₁²
                let d1 = p0 / a1.powi(2);

                // 𝛿₀ = p₀ / (a₁ (1 - ¹/₃ 𝛿₁ - 𝛿₁² - ¹³⁴/₈₁ 𝛿₁³))²
                let d0 = p0
                    / (a1 * (1.0 - d1.powi(2) - d1 * (1.0 / 3.0 + 134.0 * d1.powi(2) / 81.0)))
                        .powi(2);

                //         n₀
                // n₀" = ------
                //       1 + 𝛿₀
                kozai_mean_motion / (1.0 + d0)
            };
            if mean_motion <= 0.0 {
                Err(Error::new(
                    "the Brouwer mean motion must be positive".to_owned(),
                ))
            } else {
                Ok(propagator::Orbit {
                    inclination,
                    right_ascension,
                    eccentricity,
                    argument_of_perigee,
                    mean_anomaly,
                    mean_motion,
                })
            }
        }
    }
}

impl<'a> Constants<'a> {
    /// Initializes a new propagator from epoch quantities
    ///
    /// If the orbital elements are obtained from a TLE or OMM,
    /// the convenience function [sgp4::Constants::from_elements](struct.Constants.html#method.from_elements)
    /// can be used instead of manually mapping the `Elements` fields to the `Constants::new` parameters.
    ///
    /// # Arguments
    ///
    /// * `geopotential` - The model of Earth gravity to use in the conversion
    /// * `epoch_to_sidereal_time` - The function to use to convert the J2000 epoch to sidereal time
    /// * `epoch` - The number of years since UTC 1 January 2000 12h00 (J2000)
    /// * `drag_term` - The radiation pressure coefficient in earth radii⁻¹ (B*)
    /// * `orbit_0` - The Brouwer orbital elements at epoch
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let elements = sgp4::Elements::from_tle(
    ///     Some("ISS (ZARYA)".to_owned()),
    ///     "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992".as_bytes(),
    ///     "2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008".as_bytes(),
    /// )?;
    /// let constants = sgp4::Constants::new(
    ///     &sgp4::WGS84,
    ///     sgp4::iau_epoch_to_sidereal_time,
    ///     elements.epoch(),
    ///     elements.drag_term,
    ///     sgp4::Orbit::from_kozai_elements(
    ///         &sgp4::WGS84,
    ///         elements.inclination * (std::f64::consts::PI / 180.0),
    ///         elements.right_ascension * (std::f64::consts::PI / 180.0),
    ///         elements.eccentricity,
    ///         elements.argument_of_perigee * (std::f64::consts::PI / 180.0),
    ///         elements.mean_anomaly * (std::f64::consts::PI / 180.0),
    ///         elements.mean_motion * (std::f64::consts::PI / 720.0),
    ///     )?,
    /// )?;
    /// #     Ok(())
    /// # }
    /// ```
    pub fn new(
        geopotential: &'a Geopotential,
        epoch_to_sidereal_time: impl Fn(f64) -> f64,
        epoch: f64,
        drag_term: f64,
        orbit_0: propagator::Orbit,
    ) -> Result<Self> {
        if orbit_0.eccentricity < 0.0 || orbit_0.eccentricity >= 1.0 {
            Err(Error::new(
                "the eccentricity must be in the range [0, 1[".to_owned(),
            ))
        } else {
            // p₁ = cos I₀
            let p1 = orbit_0.inclination.cos();

            // p₂ = 1 − e₀²
            let p2 = 1.0 - orbit_0.eccentricity.powi(2);

            // k₆ = 3 p₁² - 1
            let k6 = 3.0 * p1.powi(2) - 1.0;

            // a₀" = (kₑ / n₀")²ᐟ³
            let a0 = (geopotential.ke / orbit_0.mean_motion).powf(2.0 / 3.0);

            // p₃ = a₀" (1 - e₀)
            let p3 = a0 * (1.0 - orbit_0.eccentricity);
            let (s, p6) = {
                // p₄ = aₑ (p₃ - 1)
                let p4 = geopotential.ae * (p3 - 1.0);

                // p₅ = │ 20      if p₄ < 98
                //      │ p₄ - 78 if 98 ≤ p₄ < 156
                //      │ 78      otherwise
                let p5 = if p4 < 98.0 {
                    20.0
                } else if p4 < 156.0 {
                    p4 - 78.0
                } else {
                    78.0
                };
                (
                    // s = p₅ / aₑ + 1
                    p5 / geopotential.ae + 1.0,
                    // p₆ = ((120 - p₅) / aₑ)⁴
                    ((120.0 - p5) / geopotential.ae).powi(4),
                )
            };

            // ξ = 1 / (a₀" - s)
            let xi = 1.0 / (a0 - s);

            // p₇ = p₆ ξ⁴
            let p7 = p6 * xi.powi(4);

            // η = a₀" e₀ ξ
            let eta = a0 * orbit_0.eccentricity * xi;

            // p₈ = |1 - η²|
            let p8 = (1.0 - eta.powi(2)).abs();

            // p₉ = p₇ / p₈⁷ᐟ²
            let p9 = p7 / p8.powf(3.5);

            // C₁ = B* p₉ n₀" (a₀" (1 + ³/₂ η² + e₀ η (4 + η²))
            //      + ³/₈ J₂ ξ k₆ (8 + 3 η² (8 + η²)) / p₈)
            let c1 = drag_term
                * (p9
                    * orbit_0.mean_motion
                    * (a0
                        * (1.0
                            + 1.5 * eta.powi(2)
                            + orbit_0.eccentricity * eta * (4.0 + eta.powi(2)))
                        + 0.375 * geopotential.j2 * xi / p8
                            * k6
                            * (8.0 + 3.0 * eta.powi(2) * (8.0 + eta.powi(2)))));

            // p₁₀ = (a₀" p₂)⁻²
            let p10 = 1.0 / (a0 * p2).powi(2);

            // β₀ = p₂¹ᐟ²
            let b0 = p2.sqrt();

            // p₁₁ = ³/₂ J₂ p₁₀ n₀"
            let p11 = 1.5 * geopotential.j2 * p10 * orbit_0.mean_motion;

            // p₁₂ = ¹/₂ p₁₁ J₂ p₁₀
            let p12 = 0.5 * p11 * geopotential.j2 * p10;

            // p₁₃ = - ¹⁵/₃₂ J₄ p₁₀² n₀"
            let p13 = -0.46875 * geopotential.j4 * p10.powi(2) * orbit_0.mean_motion;

            // p₁₄ = - p₁₁ p₁ + (¹/₂ p₁₂ (4 - 19 p₁²) + 2 p₁₃ (3 - 7 p₁²)) p₁
            let p14 = -p11 * p1
                + (0.5 * p12 * (4.0 - 19.0 * p1.powi(2)) + 2.0 * p13 * (3.0 - 7.0 * p1.powi(2)))
                    * p1;

            // k₁₄ = - ¹/₂ p₁₁ (1 - 5 p₁²) + ¹/₁₆ p₁₂ (7 - 114 p₁² + 395 p₁⁴)
            let k14 = -0.5 * p11 * (1.0 - 5.0 * p1.powi(2))
                + 0.0625 * p12 * (7.0 - 114.0 * p1.powi(2) + 395.0 * p1.powi(4))
                + p13 * (3.0 - 36.0 * p1.powi(2) + 49.0 * p1.powi(4));

            // p₁₅ = n₀" + ¹/₂ p₁₁ β₀ k₆ + ¹/₁₆ p₁₂ β₀ (13 - 78 p₁² + 137 p₁⁴)
            let p15 = orbit_0.mean_motion
                + 0.5 * p11 * b0 * k6
                + 0.0625 * p12 * b0 * (13.0 - 78.0 * p1.powi(2) + 137.0 * p1.powi(4));

            // C₄ = 2 B* n₀" p₉ a₀" p₂ (
            //      η (2 + ¹/₂ η²)
            //      + e₀ (¹/₂ + 2 η²)
            //      - J₂ ξ / (a p₈) (-3 k₆ (1 - 2 e₀ η + η² (³/₂ - ¹/₂ e₀ η))
            //      + ³/₄ (1 - p₁²) (2 η² - e₀ η (1 + η²)) cos 2 ω₀)
            let c4 = drag_term
                * (2.0
                    * orbit_0.mean_motion
                    * p9
                    * a0
                    * p2
                    * (eta * (2.0 + 0.5 * eta.powi(2))
                        + orbit_0.eccentricity * (0.5 + 2.0 * eta.powi(2))
                        - geopotential.j2 * xi / (a0 * p8)
                            * (-3.0
                                * k6
                                * (1.0 - 2.0 * orbit_0.eccentricity * eta
                                    + eta.powi(2) * (1.5 - 0.5 * orbit_0.eccentricity * eta))
                                + 0.75
                                    * (1.0 - p1.powi(2))
                                    * (2.0 * eta.powi(2)
                                        - orbit_0.eccentricity * eta * (1.0 + eta.powi(2)))
                                    * (2.0 * orbit_0.argument_of_perigee).cos())));

            // k₀ = - ⁷/₂ p₂ p₁₁ p₁ C₁
            let k0 = 3.5 * p2 * (-p11 * p1) * c1;

            // k₁ = ³/₂ C₁
            let k1 = 1.5 * c1;

            if orbit_0.mean_motion > 2.0 * std::f64::consts::PI / 225.0 {
                Ok(near_earth::constants(
                    geopotential,
                    drag_term,
                    orbit_0,
                    p1,
                    a0,
                    s,
                    xi,
                    eta,
                    c1,
                    c4,
                    k0,
                    k1,
                    k6,
                    k14,
                    p2,
                    p3,
                    p7,
                    p9,
                    p14,
                    p15,
                ))
            } else {
                Ok(deep_space::constants(
                    geopotential,
                    epoch_to_sidereal_time,
                    epoch,
                    orbit_0,
                    p1,
                    a0,
                    c1,
                    b0,
                    c4,
                    k0,
                    k1,
                    k14,
                    p2,
                    p14,
                    p15,
                ))
            }
        }
    }

    /// Initializes a new propagator from an `Elements` object
    ///
    /// This is the recommended method to initialize a propagator from a TLE or OMM.
    /// The WGS84 model, the IAU sidereal time expression and the accurate UTC to J2000 expression are used.
    ///
    /// # Arguments
    ///
    /// * `elements` - Orbital elements and drag term parsed from a TLE or OMM
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let constants = sgp4::Constants::from_elements(
    ///     &sgp4::Elements::from_tle(
    ///         Some("ISS (ZARYA)".to_owned()),
    ///         "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992".as_bytes(),
    ///         "2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008".as_bytes(),
    ///     )?,
    /// )?;
    /// #     Ok(())
    /// # }
    /// ```
    pub fn from_elements(elements: &Elements) -> Result<Self> {
        Constants::new(
            &WGS84,
            iau_epoch_to_sidereal_time,
            elements.epoch(),
            elements.drag_term,
            Orbit::from_kozai_elements(
                &WGS84,
                elements.inclination * (std::f64::consts::PI / 180.0),
                elements.right_ascension * (std::f64::consts::PI / 180.0),
                elements.eccentricity,
                elements.argument_of_perigee * (std::f64::consts::PI / 180.0),
                elements.mean_anomaly * (std::f64::consts::PI / 180.0),
                elements.mean_motion * (std::f64::consts::PI / 720.0),
            )?,
        )
    }

    /// Initializes a new propagator from an `Elements` object
    ///
    /// This method should be used if compatibility with the AFSPC implementation is needed.
    /// The WGS72 model, the AFSPC sidereal time expression and the AFSPC UTC to J2000 expression are used.
    ///
    /// # Arguments
    ///
    /// * `elements` - Orbital elements and drag term parsed from a TLE or OMM
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let constants = sgp4::Constants::from_elements_afspc_compatibility_mode(
    ///     &sgp4::Elements::from_tle(
    ///         Some("ISS (ZARYA)".to_owned()),
    ///         "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992".as_bytes(),
    ///         "2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008".as_bytes(),
    ///     )?,
    /// )?;
    /// #     Ok(())
    /// # }
    /// ```
    pub fn from_elements_afspc_compatibility_mode(elements: &Elements) -> Result<Self> {
        Constants::new(
            &WGS72,
            afspc_epoch_to_sidereal_time,
            elements.epoch_afspc_compatibility_mode(),
            elements.drag_term,
            Orbit::from_kozai_elements(
                &WGS72,
                elements.inclination * (std::f64::consts::PI / 180.0),
                elements.right_ascension * (std::f64::consts::PI / 180.0),
                elements.eccentricity,
                elements.argument_of_perigee * (std::f64::consts::PI / 180.0),
                elements.mean_anomaly * (std::f64::consts::PI / 180.0),
                elements.mean_motion * (std::f64::consts::PI / 720.0),
            )?,
        )
    }

    /// Returns the initial deep space resonance integrator state
    ///
    /// For most orbits, SGP4 propagation is stateless.
    /// That is, predictions at different times are always calculated from the epoch quantities.
    /// No calculations are saved by propagating to successive times sequentially.
    ///
    /// However, resonant deep space orbits (geosynchronous or Molniya) use an integrator
    /// to estimate the resonance effects of Earth gravity, with a 720 min time step. If the propagation
    /// times are monotonic, a few operations per prediction can be saved by re-using the integrator state.
    ///
    /// The high-level API `Constants::propagate` re-initializes the state with each propagation for simplicity.
    /// `Constants::initial_state` and `Constants::propagate_from_state` can be used together
    /// to speed up resonant deep space satellite propagation.
    /// For non-deep space or non-resonant orbits, their behavior is identical to `Constants::propagate`.
    ///
    /// See `Constants::propagate_from_state` for an example.
    pub fn initial_state(&self) -> Option<ResonanceState> {
        match &self.method {
            propagator::Method::NearEarth { .. } => None,
            propagator::Method::DeepSpace { resonant, .. } => match resonant {
                propagator::Resonant::No { .. } => None,
                propagator::Resonant::Yes { lambda_0, .. } => {
                    Some(ResonanceState::new(self.orbit_0.mean_motion, *lambda_0))
                }
            },
        }
    }

    /// Calculates the SGP4 position and velocity predictions
    ///
    /// This is an advanced API which results in marginally faster propagation than `Constants::propagate` in some cases
    /// (see `Constants::initial_state` for details), at the cost of added complexity for the user.
    ///
    /// The propagation times must be monotonic if the same resonance state is used repeatedly.
    /// The `afspc_compatibility_mode` makes a difference only if the satellite is on a Lyddane deep space orbit
    /// (period greater than 225 min and inclination smaller than 0.2 rad).
    ///
    /// # Arguments
    ///
    /// * `t` - The number of minutes since epoch (can be positive, negative or zero)
    /// * `state` - The deep space propagator state returned by `Constants::initial_state`
    /// * `afspc_compatibility_mode` - Set to true if compatibility with the AFSPC implementation is needed
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let elements = sgp4::Elements::from_tle(
    ///     Some("MOLNIYA 1-36".to_owned()),
    ///     "1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813".as_bytes(),
    ///     "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656".as_bytes(),
    /// )?;
    /// let constants = sgp4::Constants::from_elements(&elements)?;
    /// let mut state = constants.initial_state();
    /// for days in 0..7 {
    ///     println!("t = {} min", days * 60 * 24);
    ///     let prediction =
    ///         constants.propagate_from_state((days * 60 * 24) as f64, state.as_mut(), false)?;
    ///     println!("    r = {:?} km", prediction.position);
    ///     println!("    ṙ = {:?} km.s⁻¹", prediction.velocity);
    /// }
    /// #     Ok(())
    /// # }
    /// ```
    #[allow(clippy::many_single_char_names)]
    pub fn propagate_from_state(
        &self,
        t: f64,
        state: Option<&mut ResonanceState>,
        afspc_compatibility_mode: bool,
    ) -> Result<Prediction> {
        // p₂₂ = Ω₀ + Ω̇ t + k₀ t²
        let p22 = self.orbit_0.right_ascension + self.right_ascension_dot * t + self.k0 * t.powi(2);

        // p₂₃ = ω₀ + ω̇ t
        let p23 = self.orbit_0.argument_of_perigee + self.argument_of_perigee_dot * t;
        let (orbit, a, p32, p33, p34, p35, p36) = match &self.method {
            propagator::Method::NearEarth {
                a0,
                k2,
                k3,
                k4,
                k5,
                k6,
                high_altitude,
            } => {
                assert!(
                    state.is_none(),
                    "state must be None with a near-earth propagator",
                );
                self.near_earth_orbital_elements(
                    *a0,
                    *k2,
                    *k3,
                    *k4,
                    *k5,
                    *k6,
                    high_altitude,
                    t,
                    p22,
                    p23,
                )
            }
            propagator::Method::DeepSpace {
                eccentricity_dot,
                inclination_dot,
                solar_perturbations,
                lunar_perturbations,
                resonant,
            } => self.deep_space_orbital_elements(
                *eccentricity_dot,
                *inclination_dot,
                solar_perturbations,
                lunar_perturbations,
                resonant,
                state,
                t,
                p22,
                p23,
                afspc_compatibility_mode,
            ),
        }?;

        // p₃₇ = 1 / (a (1 - e²))
        let p37 = 1.0 / (a * (1.0 - orbit.eccentricity.powi(2)));

        // aₓₙ = e cos ω
        let axn = orbit.eccentricity * orbit.argument_of_perigee.cos();

        // aᵧₙ = e sin ω + p₃₇ p₃₂
        let ayn = orbit.eccentricity * orbit.argument_of_perigee.sin() + p37 * p32;

        // p₃₈ = M + ω + p₃₇ p₃₅ aₓₙ rem 2π
        let p38 = (orbit.mean_anomaly + orbit.argument_of_perigee + p37 * p35 * axn)
            % (2.0 * std::f64::consts::PI);

        // (E + ω)₀ = p₃₈
        let mut ew = p38;
        for _ in 0..10 {
            //             p₃₈ - aᵧₙ cos (E + ω)ᵢ + aₓₙ sin (E + ω)ᵢ - (E + ω)ᵢ
            // Δ(E + ω)ᵢ = ---------------------------------------------------
            //                   1 - cos (E + ω)ᵢ aₓₙ - sin (E + ω)ᵢ aᵧₙ
            let delta = (p38 - ayn * ew.cos() + axn * ew.sin() - ew)
                / (1.0 - ew.cos() * axn - ew.sin() * ayn);

            if delta.abs() < 1.0e-12 {
                break;
            }

            // (E + ω)ᵢ₊₁ = (E + ω)ᵢ + Δ(E + ω)ᵢ|[-0.95, 0.95]
            ew += if delta < -0.95 {
                -0.95
            } else if delta > 0.95 {
                0.95
            } else {
                delta
            };
        }

        // p₃₉ = aₓₙ² + aᵧₙ²
        let p39 = axn.powi(2) + ayn.powi(2);

        // pₗ = a (1 - p₃₉)
        let pl = a * (1.0 - p39);
        if pl < 0.0 {
            Err(Error::new("negative semi-latus rectum".to_owned()))
        } else {
            // p₄₀ = aₓₙ sin(E + ω) - aᵧₙ cos(E + ω)
            let p40 = axn * ew.sin() - ayn * ew.cos();

            // r = a (1 - aₓₙ cos(E + ω) + aᵧₙ sin(E + ω))
            let r = a * (1.0 - (axn * ew.cos() + ayn * ew.sin()));

            // ṙ = a¹ᐟ² p₄₀ / r
            let r_dot = a.sqrt() * p40 / r;

            // β = (1 - p₃₉)¹ᐟ²
            let b = (1.0 - p39).sqrt();

            // p₄₁ = p₄₀ / (1 + β)
            let p41 = p40 / (1.0 + b);

            // p₄₂ = a / r (sin(E + ω) - aᵧₙ - aₓₙ p₄₁)
            let p42 = a / r * (ew.sin() - ayn - axn * p41);

            // p₄₃ = a / r (cos(E + ω) - aₓₙ + aᵧₙ p₄₁)
            let p43 = a / r * (ew.cos() - axn + ayn * p41);

            //           p₄₂
            // u = tan⁻¹ ---
            //           p₄₃
            let u = p42.atan2(p43);

            // p₄₄ = 2 p₄₃ p₄₂
            let p44 = 2.0 * p43 * p42;

            // p₄₅ = 1 - 2 p₄₂²
            let p45 = 1.0 - 2.0 * p42.powi(2);

            // p₄₆ = (¹/₂ J₂ / pₗ) / pₗ
            let p46 = 0.5 * self.geopotential.j2 / pl / pl;

            // rₖ = r (1 - ³/₂ p₄₆ β p₃₆) + ¹/₂ (¹/₂ J₂ / pₗ) p₃₃ p₄₅
            let rk = r * (1.0 - 1.5 * p46 * b * p36)
                + 0.5 * (0.5 * self.geopotential.j2 / pl) * p33 * p45;

            // uₖ = u - ¹/₄ p₄₆ p₃₄ p₄₄
            let uk = u - 0.25 * p46 * p34 * p44;

            // Iₖ = I + ³/₂ p₄₆ cos I sin I p₄₅
            let inclination_k = orbit.inclination
                + 1.5 * p46 * orbit.inclination.cos() * orbit.inclination.sin() * p45;

            // Ωₖ = Ω + ³/₂ p₄₆ cos I p₄₄
            let right_ascension_k =
                orbit.right_ascension + 1.5 * p46 * orbit.inclination.cos() * p44;

            // ṙₖ = ṙ + n (¹/₂ J₂ / pₗ) p₃₃ / kₑ
            let rk_dot = r_dot
                - orbit.mean_motion * (0.5 * self.geopotential.j2 / pl) * p33 * p44
                    / self.geopotential.ke;

            // rḟₖ = pₗ¹ᐟ² / r + n (¹/₂ J₂ / pₗ) (p₃₃ p₄₅ + ³/₂ p₃₆) / kₑ
            let rfk_dot = pl.sqrt() / r
                + orbit.mean_motion * (0.5 * self.geopotential.j2 / pl) * (p33 * p45 + 1.5 * p36)
                    / self.geopotential.ke;

            // u₀ = - sin Ωₖ cos Iₖ sin uₖ + cos Ωₖ cos uₖ
            let u0 = -right_ascension_k.sin() * inclination_k.cos() * uk.sin()
                + right_ascension_k.cos() * uk.cos();
            // u₁ = cos Ωₖ cos Iₖ sin uₖ + sin Ωₖ cos uₖ
            let u1 = right_ascension_k.cos() * inclination_k.cos() * uk.sin()
                + right_ascension_k.sin() * uk.cos();
            // u₂ = sin Iₖ sin uₖ
            let u2 = inclination_k.sin() * uk.sin();
            Ok(Prediction {
                position: [
                    // r₀ = rₖ u₀ aₑ
                    rk * u0 * self.geopotential.ae,
                    // r₁ = rₖ u₁ aₑ
                    rk * u1 * self.geopotential.ae,
                    // r₂ = rₖ u₂ aₑ
                    rk * u2 * self.geopotential.ae,
                ],
                velocity: [
                    // ṙ₀ = (ṙₖ u₀ + rḟₖ (- sin Ωₖ cos Iₖ cos uₖ - cos Ωₖ sin uₖ)) aₑ kₑ / 60
                    (rk_dot * u0
                        + rfk_dot
                            * (-right_ascension_k.sin() * inclination_k.cos() * uk.cos()
                                - right_ascension_k.cos() * uk.sin()))
                        * (self.geopotential.ae * self.geopotential.ke / 60.0),
                    // ṙ₁ = (ṙₖ u₁ + rḟₖ (cos Ωₖ cos Iₖ cos uₖ - sin Ωₖ sin uₖ)) aₑ kₑ / 60
                    (rk_dot * u1
                        + rfk_dot
                            * (right_ascension_k.cos() * inclination_k.cos() * uk.cos()
                                - right_ascension_k.sin() * uk.sin()))
                        * (self.geopotential.ae * self.geopotential.ke / 60.0),
                    // ṙ₂ = (ṙₖ u₂ + rḟₖ (sin Iₖ cos uₖ)) aₑ kₑ / 60
                    (rk_dot * u2 + rfk_dot * (inclination_k.sin() * uk.cos()))
                        * (self.geopotential.ae * self.geopotential.ke / 60.0),
                ],
            })
        }
    }

    /// Calculates the SGP4 position and velocity predictions
    ///
    /// This is the recommended method to propagate epoch orbital elements.
    ///
    /// # Arguments
    /// `t` - The number of minutes since epoch (can be positive, negative or zero)
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let constants = sgp4::Constants::from_elements_afspc_compatibility_mode(
    ///     &sgp4::Elements::from_tle(
    ///         Some("ISS (ZARYA)".to_owned()),
    ///         "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992".as_bytes(),
    ///         "2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008".as_bytes(),
    ///     )?,
    /// )?;
    /// let prediction = constants.propagate(60.0 * 24.0);
    /// #     Ok(())
    /// # }
    /// ```
    pub fn propagate(&self, t: f64) -> Result<Prediction> {
        self.propagate_from_state(t, self.initial_state().as_mut(), false)
    }

    /// Calculates the SGP4 position and velocity predictions
    ///
    /// This method should be used if compatibility with the AFSPC implementation is needed.
    /// Its behavior is different from `Constants::propagate`
    /// only if the satellite is on a Lyddane deep space orbit
    /// (period greater than 225 min and inclination smaller than 0.2 rad).
    ///
    /// # Arguments
    /// `t` - The number of minutes since epoch (can be positive, negative or zero)
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let constants = sgp4::Constants::from_elements_afspc_compatibility_mode(
    ///     &sgp4::Elements::from_tle(
    ///         Some("ISS (ZARYA)".to_owned()),
    ///         "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992".as_bytes(),
    ///         "2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008".as_bytes(),
    ///     )?,
    /// )?;
    /// let prediction = constants.propagate_afspc_compatibility_mode(60.0 * 24.0);
    /// #     Ok(())
    /// # }
    /// ```
    pub fn propagate_afspc_compatibility_mode(&self, t: f64) -> Result<Prediction> {
        self.propagate_from_state(t, self.initial_state().as_mut(), true)
    }
}
