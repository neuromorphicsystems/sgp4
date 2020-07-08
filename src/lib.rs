mod deep_space;
pub mod model;
mod near_earth;
mod propagator;
mod third_body;
pub mod tle;

pub use propagator::Constants;
pub use propagator::Error;
pub use propagator::Orbit;
pub use propagator::Prediction;
pub use propagator::Result;

impl Orbit {
    pub fn from_kozai_elements(
        geopotential: &model::Geopotential,
        inclination: f64,
        right_ascension: f64,
        eccentricity: f64,
        argument_of_perigee: f64,
        mean_anomaly: f64,
        kozai_mean_motion: f64,
    ) -> Result<Self> {
        if kozai_mean_motion <= 0.0 {
            Err(Error::new("the Kozai mean motion must be positive"))
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
                Err(Error::new("the Brouwer mean motion must be positive"))
            } else {
                Ok(Orbit {
                    inclination: inclination,
                    right_ascension: right_ascension,
                    eccentricity: eccentricity,
                    argument_of_perigee: argument_of_perigee,
                    mean_anomaly: mean_anomaly,
                    mean_motion: mean_motion,
                })
            }
        }
    }
}

// geopotential: the gravity model to use in calculations
// epoch: years since UTC 1 January 2000 12h00 t₀
// drag_term: the radiation pressure coefficient B*, in earth radii⁻¹
// inclination_0: the angle between the equator and the orbit plane I₀, in rad
// right_ascension: the angle between vernal equinox and the point where
//                  the orbit crosses the equatorial plane Ω₀, in rad
// eccentricity_0: the shape of the orbit e₀
// argument_of_perigee: the angle between the ascending node and the orbit's
//                      point of closest approach to the earth ω₀, in rad
// mean_anomaly: the angle of the satellite location measured from perigee M₀, in rad
// mean_motion: mean number of orbits per day (Kozai mean motion) n₀, in rad.min⁻¹
impl<'a> Constants<'a> {
    pub fn new(
        geopotential: &'a model::Geopotential,
        epoch_to_sidereal_time: impl Fn(f64) -> f64,
        epoch: f64,
        drag_term: f64,
        orbit_0: Orbit,
    ) -> Result<Self> {
        if orbit_0.eccentricity < 0.0 || orbit_0.eccentricity >= 1.0 {
            Err(Error::new("the eccentricity must be in the range [0, 1["))
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

            if orbit_0.mean_motion > 2.0 * model::PI / 225.0 {
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

    pub fn from_tle(tle: &tle::Tle) -> Result<Self> {
        Constants::new(
            &model::WGS84,
            model::iau_epoch_to_sidereal_time,
            tle.epoch(),
            tle.drag_term,
            Orbit::from_kozai_elements(
                &model::WGS72,
                tle.inclination * (model::PI / 180.0),
                tle.right_ascension * (model::PI / 180.0),
                tle.eccentricity,
                tle.argument_of_perigee * (model::PI / 180.0),
                tle.mean_anomaly * (model::PI / 180.0),
                tle.mean_motion * (model::PI / 720.0),
            )?,
        )
    }

    pub fn from_tle_afspc_compatibility_mode(tle: &tle::Tle) -> Result<Self> {
        Constants::new(
            &model::WGS72,
            model::afspc_epoch_to_sidereal_time,
            tle.epoch_afspc_compatibility_mode(),
            tle.drag_term,
            Orbit::from_kozai_elements(
                &model::WGS72,
                tle.inclination * (model::PI / 180.0),
                tle.right_ascension * (model::PI / 180.0),
                tle.eccentricity,
                tle.argument_of_perigee * (model::PI / 180.0),
                tle.mean_anomaly * (model::PI / 180.0),
                tle.mean_motion * (model::PI / 720.0),
            )?,
        )
    }

    pub fn initial_state(&self) -> Option<deep_space::ResonanceState> {
        match &self.method {
            propagator::Method::NearEarth { .. } => None,
            propagator::Method::DeepSpace { resonant, .. } => match resonant {
                propagator::Resonant::No { .. } => None,
                propagator::Resonant::Yes { lambda_0, .. } => Some(
                    deep_space::ResonanceState::new(self.orbit_0.mean_motion, *lambda_0),
                ),
            },
        }
    }

    pub fn propagate_from_state(
        &self,
        t: f64,
        state: Option<&mut deep_space::ResonanceState>,
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
        let p38 = (orbit.mean_anomaly + orbit.argument_of_perigee + p37 * p35 * axn) % (2.0 * model::PI);

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
            Err(Error::new("negative semi-latus rectum"))
        } else {
            // p₄₀ = aₓₙ cos(E + ω) + aᵧₙ sin(E + ω)
            let p40 = axn * ew.cos() + ayn * ew.sin();

            // p₄₁ = aₓₙ sin(E + ω) - aᵧₙ cos(E + ω)
            let p41 = axn * ew.sin() - ayn * ew.cos();

            // r = a (1 - p₄₀)
            let r = a * (1.0 - p40);

            // ṙ = a¹ᐟ² p₄₁ / r
            let r_dot = a.sqrt() * p41 / r;

            // β = (1 - p₃₉)¹ᐟ²
            let b = (1.0 - p39).sqrt();

            // p₄₂ = p₄₁ / (1 + β)
            let p42 = p41 / (1.0 + b);

            // p₄₃ = a / r (sin(E + ω) - aᵧₙ - aₓₙ p₄₂)
            let p43 = a / r * (ew.sin() - ayn - axn * p42);

            // p₄₄ = a / r (cos(E + ω) - aₓₙ + aᵧₙ p₄₂)
            let p44 = a / r * (ew.cos() - axn + ayn * p42);

            //           p₄₃
            // u = tan⁻¹ ---
            //           p₄₄
            let u = p43.atan2(p44);

            // p₄₅ = 2 p₄₄ p₄₃
            let p45 = 2.0 * p44 * p43;

            // p₄₆ = 1 - 2 p₄₃²
            let p46 = 1.0 - 2.0 * p43.powi(2);

            // p₄₇ = (¹/₂ J₂ / pₗ) / pₗ
            let p47 = 0.5 * self.geopotential.j2 / pl / pl;

            // rₖ = r (1 - ³/₂ p₄₇ β p₃₆) + ¹/₂ (¹/₂ J₂ / pₗ) p₃₃ p₄₆
            let rk = r * (1.0 - 1.5 * p47 * b * p36)
                + 0.5 * (0.5 * self.geopotential.j2 / pl) * p33 * p46;

            // uₖ = u - ¹/₄ p₄₇ p₃₄ p₄₅
            let uk = u - 0.25 * p47 * p34 * p45;

            // Iₖ = I + ³/₂ p₄₇ cos I sin I p₄₆
            let inclination_k = orbit.inclination
                + 1.5 * p47 * orbit.inclination.cos() * orbit.inclination.sin() * p46;

            // Ωₖ = Ω + ³/₂ p₄₇ cos I p₄₅
            let right_ascension_k =
                orbit.right_ascension + 1.5 * p47 * orbit.inclination.cos() * p45;

            // ṙₖ = ṙ + n (¹/₂ J₂ / pₗ) p₃₃ / kₑ
            let rk_dot = r_dot
                - orbit.mean_motion * (0.5 * self.geopotential.j2 / pl) * p33 * p45
                    / self.geopotential.ke;

            // rḟₖ = pₗ¹ᐟ² / r + n (¹/₂ J₂ / pₗ) (p₃₃ p₄₆ + ³/₂ p₃₆) / kₑ
            let rfk_dot = pl.sqrt() / r
                + orbit.mean_motion * (0.5 * self.geopotential.j2 / pl) * (p33 * p46 + 1.5 * p36)
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

    pub fn propagate(&self, t: f64) -> Result<Prediction> {
        self.propagate_from_state(t, self.initial_state().as_mut(), false)
    }

    pub fn propagate_afspc_compatibility_mode(&self, t: f64) -> Result<Prediction> {
        self.propagate_from_state(t, self.initial_state().as_mut(), true)
    }
}
