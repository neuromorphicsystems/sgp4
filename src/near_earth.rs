use crate::gp;
use crate::model;
use crate::propagator;

#[allow(clippy::too_many_arguments)]
pub(crate) fn constants(
    geopotential: &model::Geopotential,
    drag_term: f64,
    orbit_0: propagator::Orbit,
    p1: f64,
    a0: f64,
    s: f64,
    xi: f64,
    eta: f64,
    c1: f64,
    c4: f64,
    k0: f64,
    k1: f64,
    k6: f64,
    k14: f64,
    p2: f64,
    p3: f64,
    p7: f64,
    p9: f64,
    p14: f64,
    p15: f64,
) -> propagator::Constants {
    propagator::Constants {
        geopotential,

        // Ω̇ = p₁₄
        right_ascension_dot: p14,

        // ω̇ = k₁₄
        argument_of_perigee_dot: k14,

        // Ṁ = p₁₅
        mean_anomaly_dot: p15,
        c1,
        c4,
        k0,
        k1,
        method: propagator::Method::NearEarth {
            a0,

            //        1 J₃
            // k₂ = - - -- sin I₀
            //        2 J₂
            k2: -0.5 * (geopotential.j3 / geopotential.j2) * orbit_0.inclination.sin(),

            // k₃ = 1 - p₁²
            k3: 1.0 - p1.powi(2),

            // k₄ = 7 p₁² - 1
            k4: 7.0 * p1.powi(2) - 1.0,

            //      │   1 J₃        3 + 5 p₁
            // k₅ = │ - - -- sin I₀ --------    if |1 + p₁| > 1.5 × 10⁻¹²
            //      │   4 J₂         1 + p₁
            //      │   1 J₃         3 + 5 p₁
            //      │ - - -- sin I₀ ----------- otherwise
            //      │   4 J₂        1.5 × 10⁻¹²
            k5: if (1.0 + p1).abs() > 1.5e-12 {
                -0.25
                    * (geopotential.j3 / geopotential.j2)
                    * orbit_0.inclination.sin()
                    * (3.0 + 5.0 * p1)
                    / (1.0 + p1)
            } else {
                -0.25
                    * (geopotential.j3 / geopotential.j2)
                    * orbit_0.inclination.sin()
                    * (3.0 + 5.0 * p1)
                    / 1.5e-12
            },
            k6,

            // p₃ < 220 / (aₑ + 1)
            high_altitude: if p3 < 220.0 / geopotential.ae + 1.0 {
                propagator::HighAltitude::No {}
            } else {
                // D₂ = 4 a₀" ξ C₁²
                let d2 = 4.0 * a0 * xi * c1.powi(2);

                // p₁₆ = D₂ ξ C₁ / 3
                let p16 = d2 * xi * c1 / 3.0;

                // D₃ = (17 a + s) p₁₆
                let d3 = (17.0 * a0 + s) * p16;

                // D₄ = ¹/₂ p₁₆ a₀" ξ (221 a₀" + 31 s) C₁
                let d4 = 0.5 * p16 * a0 * xi * (221.0 * a0 + 31.0 * s) * c1;

                propagator::HighAltitude::Yes {
                    // C₅ = 2 B* p₉ a₀" p₂ (1 + 2.75 (η² + η e₀) + e₀ η³)
                    c5: drag_term
                        * (2.0
                            * p9
                            * a0
                            * p2
                            * (1.0
                                + 2.75 * (eta.powi(2) + eta * orbit_0.eccentricity)
                                + eta * orbit_0.eccentricity * eta.powi(2))),
                    d2,
                    d3,
                    d4,
                    eta,

                    // k₇ = sin M₀
                    k7: orbit_0.mean_anomaly.sin(),

                    // k₈ = D₂ + 2 C₁²
                    k8: d2 + 2.0 * c1.powi(2),

                    // k₉ = ¹/₄ (3 D₃ + C₁ (12 D₂ + 10 C₁²))
                    k9: 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1.powi(2))),

                    // k₁₀ = ¹/₅ (3 D₄ + 12 C₁ D₃ + 6 D₂² + 15 C₁² (2 D₂ + C₁²))
                    k10: 0.2
                        * (3.0 * d4
                            + 12.0 * c1 * d3
                            + 6.0 * d2.powi(2)
                            + 15.0 * c1.powi(2) * (2.0 * d2 + c1.powi(2))),

                    elliptic: if orbit_0.eccentricity > 1.0e-4 {
                        propagator::Elliptic::Yes {
                            // k₁₁ = (1 + η cos M₀)³
                            k11: (1.0 + eta * orbit_0.mean_anomaly.cos()).powi(3),

                            //                     J₃ p₇ ξ  n₀" sin I₀
                            // k₁₂ = - 2 B* cos ω₀ -- ----------------
                            //                     J₂        e₀
                            k12: drag_term
                                * (-2.0
                                    * p7
                                    * xi
                                    * (geopotential.j3 / geopotential.j2)
                                    * orbit_0.mean_motion
                                    * orbit_0.inclination.sin()
                                    / orbit_0.eccentricity)
                                * orbit_0.argument_of_perigee.cos(),

                            //         2 p₇ B*
                            // k₁₃ = - - -----
                            //         3 e₀ η
                            k13: -2.0 / 3.0 * p7 * drag_term / (orbit_0.eccentricity * eta),
                        }
                    } else {
                        propagator::Elliptic::No {}
                    },
                }
            },
        },
        orbit_0,
    }
}

impl<'a> propagator::Constants<'a> {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn near_earth_orbital_elements(
        &self,
        a0: f64,
        k2: f64,
        k3: f64,
        k4: f64,
        k5: f64,
        k6: f64,
        high_altitude: &propagator::HighAltitude,
        t: f64,
        p22: f64,
        p23: f64,
    ) -> gp::Result<(propagator::Orbit, f64, f64, f64, f64, f64, f64)> {
        // p₂₄ = M₀ + Ṁ t
        let p24 = self.orbit_0.mean_anomaly + self.mean_anomaly_dot * t;
        let (argument_of_perigee, mean_anomaly, a, p27) = match high_altitude {
            propagator::HighAltitude::No {} => (
                // ω = p₂₃
                p23,
                // M = p₂₄ + n₀" k₁ t²
                p24 + self.orbit_0.mean_motion * self.k1 * t.powi(2),
                // a = a₀" (1 - C₁ t)²
                a0 * (1.0 - self.c1 * t).powi(2),
                // p₂₇ = e₀ - C₄ t
                self.orbit_0.eccentricity - self.c4 * t,
            ),
            propagator::HighAltitude::Yes {
                c5,
                d2,
                d3,
                d4,
                eta,
                k7,
                k8,
                k9,
                k10,
                elliptic,
            } => {
                // ω = │ p₂₃ - p₂₅ if e₀ > 10⁻⁴
                //     │ p₂₃       otherwise
                // p₂₆ = │ p₂₄ + p₂₅ if e₀ > 10⁻⁴
                //       │ p₂₄       otherwise
                let (argument_of_perigee, p26) = match elliptic {
                    propagator::Elliptic::Yes { k11, k12, k13 } => {
                        // p₂₅ = k₁₃ ((1 + η cos p₂₄)³ - k₁₁) + k₁₂ t
                        let p25 = k13 * ((1.0 + eta * p24.cos()).powi(3) - k11) + k12 * t;
                        (p23 - p25, p24 + p25)
                    }
                    propagator::Elliptic::No {} => (p23, p24),
                };
                (
                    argument_of_perigee,
                    // M = p₂₆ + n₀" (k₁ t² + k₈ t³ + t⁴ (k₉ + t k₁₀)
                    p26 + self.orbit_0.mean_motion
                        * (self.k1 * t.powi(2) + k8 * t.powi(3) + t.powi(4) * (k9 + t * k10)),
                    // a = a₀" (1 - C₁ t - D₂ t² - D₃ t³ - D₄ t⁴)²
                    a0 * (1.0 - self.c1 * t - d2 * t.powi(2) - d3 * t.powi(3) - d4 * t.powi(4))
                        .powi(2),
                    // p₂₇ = e₀ - (C₄ t + C₅ (sin p₂₆ - k₇))
                    self.orbit_0.eccentricity - (self.c4 * t + c5 * (p26.sin() - k7)),
                )
            }
        };
        if !(-0.001..1.0).contains(&p27) {
            Err(gp::Error::new("diverging eccentricity".to_owned()))
        } else {
            // e = │ 10⁻⁶ if p₂₇ < 10⁻⁶
            //     │ p₂₇  otherwise
            let eccentricity = p27.max(1.0e-6);
            Ok((
                propagator::Orbit {
                    // I = I₀
                    inclination: self.orbit_0.inclination,

                    // Ω = p₂₂
                    right_ascension: p22,
                    eccentricity,
                    argument_of_perigee,
                    mean_anomaly,

                    // n = kₑ / a³ᐟ²
                    mean_motion: self.geopotential.ke / a.powf(1.5),
                },
                a,
                // p₃₂ = k₂
                k2,
                // p₃₃ = k₃
                k3,
                // p₃₄ = k₄
                k4,
                // p₃₅ = k₅
                k5,
                // p₃₆ = k₆
                k6,
            ))
        }
    }
}
