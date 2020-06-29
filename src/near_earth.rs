use crate::model;
use crate::propagator;

pub fn constants<'a>(
    geopotential: &'a model::Geopotential,
    drag_term: f64,
    orbit_0: propagator::Orbit,
    p0: f64,
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
    p1: f64,
    p3: f64,
    p6: f64,
    p8: f64,
    p13: f64,
    p14: f64,
) -> propagator::Constants<'a> {
    propagator::Constants {
        geopotential: geopotential,
        drag_term: drag_term,

        // Ω̇ = p₁₃
        right_ascension_dot: p13,

        // ω̇ = k₁₄
        argument_of_perigee_dot: k14,

        // Ṁ = p₁₄
        mean_anomaly_dot: p14,
        c1: c1,
        c4: c4,
        k0: k0,
        k1: k1,
        method: propagator::Method::NearEarth {
            a0: a0,

            //        1 J₃
            // k₂ = - - -- sin I₀
            //        2 J₂
            k2: -0.5 * (geopotential.j3 / geopotential.j2) * orbit_0.inclination.sin(),

            // k₃ = 1 - p₀²
            k3: 1.0 - p0.powi(2),

            // k₄ = 7 p₀² - 1
            k4: 7.0 * p0.powi(2) - 1.0,

            //      │   1 J₃        3 + 5 p₀
            // k₅ = │ - - -- sin I₀ --------    if |1 + p₀| > 1.5 × 10⁻¹²
            //      │   4 J₂         1 + p₀
            //      │   1 J₃         3 + 5 p₀
            //      │ - - -- sin I₀ ----------- otherwise
            //      │   4 J₂        1.5 × 10⁻¹²
            k5: if (1.0 + p0).abs() > 1.5e-12 {
                -0.25
                    * (geopotential.j3 / geopotential.j2)
                    * orbit_0.inclination.sin()
                    * (3.0 + 5.0 * p0)
                    / (1.0 + p0)
            } else {
                -0.25
                    * (geopotential.j3 / geopotential.j2)
                    * orbit_0.inclination.sin()
                    * (3.0 + 5.0 * p0)
                    / 1.5e-12
            },
            k6: k6,
            full: if p3 < 220.0 / geopotential.ae + 1.0 {
                propagator::Full::No {}
            } else {
                // D₂ = 4 a₀" ξ C₁²
                let d2 = 4.0 * a0 * xi * c1.powi(2);

                // p₁₅ = D₂ ξ C₁ / 3
                let p15 = d2 * xi * c1 / 3.0;

                // D₃ = (17 a + s) p₁₅
                let d3 = (17.0 * a0 + s) * p15;

                // D₄ = 0.5 p₁₅ a₀" ξ (221 a₀" + 31 s) C₁;
                let d4 = 0.5 * p15 * a0 * xi * (221.0 * a0 + 31.0 * s) * c1;

                propagator::Full::Yes {
                    // C₅ = 2 p₈ a₀" p₁ (1 + 2.75 (η² + η e₀) + e₀ η³)
                    c5: 2.0
                        * p8
                        * a0
                        * p1
                        * (1.0
                            + 2.75 * (eta.powi(2) + eta * orbit_0.eccentricity)
                            + eta * orbit_0.eccentricity * eta.powi(2)),
                    d2: d2,
                    d3: d3,
                    d4: d4,
                    eta: eta,

                    // k₇ = (1 + η cos M₀)³
                    k7: (1.0 + eta * orbit_0.mean_anomaly.cos()).powi(3),

                    // k₈ = sin M₀
                    k8: orbit_0.mean_anomaly.sin(),

                    // k₉ = D₂ + 2 C₁²
                    k9: d2 + 2.0 * c1.powi(2),

                    // k₁₀ = ¹/₄ (3 D₃ + C₁ (12 D₂ + 10 C₁²))
                    k10: 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1.powi(2))),

                    // k₁₁ = ¹/₅ (3 D₄ + 12 C₁ D₃ + 6 D₂² + 15 C₁² (2 D₂ + C₁²))
                    k11: 0.2
                        * (3.0 * d4
                            + 12.0 * c1 * d3
                            + 6.0 * d2.powi(2)
                            + 15.0 * c1.powi(2) * (2.0 * d2 + c1.powi(2))),

                    elliptic: if orbit_0.eccentricity > 1.0e-4 {
                        propagator::Elliptic::Yes {
                            //                    J₃ p₆ ξ  n₀" sin I₀
                            // k₁₂ = - 2 B* cos ω₀ -- ----------------
                            //                    J₂        e₀
                            k12: drag_term
                                * (-2.0
                                    * p6
                                    * xi
                                    * (geopotential.j3 / geopotential.j2)
                                    * orbit_0.mean_motion
                                    * orbit_0.inclination.sin()
                                    / orbit_0.eccentricity)
                                * orbit_0.argument_of_perigee.cos(),

                            //         2 p₆ B*
                            // k₁₃ = - - -----
                            //         3 e₀ η
                            k13: -2.0 / 3.0 * p6 * drag_term / (orbit_0.eccentricity * eta),
                        }
                    } else {
                        propagator::Elliptic::No {}
                    },
                }
            },
        },
        orbit_0: orbit_0,
    }
}

impl<'a> propagator::Constants<'a> {
    pub fn near_earth_orbital_elements(
        &self,
        a0: f64,
        k2: f64,
        k3: f64,
        k4: f64,
        k5: f64,
        k6: f64,
        full: &propagator::Full,
        t: f64,
        p21: f64,
        p22: f64,
    ) -> propagator::Result<(propagator::Orbit, f64, f64, f64, f64, f64, f64, f64)> {
        // p₂₃ = M₀ + Ṁ t
        let p23 = self.orbit_0.mean_anomaly + self.mean_anomaly_dot * t;
        let (argument_of_perigee, mean_anomaly, a, l, p25) = match full {
            propagator::Full::No {} => (
                // ω = p₂₂
                p22,
                // M = p₂₃
                p23,
                // a = a₀" (1 - C₁ t)²
                a0 * (1.0 - self.c1 * t).powi(2),
                // 𝕃 = p₂₃ + n₀" k₁ t²
                p23 + self.orbit_0.mean_motion * self.k1 * t.powi(2),
                // p₂₅ = e₀ - B* C₄ t
                self.orbit_0.eccentricity - self.drag_term * self.c4 * t,
            ),
            propagator::Full::Yes {
                c5,
                d2,
                d3,
                d4,
                eta,
                k7,
                k8,
                k9,
                k10,
                k11,
                elliptic,
            } => {
                // ω = │ p₂₂ - p₂₄ if e₀ > 10⁻⁴
                //     │ p₂₂       otherwise
                // M = │ p₂₃ + p₂₄ if e₀ > 10⁻⁴
                //     │ p₂₃       otherwise
                let (argument_of_perigee, mean_anomaly) = match elliptic {
                    propagator::Elliptic::Yes { k12, k13 } => {
                        // p₂₄ = k₁₃ ((1 + η cos p₂₃)³ - k₇) + k₁₂ t
                        let p24 = k13 * ((1.0 + eta * p23.cos()).powi(3) - k7) + k12 * t;
                        (p22 - p24, p23 + p24)
                    }
                    propagator::Elliptic::No {} => (p22, p23),
                };
                (
                    argument_of_perigee,
                    mean_anomaly,
                    // a = a₀" (1 - C₁ t - D₂ t² - D₃ t³ - D₄ t⁴)²
                    a0 * (1.0 - self.c1 * t - d2 * t.powi(2) - d3 * t.powi(3) - d4 * t.powi(4))
                        .powi(2),
                    // 𝕃 = M + n₀" (k₁ t² + k₉ t³ + t⁴ (k₁₀ + t k₁₁)
                    mean_anomaly
                        + self.orbit_0.mean_motion
                            * (self.k1 * t.powi(2) + k9 * t.powi(3) + t.powi(4) * (k10 + t * k11)),
                    // p₂₅ = e₀ - (B* C₄ t + B* C₅ (sin M - k₈))
                    self.orbit_0.eccentricity
                        - (self.drag_term * self.c4 * t
                            + self.drag_term * c5 * (mean_anomaly.sin() - k8)),
                )
            }
        };
        if p25 >= 1.0 || p25 < -0.001 {
            Err(propagator::Error::new("diverging eccentricity"))
        } else {
            // e = │ 10⁻⁶ if p₂₅ < 10⁻⁶
            //     │ p₂₅  otherwise
            let eccentricity = p25.max(1.0e-6);
            Ok((
                propagator::Orbit {
                    // I = I₀
                    inclination: self.orbit_0.inclination,

                    // Ω = p₂₁
                    right_ascension: p21,
                    eccentricity: eccentricity,
                    argument_of_perigee: argument_of_perigee,
                    mean_anomaly: mean_anomaly,

                    // n = kₑ / a³ᐟ²
                    mean_motion: self.geopotential.ke / a.powf(1.5),
                },
                a,
                l,
                // p₃₀ = k₂
                k2,
                // p₃₁ = k₃
                k3,
                // p₃₂ = k₄
                k4,
                // p₃₃ = k₅
                k5,
                // p₃₄ = k₆
                k6,
            ))
        }
    }
}
