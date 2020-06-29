use crate::model;

pub struct Perturbations {
    kx0: f64,
    kx1: f64,
    kx2: f64,
    kx3: f64,
    kx4: f64,
    kx5: f64,
    kx6: f64,
    kx7: f64,
    kx8: f64,
    kx9: f64,
    kx10: f64,
    kx11: f64,
    third_body_mean_anomaly_0: f64,
}

pub struct Dots {
    pub inclination: f64,
    pub right_ascension: f64,
    pub eccentricity: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
}

// inclination_0: the angle between the equator and the orbit plane i₀, in rad
// right_ascension: the angle between vernal equinox and the point where
//                  the orbit crosses the equatorial plane Ω₀, in rad
// eccentricity_0: the shape of the orbit e₀
// argument_of_perigee: the angle between the ascending node and the orbit's
//                      point of closest approach to the earth ω₀, in rad
// n0: mean number of orbits per day (Brouwer mean motion) n₀", in rad.min⁻¹
// third_body_inclination_0_sine: the sine of the third body's inclination_0 sin Iₓ
// third_body_inclination_0_cosine: the cosine of the third body's inclination_0 cos Iₓ
// delta_right_ascension_sine: the sine of the third body's relative
//                             right ascension of the ascending node sin(Ω₀ - Ωₓ)
// delta_right_ascension_cosine: the cosine of the third body's relative
//                               right ascension of the ascending node cos(Ω₀ - Ωₓ)
// eccentricity_0: the third body's eccentricity_0 eₓ
// third_body_argument_of_perigee_sine: the sine of the third body's argument of perigee sin ωₓ
// third_body_argument_of_perigee_cosine: the cosine of the third body's argument of perigee cos ωₓ
// third_body_perturbation_coefficient: linear scaling factor Cₓ, in rad.min⁻¹
// third_body_mean_motion: the third body's mean motion nₓ
// p0: the constant p₀ = 1 − e₀²
// b0: the constant β₀ = p₀¹ᐟ²
pub fn perturbations_and_dots(
    inclination_0: f64,
    eccentricity_0: f64,
    argument_of_perigee_0: f64,
    n0: f64,
    third_body_inclination_0_sine: f64,
    third_body_inclination_0_cosine: f64,
    delta_right_ascension_sine: f64,
    delta_right_ascension_cosine: f64,
    third_body_eccentricity_0: f64,
    third_body_argument_of_perigee_sine: f64,
    third_body_argument_of_perigee_cosine: f64,
    third_body_perturbation_coefficient: f64,
    third_body_mean_motion: f64,
    third_body_mean_anomaly_0: f64,
    p0: f64,
    b0: f64,
) -> (Perturbations, Dots) {
    // a₁ = cos ωₓ cos(Ω₀ - Ωₓ) + sin ωₓ cos Iₓ sin(Ω₀ - Ωₓ)
    let a1 = third_body_argument_of_perigee_cosine * delta_right_ascension_cosine
        + third_body_argument_of_perigee_sine
            * third_body_inclination_0_cosine
            * delta_right_ascension_sine;

    // a₃ = - sin ωₓ cos(Ω₀ - Ωₓ) + cos ωₓ cos Iₓ sin(Ω₀ - Ωₓ)
    let a3 = -third_body_argument_of_perigee_sine * delta_right_ascension_cosine
        + third_body_argument_of_perigee_cosine
            * third_body_inclination_0_cosine
            * delta_right_ascension_sine;

    // a₇ = - cos ωₓ sin(Ω₀ - Ωₓ) + sin ωₓ cos Iₓ cos(Ω₀ - Ωₓ)
    let a7 = -third_body_argument_of_perigee_cosine * delta_right_ascension_sine
        + third_body_argument_of_perigee_sine
            * third_body_inclination_0_cosine
            * delta_right_ascension_cosine;

    // a₈ = sin ωₓ sin Iₓ
    let a8 = third_body_argument_of_perigee_sine * third_body_inclination_0_sine;

    // a₉ = sin ωₓ sin(Ω₀ - Ωₓ) + cos ωₓ cos Iₓ cos(Ω₀ - Ωₓ)
    let a9 = third_body_argument_of_perigee_sine * delta_right_ascension_sine
        + third_body_argument_of_perigee_cosine
            * third_body_inclination_0_cosine
            * delta_right_ascension_cosine;

    // a₁₀ = cos ωₓ sin Iₓ
    let a10 = third_body_argument_of_perigee_cosine * third_body_inclination_0_sine;

    // a₂ = a₇ cos i₀ + a₈ sin i₀
    let a2 = inclination_0.cos() * a7 + inclination_0.sin() * a8;

    // a₄ = a₉ cos i₀ + a₁₀ sin i₀
    let a4 = inclination_0.cos() * a9 + inclination_0.sin() * a10;

    // a₅ = - a₇ sin i₀ + a₈ cos i₀
    let a5 = -inclination_0.sin() * a7 + inclination_0.cos() * a8;

    // a₆ = - a₉ sin i₀ + a₁₀ cos i₀
    let a6 = -inclination_0.sin() * a9 + inclination_0.cos() * a10;

    // X₁ = a₁ cos ω₀ + a₂ sin ω₀
    let x1 = a1 * argument_of_perigee_0.cos() + a2 * argument_of_perigee_0.sin();

    // X₂ = a₃ cos ω₀ + a₄ sin ω₀
    let x2 = a3 * argument_of_perigee_0.cos() + a4 * argument_of_perigee_0.sin();

    // X₃ = - a₁ sin ω₀ + a₂ cos ω₀
    let x3 = -a1 * argument_of_perigee_0.sin() + a2 * argument_of_perigee_0.cos();

    // X₄ = - a₃ sin ω₀ + a₄ cos ω₀
    let x4 = -a3 * argument_of_perigee_0.sin() + a4 * argument_of_perigee_0.cos();

    // X₅ = a₅ sin ω₀
    let x5 = a5 * argument_of_perigee_0.sin();

    // X₆ = a₆ sin ω₀
    let x6 = a6 * argument_of_perigee_0.sin();

    // X₇ = a₅ cos ω₀
    let x7 = a5 * argument_of_perigee_0.cos();

    // X₈ = a₆ cos ω₀
    let x8 = a6 * argument_of_perigee_0.cos();

    // Z₃₁ = 12 X₁² - 3 X₃²
    let z31 = 12.0 * x1.powi(2) - 3.0 * x3.powi(2);

    // Z₃₂ = 24 X₁ X₂ - 6 X₃ X₄
    let z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;

    // Z₃₃ = 12 X₂² - 3 X₄²
    let z33 = 12.0 * x2.powi(2) - 3.0 * x4.powi(2);

    // Z₁₁ = -6.0 a₁ a₅ + e₀² (-24 X₁ X₇ - 6 X₃ X₅)
    let z11 = -6.0 * a1 * a5 + eccentricity_0.powi(2) * (-24.0 * x1 * x7 - 6.0 * x3 * x5);

    // Z₁₂ = -6 (a₁ a₆ + a₃ a₅) + e₀² (-24 (X₂ X₇ + X₁ X₈) - 6 (X₃ X₆ + X₄ X₅))
    let z12 = -6.0 * (a1 * a6 + a3 * a5)
        + eccentricity_0.powi(2) * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));

    // Z₁₃ = -6 a₃ a₆ + e₀² (-24 X₂ X₈ - 6 X₄ X₆)
    let z13 = -6.0 * a3 * a6 + eccentricity_0.powi(2) * (-24.0 * x2 * x8 - 6.0 * x4 * x6);

    // Z₂₁ = 6 a₂ a₅ + e₀² (24.0 X₁ X₅ - 6 X₃ X₇)
    let z21 = 6.0 * a2 * a5 + eccentricity_0.powi(2) * (24.0 * x1 * x5 - 6.0 * x3 * x7);

    // Z₂₂ = 6 (a₄ a₅ + a₂ a₆) + e₀² (24 (X₂ X₅ + X₁ X₆) - 6 (X₄ X₇ + X₃ X₈))
    let z22 = 6.0 * (a4 * a5 + a2 * a6)
        + eccentricity_0.powi(2) * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));

    // Z₂₃ = 6 a₄ a₆ + e₀² (24 X₂ X₆ - 6 X₄ X₈)
    let z23 = 6.0 * a4 * a6 + eccentricity_0.powi(2) * (24.0 * x2 * x6 - 6.0 * x4 * x8);

    // Z₁ = 2 (3 (a₁² + a₂²) + Z₃₁ e₀²) + p₀ Z₃₁
    let z1 = (3.0 * (a1.powi(2) + a2.powi(2)) + z31 * eccentricity_0.powi(2)) * 2.0 + p0 * z31;

    // Z₂ = 2 (6 (a₁ a₃ + a₂ a₄) + Z₃₂ e₀²) + p₀ Z₃₂
    let z2 = (6.0 * (a1 * a3 + a2 * a4) + z32 * eccentricity_0.powi(2)) * 2.0 + p0 * z32;

    // Z₃ = 2 (3 (a₃² + a₄²) + Z₃₃ e₀²) + p₀ Z₃₃
    let z3 = (3.0 * (a3.powi(2) + a4.powi(2)) + z33 * eccentricity_0.powi(2)) * 2.0 + p0 * z33;

    // lₓ₀ = Cₓ / n₀"
    let lx0 = third_body_perturbation_coefficient / n0;

    //         1 lₓ₀
    // lₓ₁ = - - ---
    //         2 β₀
    let lx1 = -0.5 * lx0 / b0;

    // lₓ₂ = lₓ₀ β₀
    let lx2 = lx0 * b0;

    // lₓ₃ = -15 e₀ lₓ₂
    let lx3 = -15.0 * eccentricity_0 * lx2;

    // Ω̇ₓ = │ 0                             if i₀ < 5.2359877 × 10⁻²
    //      │                               or i₀ > π - 5.2359877 × 10⁻²
    //      │ - nₓ lₓ₁ (Z₂₁ + Z₂₃) / sin i₀ otherwise
    let third_body_right_ascension_dot =
        if inclination_0 < 5.2359877e-2 || inclination_0 > model::PI - 5.2359877e-2 {
            0.0
        } else {
            -third_body_mean_motion * lx1 * (z21 + z23) / inclination_0.sin()
        };
    (
        Perturbations {
            // kₓ₀ = 2 lₓ₃ (X₂ X₃ + X₁ X₄)
            kx0: 2.0 * lx3 * (x2 * x3 + x1 * x4),

            // kₓ₁ = 2 lₓ₃ (X₂ X₄ - X₁ X₃)
            kx1: 2.0 * lx3 * (x2 * x4 - x1 * x3),

            // kₓ₂ = 2 lₓ₁ Z₁₂
            kx2: 2.0 * lx1 * z12,

            // kₓ₃ = 2 lₓ₁ (Z₁₃ - Z₁₁)
            kx3: 2.0 * lx1 * (z13 - z11),

            // kₓ₄ = -2 lₓ₀ Z₂
            kx4: -2.0 * lx0 * z2,

            // kₓ₅ = - 2 lₓ₀ (Z₃ - Z₁)
            kx5: -2.0 * lx0 * (z3 - z1),

            // kₓ₆ = - 2 lₓ₀ (- 21 - 9 e₀²) eₓ
            kx6: -2.0 * lx0 * (-21.0 - 9.0 * eccentricity_0.powi(2)) * third_body_eccentricity_0,

            // kₓ₇ = 2 lₓ₂ Z₃₂
            kx7: 2.0 * lx2 * z32,

            // kₓ₈ = 2 lₓ₂ (Z₃₃ - Z₃₁)
            kx8: 2.0 * lx2 * (z33 - z31),

            // kₓ₉ = - 18 lₓ₂ eₓ
            kx9: -18.0 * lx2 * third_body_eccentricity_0,

            // kₓ₁₀ = - 2 lₓ₁ Z₂₂
            kx10: -2.0 * lx1 * z22,

            // kₓ₁₁ = - 2 lₓ₁ (Z₂₃ - Z₂₁)
            kx11: -2.0 * lx1 * (z23 - z21),
            third_body_mean_anomaly_0: third_body_mean_anomaly_0,
        },
        Dots {
            // İₓ = lₓ₁ nₓ (Z₁₁ + Z₁₃)
            inclination: lx1 * third_body_mean_motion * (z11 + z13),
            right_ascension: third_body_right_ascension_dot,

            // ėₓ = lₓ₃ nₓ (X₁ X₃ + X₂ X₄)
            eccentricity: lx3 * third_body_mean_motion * (x1 * x3 + x2 * x4),

            // ω̇ₓ = lₓ₂ nₓ (Z₃₁ + Z₃₃ - 6) - cos i₀ Ω̇ₓ
            argument_of_perigee: lx2 * third_body_mean_motion * (z31 + z33 - 6.0)
                - inclination_0.cos() * third_body_right_ascension_dot,

            // Ṁₓ = - nₓ lₓ₀ (Z₁ + Z₃ - 14 - 6 e₀²)
            mean_anomaly: -third_body_mean_motion
                * lx0
                * (z1 + z3 - 14.0 - 6.0 * eccentricity_0.powi(2)),
        },
    )
}

impl Perturbations {
    pub fn long_period_periodic_effects(
        &self,
        third_body_eccentricity_0: f64,
        third_body_mean_motion: f64,
        t: f64,
    ) -> (f64, f64, f64, f64, f64) {
        // Mₓ = Mₓ₀ + nₓ t
        let third_body_mean_anomaly = self.third_body_mean_anomaly_0 + third_body_mean_motion * t;

        // fₓ = Mₓ + 2 eₓ sin Mₓ
        let fx = third_body_mean_anomaly
            + 2.0 * third_body_eccentricity_0 * third_body_mean_anomaly.sin();

        // f₂ = ¹/₂ sin²fₓ - ¹/₄
        let f2 = 0.5 * fx.sin().powi(2) - 0.25;

        // f₃ = - ¹/₂ sin fₓ cos fₓ
        let f3 = -0.5 * fx.sin() * fx.cos();
        (
            // δeₓ = kₓ₀ f₂ + kₓ₁ f₃
            self.kx0 * f2 + self.kx1 * f3,
            // δIₓ = kₓ₂ f₂ + kₓ₃ f₃
            self.kx2 * f2 + self.kx3 * f3,
            // δMₓ = kₓ₄ f₂ + kₓ₅ f₃ + kₓ₆ sin fₓ
            self.kx4 * f2 + self.kx5 * f3 + self.kx6 * fx.sin(),
            // lₓ₄ = kₓ₇ f₂ + kₓ₈ f₃ + kₓ₉ sin fₓ
            self.kx7 * f2 + self.kx8 * f3 + self.kx9 * fx.sin(),
            // lₓ₅ = kₓ₁₀ f₂ + kₓ₁₁ f₃
            self.kx10 * f2 + self.kx11 * f3,
        )
    }
}
