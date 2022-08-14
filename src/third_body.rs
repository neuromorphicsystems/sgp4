use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub(crate) struct Perturbations {
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

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields, rename_all = "camelCase")]
pub(crate) struct Dots {
    pub(crate) inclination: f64,
    pub(crate) right_ascension: f64,
    pub(crate) eccentricity: f64,
    pub(crate) argument_of_perigee: f64,
    pub(crate) mean_anomaly: f64,
}

// inclination_0: the angle between the equator and the orbit plane i₀, in rad
// right_ascension: the angle between vernal equinox and the point where
//                  the orbit crosses the equatorial plane Ω₀, in rad
// eccentricity_0: the shape of the orbit e₀
// argument_of_perigee: the angle between the ascending node and the orbit's
//                      point of closest approach to the earth ω₀, in rad
// n0: mean number of orbits per day (Brouwer mean motion) n₀", in rad.min⁻¹
// third_body_inclination_sine: the sine of the third body's inclination_0 sin Iₓ
// third_body_inclination_cosine: the cosine of the third body's inclination_0 cos Iₓ
// delta_right_ascension_sine: the sine of the third body's relative
//                             right ascension of the ascending node sin(Ω₀ - Ωₓ)
// delta_right_ascension_cosine: the cosine of the third body's relative
//                               right ascension of the ascending node cos(Ω₀ - Ωₓ)
// eccentricity_0: the third body's eccentricity_0 eₓ
// third_body_argument_of_perigee_sine: the sine of the third body's argument of perigee sin ωₓ
// third_body_argument_of_perigee_cosine: the cosine of the third body's argument of perigee cos ωₓ
// third_body_perturbation_coefficient: linear scaling factor Cₓ, in rad.min⁻¹
// third_body_mean_motion: the third body's mean motion nₓ
// p1: the constant p₁ = 1 − e₀²
// b0: the constant β₀ = p₁¹ᐟ²
#[allow(clippy::too_many_arguments)]
pub(crate) fn perturbations_and_dots(
    inclination_0: f64,
    eccentricity_0: f64,
    argument_of_perigee_0: f64,
    n0: f64,
    third_body_inclination_sine: f64,
    third_body_inclination_cosine: f64,
    delta_right_ascension_sine: f64,
    delta_right_ascension_cosine: f64,
    third_body_eccentricity: f64,
    third_body_argument_of_perigee_sine: f64,
    third_body_argument_of_perigee_cosine: f64,
    third_body_perturbation_coefficient: f64,
    third_body_mean_motion: f64,
    third_body_mean_anomaly_0: f64,
    p1: f64,
    b0: f64,
) -> (Perturbations, Dots) {
    // aₓ₁ = cos ωₓ cos(Ω₀ - Ωₓ) + sin ωₓ cos Iₓ sin(Ω₀ - Ωₓ)
    let ax1 = third_body_argument_of_perigee_cosine * delta_right_ascension_cosine
        + third_body_argument_of_perigee_sine
            * third_body_inclination_cosine
            * delta_right_ascension_sine;

    // aₓ₃ = - sin ωₓ cos(Ω₀ - Ωₓ) + cos ωₓ cos Iₓ sin(Ω₀ - Ωₓ)
    let ax3 = -third_body_argument_of_perigee_sine * delta_right_ascension_cosine
        + third_body_argument_of_perigee_cosine
            * third_body_inclination_cosine
            * delta_right_ascension_sine;

    // aₓ₇ = - cos ωₓ sin(Ω₀ - Ωₓ) + sin ωₓ cos Iₓ cos(Ω₀ - Ωₓ)
    let ax7 = -third_body_argument_of_perigee_cosine * delta_right_ascension_sine
        + third_body_argument_of_perigee_sine
            * third_body_inclination_cosine
            * delta_right_ascension_cosine;

    // aₓ₈ = sin ωₓ sin Iₓ
    let ax8 = third_body_argument_of_perigee_sine * third_body_inclination_sine;

    // aₓ₉ = sin ωₓ sin(Ω₀ - Ωₓ) + cos ωₓ cos Iₓ cos(Ω₀ - Ωₓ)
    let ax9 = third_body_argument_of_perigee_sine * delta_right_ascension_sine
        + third_body_argument_of_perigee_cosine
            * third_body_inclination_cosine
            * delta_right_ascension_cosine;

    // aₓ₁₀ = cos ωₓ sin Iₓ
    let ax10 = third_body_argument_of_perigee_cosine * third_body_inclination_sine;

    // aₓ₂ = aₓ₇ cos I₀ + aₓ₈ sin I₀
    let ax2 = inclination_0.cos() * ax7 + inclination_0.sin() * ax8;

    // aₓ₄ = aₓ₉ cos I₀ + aₓ₁₀ sin I₀
    let ax4 = inclination_0.cos() * ax9 + inclination_0.sin() * ax10;

    // aₓ₅ = - aₓ₇ sin I₀ + aₓ₈ cos I₀
    let ax5 = -inclination_0.sin() * ax7 + inclination_0.cos() * ax8;

    // aₓ₆ = - aₓ₉ sin I₀ + aₓ₁₀ cos I₀
    let ax6 = -inclination_0.sin() * ax9 + inclination_0.cos() * ax10;

    // Xₓ₁ = aₓ₁ cos ω₀ + aₓ₂ sin ω₀
    let xx1 = ax1 * argument_of_perigee_0.cos() + ax2 * argument_of_perigee_0.sin();

    // Xₓ₂ = aₓ₃ cos ω₀ + aₓ₄ sin ω₀
    let xx2 = ax3 * argument_of_perigee_0.cos() + ax4 * argument_of_perigee_0.sin();

    // Xₓ₃ = - aₓ₁ sin ω₀ + aₓ₂ cos ω₀
    let xx3 = -ax1 * argument_of_perigee_0.sin() + ax2 * argument_of_perigee_0.cos();

    // Xₓ₄ = - aₓ₃ sin ω₀ + aₓ₄ cos ω₀
    let xx4 = -ax3 * argument_of_perigee_0.sin() + ax4 * argument_of_perigee_0.cos();

    // Xₓ₅ = aₓ₅ sin ω₀
    let xx5 = ax5 * argument_of_perigee_0.sin();

    // Xₓ₆ = aₓ₆ sin ω₀
    let xx6 = ax6 * argument_of_perigee_0.sin();

    // Xₓ₇ = aₓ₅ cos ω₀
    let xx7 = ax5 * argument_of_perigee_0.cos();

    // Xₓ₈ = aₓ₆ cos ω₀
    let xx8 = ax6 * argument_of_perigee_0.cos();

    // Zₓ₃₁ = 12 Xₓ₁² - 3 Xₓ₃²
    let zx31 = 12.0 * xx1.powi(2) - 3.0 * xx3.powi(2);

    // Zₓ₃₂ = 24 Xₓ₁ Xₓ₂ - 6 Xₓ₃ Xₓ₄
    let zx32 = 24.0 * xx1 * xx2 - 6.0 * xx3 * xx4;

    // Zₓ₃₃ = 12 Xₓ₂² - 3 Xₓ₄²
    let zx33 = 12.0 * xx2.powi(2) - 3.0 * xx4.powi(2);

    // Zₓ₁₁ = - 6 aₓ₁ aₓ₅ + e₀² (- 24 Xₓ₁ Xₓ₇ - 6 Xₓ₃ Xₓ₅)
    let zx11 = -6.0 * ax1 * ax5 + eccentricity_0.powi(2) * (-24.0 * xx1 * xx7 - 6.0 * xx3 * xx5);

    // Zₓ₁₃ = - 6 aₓ₃ aₓ₆ + e₀² (-24 Xₓ₂ Xₓ₈ - 6 Xₓ₄ Xₓ₆)
    let zx13 = -6.0 * ax3 * ax6 + eccentricity_0.powi(2) * (-24.0 * xx2 * xx8 - 6.0 * xx4 * xx6);

    // Zₓ₂₁ = 6 aₓ₂ aₓ₅ + e₀² (24.0 Xₓ₁ Xₓ₅ - 6 Xₓ₃ Xₓ₇)
    let zx21 = 6.0 * ax2 * ax5 + eccentricity_0.powi(2) * (24.0 * xx1 * xx5 - 6.0 * xx3 * xx7);

    // Zₓ₂₃ = 6 aₓ₄ aₓ₆ + e₀² (24 Xₓ₂ Xₓ₆ - 6 Xₓ₄ Xₓ₈)
    let zx23 = 6.0 * ax4 * ax6 + eccentricity_0.powi(2) * (24.0 * xx2 * xx6 - 6.0 * xx4 * xx8);

    // Zₓ₁ = 2 (3 (aₓ₁² + aₓ₂²) + Zₓ₃₁ e₀²) + p₁ Zₓ₃₁
    let zx1 = (3.0 * (ax1.powi(2) + ax2.powi(2)) + zx31 * eccentricity_0.powi(2)) * 2.0 + p1 * zx31;

    // Zₓ₃ = 2 (3 (aₓ₃² + aₓ₄²) + Zₓ₃₃ e₀²) + p₁ Zₓ₃₃
    let zx3 = (3.0 * (ax3.powi(2) + ax4.powi(2)) + zx33 * eccentricity_0.powi(2)) * 2.0 + p1 * zx33;

    // pₓ₀ = Cₓ / n₀"
    let px0 = third_body_perturbation_coefficient / n0;

    //         1 pₓ₀
    // pₓ₁ = - - ---
    //         2 β₀
    let px1 = -0.5 * px0 / b0;

    // pₓ₂ = pₓ₀ β₀
    let px2 = px0 * b0;

    // pₓ₃ = - 15 e₀ pₓ₂
    let px3 = -15.0 * eccentricity_0 * px2;

    // Ω̇ₓ = │ 0                               if I₀ < 5.2359877 × 10⁻²
    //      │                                 or I₀ > π - 5.2359877 × 10⁻²
    //      │ - nₓ pₓ₁ (Zₓ₂₁ + Zₓ₂₃) / sin I₀ otherwise
    let third_body_right_ascension_dot =
        if !(5.2359877e-2..=std::f64::consts::PI - 5.2359877e-2).contains(&inclination_0) {
            0.0
        } else {
            -third_body_mean_motion * px1 * (zx21 + zx23) / inclination_0.sin()
        };
    (
        Perturbations {
            // kₓ₀ = 2 pₓ₃ (Xₓ₂ Xₓ₃ + Xₓ₁ Xₓ₄)
            kx0: 2.0 * px3 * (xx2 * xx3 + xx1 * xx4),

            // kₓ₁ = 2 pₓ₃ (Xₓ₂ Xₓ₄ - Xₓ₁ Xₓ₃)
            kx1: 2.0 * px3 * (xx2 * xx4 - xx1 * xx3),

            // kₓ₂ = 2 pₓ₁ (- 6 (aₓ₁ aₓ₆ + aₓ₃ aₓ₅) + e₀² (- 24 (Xₓ₂ Xₓ₇ + Xₓ₁ Xₓ₈) - 6 (Xₓ₃ Xₓ₆ + Xₓ₄ Xₓ₅)))
            kx2: 2.0
                * px1
                * (-6.0 * (ax1 * ax6 + ax3 * ax5)
                    + eccentricity_0.powi(2)
                        * (-24.0 * (xx2 * xx7 + xx1 * xx8) - 6.0 * (xx3 * xx6 + xx4 * xx5))),

            // kₓ₃ = 2 pₓ₁ (Zₓ₁₃ - Zₓ₁₁)
            kx3: 2.0 * px1 * (zx13 - zx11),

            // kₓ₄ = - 2 pₓ₀ (2 (6 (aₓ₁ aₓ₃ + aₓ₂ aₓ₄) + Zₓ₃₂ e₀²) + p₁ Zₓ₃₂)
            kx4: -2.0
                * px0
                * ((6.0 * (ax1 * ax3 + ax2 * ax4) + zx32 * eccentricity_0.powi(2)) * 2.0
                    + p1 * zx32),

            // kₓ₅ = - 2 pₓ₀ (Zₓ₃ - Zₓ₁)
            kx5: -2.0 * px0 * (zx3 - zx1),

            // kₓ₆ = - 2 pₓ₀ (- 21 - 9 e₀²) eₓ
            kx6: -2.0 * px0 * (-21.0 - 9.0 * eccentricity_0.powi(2)) * third_body_eccentricity,

            // kₓ₇ = 2 pₓ₂ Zₓ₃₂
            kx7: 2.0 * px2 * zx32,

            // kₓ₈ = 2 pₓ₂ (Zₓ₃₃ - Zₓ₃₁)
            kx8: 2.0 * px2 * (zx33 - zx31),

            // kₓ₉ = - 18 pₓ₂ eₓ
            kx9: -18.0 * px2 * third_body_eccentricity,

            // kₓ₁₀ = - 2 pₓ₁ (6 (aₓ₄ aₓ₅ + aₓ₂ aₓ₆) + e₀² (24 (Xₓ₂ Xₓ₅ + Xₓ₁ Xₓ₆) - 6 (Xₓ₄ Xₓ₇ + Xₓ₃ Xₓ₈)))
            kx10: -2.0
                * px1
                * (6.0 * (ax4 * ax5 + ax2 * ax6)
                    + eccentricity_0.powi(2)
                        * (24.0 * (xx2 * xx5 + xx1 * xx6) - 6.0 * (xx4 * xx7 + xx3 * xx8))),

            // kₓ₁₁ = - 2 pₓ₁ (Zₓ₂₃ - Zₓ₂₁)
            kx11: -2.0 * px1 * (zx23 - zx21),
            third_body_mean_anomaly_0,
        },
        Dots {
            // İₓ = pₓ₁ nₓ (Zₓ₁₁ + Zₓ₁₃)
            inclination: px1 * third_body_mean_motion * (zx11 + zx13),
            right_ascension: third_body_right_ascension_dot,

            // ėₓ = pₓ₃ nₓ (Xₓ₁ Xₓ₃ + Xₓ₂ Xₓ₄)
            eccentricity: px3 * third_body_mean_motion * (xx1 * xx3 + xx2 * xx4),

            // ω̇ₓ = pₓ₂ nₓ (Zₓ₃₁ + Zₓ₃₃ - 6) - cos I₀ Ω̇ₓ
            argument_of_perigee: px2 * third_body_mean_motion * (zx31 + zx33 - 6.0)
                - inclination_0.cos() * third_body_right_ascension_dot,

            // Ṁₓ = - nₓ pₓ₀ (Zₓ₁ + Zₓ₃ - 14 - 6 e₀²)
            mean_anomaly: -third_body_mean_motion
                * px0
                * (zx1 + zx3 - 14.0 - 6.0 * eccentricity_0.powi(2)),
        },
    )
}

impl Perturbations {
    pub(crate) fn long_period_periodic_effects(
        &self,
        third_body_eccentricity: f64,
        third_body_mean_motion: f64,
        t: f64,
    ) -> (f64, f64, f64, f64, f64) {
        // Mₓ = Mₓ₀ + nₓ t
        let third_body_mean_anomaly = self.third_body_mean_anomaly_0 + third_body_mean_motion * t;

        // fₓ = Mₓ + 2 eₓ sin Mₓ
        let fx =
            third_body_mean_anomaly + 2.0 * third_body_eccentricity * third_body_mean_anomaly.sin();

        // Fₓ₂ = ¹/₂ sin²fₓ - ¹/₄
        let fx2 = 0.5 * fx.sin().powi(2) - 0.25;

        // Fₓ₃ = - ¹/₂ sin fₓ cos fₓ
        let fx3 = -0.5 * fx.sin() * fx.cos();
        (
            // δeₓ = kₓ₀ Fₓ₂ + kₓ₁ Fₓ₃
            self.kx0 * fx2 + self.kx1 * fx3,
            // δIₓ = kₓ₂ Fₓ₂ + kₓ₃ Fₓ₃
            self.kx2 * fx2 + self.kx3 * fx3,
            // δMₓ = kₓ₄ Fₓ₂ + kₓ₅ Fₓ₃ + kₓ₆ sin fₓ
            self.kx4 * fx2 + self.kx5 * fx3 + self.kx6 * fx.sin(),
            // pₓ₄ = kₓ₇ Fₓ₂ + kₓ₈ Fₓ₃ + kₓ₉ sin fₓ
            self.kx7 * fx2 + self.kx8 * fx3 + self.kx9 * fx.sin(),
            // pₓ₅ = kₓ₁₀ Fₓ₂ + kₓ₁₁ Fₓ₃
            self.kx10 * fx2 + self.kx11 * fx3,
        )
    }
}
