/// Verify that gsto (Greenwich Sidereal Time at epoch) matches the
/// Vallado C++ reference implementation.
///
/// Before this fix, `from_elements_afspc_compatibility_mode` used
/// `afspc_epoch_to_sidereal_time` for gsto. However, the Vallado C++
/// always uses `gstime_SGP4` (the IAU polynomial) — the AFSPC formula
/// is computed in initl but unconditionally overwritten (SGP4.cpp:1273):
///
///     gsto = gstime_SGP4(epoch + 2433281.5);
///
/// At t=0 the gsto difference cancels for resonant deep-space orbits
/// (subtracted into lambda_0 at init, added back at propagation).
/// These tests use long propagation times where the difference
/// accumulates and becomes measurable.

fn check(line1: &str, line2: &str, tsince: f64, ref_pos: [f64; 3], tol_km: f64, label: &str) {
    let tle = format!("{line1}\n{line2}\n");
    let elements = sgp4::parse_2les(&tle).unwrap();
    let constants =
        sgp4::Constants::from_elements_afspc_compatibility_mode(&elements[0]).unwrap();
    let pred = constants
        .propagate_afspc_compatibility_mode(sgp4::MinutesSinceEpoch(tsince))
        .unwrap();
    for i in 0..3 {
        let delta = (pred.position[i] - ref_pos[i]).abs();
        assert!(
            delta < tol_km,
            "{label} t={tsince} component {}: delta {delta:.3e} km exceeds {tol_km:.0e} km\n  \
             got:      {:.15}\n  expected: {:.15}",
            ["x", "y", "z"][i],
            pred.position[i],
            ref_pos[i]
        );
    }
}

#[test]
fn gsto_geo_resonant_long_propagation() {
    // GEO-resonant satellite (29238) propagated ~1 year (500,000 min).
    // The ~1e-10 rad gsto error accumulates at GEO altitude to ~3.7 mm
    // with the old AFSPC formula, exceeding the 1 mm tolerance.
    //
    // Reference: Vallado SGP4.cpp (v2020-07-13), clang -O3, WGS72, opsmode 'a'.
    check(
        "1 29238U 06022G   06177.28732010  .00000104  00000-0  10000-3 0   886",
        "2 29238   0.0536 356.7318 0004925 230.4823 326.2002  1.00271289 11990",
        500000.0,
        [33275.035513542701665, -25842.368923243404424, -7.464642954329388],
        1.0e-6, // 1 mm
        "sat 29238 (GEO resonant)",
    );
}