fn main() -> sgp4::Result<()> {
    let elements: sgp4::Elements = serde_json::from_str(
        r#"{
            "OBJECT_NAME": "ISS (ZARYA)",
            "OBJECT_ID": "1998-067A",
            "EPOCH": "2020-07-12T21:16:01.000416",
            "MEAN_MOTION": 15.49507896,
            "ECCENTRICITY": 0.0001413,
            "INCLINATION": 51.6461,
            "RA_OF_ASC_NODE": 221.2784,
            "ARG_OF_PERICENTER": 89.1723,
            "MEAN_ANOMALY": 280.4612,
            "EPHEMERIS_TYPE": 0,
            "CLASSIFICATION_TYPE": "U",
            "NORAD_CAT_ID": 25544,
            "ELEMENT_SET_NO": 999,
            "REV_AT_EPOCH": 23600,
            "BSTAR": -3.1515e-5,
            "MEAN_MOTION_DOT": -2.218e-5,
            "MEAN_MOTION_DDOT": 0
        }"#,
    )?;
    let constants = sgp4::Constants::from_elements(&elements)?;
    for hours in 0..24 {
        println!("t = {} min", hours * 60);
        let prediction = constants.propagate((hours * 60) as f64)?;
        println!("    r = {:?} km", prediction.position);
        println!("    ṙ = {:?} km.s⁻¹", prediction.velocity);
    }
    Ok(())
}
