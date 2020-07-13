fn main() -> sgp4::Result<()> {
    let elements = sgp4::Elements::from_tle(
        Some("ISS (ZARYA)".to_owned()),
        "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992".as_bytes(),
        "2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008".as_bytes(),
    )?;
    let constants = sgp4::Constants::from_elements_afspc_compatibility_mode(&elements)?;
    for hours in 0..24 {
        println!("t = {} min", hours * 60);
        let prediction = constants.propagate_afspc_compatibility_mode((hours * 60) as f64)?;
        println!("    r = {:?} km", prediction.position);
        println!("    ṙ = {:?} km.s⁻¹", prediction.velocity);
    }
    Ok(())
}
