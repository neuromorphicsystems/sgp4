fn main() -> sgp4::Result<()> {
    let elements = sgp4::Elements::from_tle(
        Some("MOLNIYA 1-36".to_owned()),
        "1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813".as_bytes(),
        "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656".as_bytes(),
    )?;
    let constants = sgp4::Constants::from_elements(&elements)?;
    let mut state = constants.initial_state();
    for days in 0..7 {
        println!("t = {} min", days * 60 * 24);
        let prediction =
            constants.propagate_from_state((days * 60 * 24) as f64, state.as_mut(), false)?;
        println!("    r = {:?} km", prediction.position);
        println!("    ṙ = {:?} km.s⁻¹", prediction.velocity);
    }
    Ok(())
}
