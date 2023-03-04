fn main() -> anyhow::Result<()> {
    let response = ureq::get("https://celestrak.com/NORAD/elements/gp.php")
        .query("GROUP", "galileo")
        .query("FORMAT", "json")
        .call()?;
    let elements_vec: Vec<sgp4::Elements> = response.into_json()?;
    for elements in &elements_vec {
        println!("{}", elements.object_name.as_ref().unwrap());
        let constants = sgp4::Constants::from_elements(elements)?;
        for hours in &[12, 24] {
            println!("    t = {} min", hours * 60);
            let prediction = constants.propagate((hours * 60) as f64)?;
            println!("        r = {:?} km", prediction.position);
            println!("        ṙ = {:?} km.s⁻¹", prediction.velocity);
        }
    }
    Ok(())
}
