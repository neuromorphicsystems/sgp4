fn main() -> sgp4::Result<()> {
    let response = ureq::post("https://www.space-track.org/ajaxauth/login").send_form(&[
        ("identity", "your-username"),
        ("password", "your-password"),
        ("query", "https://www.space-track.org/basicspacedata/query/class/gp/orderby/LAUNCH_DATE desc/limit/20/emptyresult/show"),
    ]);
    if response.error() {
        Err(sgp4::Error::new(format!(
            "network error {}: {}",
            response.status(),
            response.into_string()?
        )))
    } else {
        let elements_group: Vec<sgp4::Elements> = response.into_json_deserialize()?;
        for elements in &elements_group {
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
}
