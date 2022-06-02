fn main() -> sgp4::Result<()> {
    match ureq::post("https://www.space-track.org/ajaxauth/login").send_form(&[
        ("identity", "your-username"),
        ("password", "your-password"),
        ("query", "https://www.space-track.org/basicspacedata/query/class/gp/orderby/LAUNCH_DATE desc/limit/20/emptyresult/show"),
    ]) {
        Ok(response) => {
            let elements_group: Vec<sgp4::Elements> = response.into_json()?;
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
        },
        Err(ureq::Error::Status(code, response)) => Err(sgp4::Error::new(format!(
            "network error {}: {}",
            code,
            response.into_string()?
        ))),
        Err(error) => Err(sgp4::Error::new(error.to_string()))
    }
}
