fn main() -> anyhow::Result<()> {
    let response = ureq::post("https://www.space-track.org/ajaxauth/login").send_form(&[
        ("identity", "your-username"),
        ("password", "your-password"),
        ("query", "https://www.space-track.org/basicspacedata/query/class/gp/decay_date/null-val/epoch/<now-30/orderby/norad_cat_id/limit/20/format/json"),
    ])?;
    let elements_vec: Vec<sgp4::Elements> = response.into_json()?;
    for elements in &elements_vec {
        println!("{}", elements.object_name.as_ref().unwrap());
        let constants = sgp4::Constants::from_elements(elements)?;
        for hours in &[12, 24] {
            println!("    t = {} min", hours * 60);
            let prediction = constants.propagate(sgp4::MinutesSinceEpoch((hours * 60) as f64))?;
            println!("        r = {:?} km", prediction.position);
            println!("        ṙ = {:?} km.s⁻¹", prediction.velocity);
        }
    }
    Ok(())
}
