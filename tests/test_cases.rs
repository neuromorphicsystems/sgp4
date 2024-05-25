#[derive(serde::Serialize, serde::Deserialize)]
#[serde(untagged)]
pub enum State {
    Ok {
        time: f64,
        position: [f64; 3],
        velocity: [f64; 3],
        date: Option<toml::value::Datetime>,
    },
    Err {
        time: f64,
        error: String,
    },
}

#[derive(serde::Serialize, serde::Deserialize)]
pub struct TestCase {
    pub line1: String,
    pub line2: String,
    pub states: Vec<State>,
}

#[derive(serde::Serialize, serde::Deserialize)]
pub struct TestCases {
    pub list: Vec<TestCase>,
}
