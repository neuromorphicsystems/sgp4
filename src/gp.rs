use chrono::{Datelike, Timelike};

#[derive(Debug, Clone)]
pub struct Error {
    message: String,
}

impl Error {
    pub fn new(message: String) -> Error {
        Error { message: message }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "{}", self.message)
    }
}

impl From<std::str::Utf8Error> for Error {
    fn from(error: std::str::Utf8Error) -> Self {
        Error::new(error.to_string())
    }
}

impl From<std::num::ParseIntError> for Error {
    fn from(error: std::num::ParseIntError) -> Self {
        Error::new(error.to_string())
    }
}

impl From<std::num::ParseFloatError> for Error {
    fn from(error: std::num::ParseFloatError) -> Self {
        Error::new(error.to_string())
    }
}

impl From<serde_json::Error> for Error {
    fn from(error: serde_json::Error) -> Self {
        Error::new(error.to_string())
    }
}

trait DecimalPointAssumedRepresentation {
    fn parse_decimal_point_assumed(&self) -> Result<f64>;
}

impl DecimalPointAssumedRepresentation for [u8] {
    fn parse_decimal_point_assumed(&self) -> Result<f64> {
        let trimmed = std::str::from_utf8(self)?.trim_start();
        if trimmed.starts_with("-") {
            Ok(format!("-.{}", &trimmed[1..]).parse::<f64>()?)
        } else {
            Ok(format!(".{}", trimmed).parse::<f64>()?)
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

#[derive(serde::Serialize, serde::Deserialize)]
pub enum Classification {
    #[serde(rename = "U")]
    Unclassified,

    #[serde(rename = "C")]
    Classified,

    #[serde(rename = "S")]
    Secret,
}

#[derive(serde::Serialize, serde::Deserialize)]
pub struct Elements {
    #[serde(rename = "OBJECT_NAME")]
    pub object_name: Option<String>,

    #[serde(rename = "OBJECT_ID")]
    pub international_designator: Option<String>,

    #[serde(rename = "NORAD_CAT_ID")]
    pub norad_id: u64,

    #[serde(rename = "CLASSIFICATION_TYPE")]
    pub classification: Classification,

    #[serde(rename = "EPOCH")]
    pub datetime: chrono::naive::NaiveDateTime,

    #[serde(rename = "MEAN_MOTION_DOT")]
    pub mean_motion_dot: f64,

    #[serde(rename = "MEAN_MOTION_DDOT")]
    pub mean_motion_ddot: f64,

    #[serde(rename = "BSTAR")]
    pub drag_term: f64,

    #[serde(rename = "ELEMENT_SET_NO")]
    pub element_set_number: u64,

    #[serde(rename = "INCLINATION")]
    pub inclination: f64,

    #[serde(rename = "RA_OF_ASC_NODE")]
    pub right_ascension: f64,

    #[serde(rename = "ECCENTRICITY")]
    pub eccentricity: f64,

    #[serde(rename = "ARG_OF_PERICENTER")]
    pub argument_of_perigee: f64,

    #[serde(rename = "MEAN_ANOMALY")]
    pub mean_anomaly: f64,

    #[serde(rename = "MEAN_MOTION")]
    pub mean_motion: f64,

    #[serde(rename = "REV_AT_EPOCH")]
    pub revolution_number: u64,

    #[serde(rename = "EPHEMERIS_TYPE")]
    pub ephemeris_type: u8,
}

impl Elements {
    pub fn from_tle(object_name: Option<String>, line1: &[u8], line2: &[u8]) -> Result<Elements> {
        if line1.len() != 69 {
            return Err(Error::new("line 1 must have 69 characters".to_owned()));
        }
        if line2.len() != 69 {
            return Err(Error::new("line 2 must have 69 characters".to_owned()));
        }
        if line1[0] != b'1' {
            return Err(Error::new(
                "line 1 must start with the character '1'".to_owned(),
            ));
        }
        if line2[0] != b'2' {
            return Err(Error::new(
                "line 2 must start with the character '2'".to_owned(),
            ));
        }
        for index in [1, 8, 17, 32, 43, 52, 61, 63].iter() {
            if line1[*index] != b' ' {
                return Err(Error::new(format!(
                    "line 1:{} must be a space character",
                    index + 1
                )));
            }
        }
        for index in [1, 7, 16, 25, 33, 42, 51].iter() {
            if line2[*index] != b' ' {
                return Err(Error::new(format!(
                    "line 2:{} must be a space character",
                    index + 1
                )));
            }
        }
        let norad_id = std::str::from_utf8(&line1[2..7])?.parse::<u64>()?;
        if norad_id != std::str::from_utf8(&line2[2..7])?.parse::<u64>()? {
            return Err(Error::new(
                "line 1 and 2 have different satellite numbers".to_owned(),
            ));
        }
        for line in &[line1, line2] {
            if (line[..68]
                .iter()
                .fold(0, |accumulator, character| match character {
                    b'-' => accumulator + 1,
                    character if character >= &b'0' && character <= &b'9' => {
                        accumulator + (character - b'0') as u16
                    }
                    _ => accumulator,
                })
                % 10) as u8
                != line[68] - b'0'
            {
                return Err(Error::new("bad checksum".to_owned()));
            }
        }
        Ok(Elements {
            object_name: object_name,
            norad_id: norad_id,
            classification: match line1[7] {
                b'U' => Classification::Unclassified,
                b'C' => Classification::Classified,
                b'S' => Classification::Secret,
                _ => return Err(Error::new("unknown classification".to_owned())),
            },
            international_designator: if line1[9..17]
                .iter()
                .all(|character| *character == ' ' as u8)
            {
                None
            } else {
                Some(format!(
                    "{}-{}",
                    match std::str::from_utf8(&line1[9..11])?.parse::<u8>()? {
                        launch_year if launch_year < 57 => 2000 + launch_year as u16,
                        launch_year => 1900 + launch_year as u16,
                    },
                    std::str::from_utf8(&line1[11..17])?.trim()
                ))
            },
            datetime: {
                let day = std::str::from_utf8(&line1[20..32])?
                    .trim_start()
                    .parse::<f64>()?;
                let seconds = day.fract() * (24.0 * 60.0 * 60.0);
                chrono::NaiveDate::from_yo(
                    match std::str::from_utf8(&line1[18..20])?.parse::<u8>()? {
                        year if year < 57 => year as i32 + 2000,
                        year => year as i32 + 1900,
                    },
                    day as u32,
                )
                .and_time(chrono::NaiveTime::from_num_seconds_from_midnight(
                    seconds as u32,
                    (seconds.fract() * 1e9).round() as u32,
                ))
            },
            mean_motion_dot: std::str::from_utf8(&line1[33..43])?.trim_start().parse()?,
            mean_motion_ddot: line1[44..50].parse_decimal_point_assumed()?
                * 10.0_f64.powi(std::str::from_utf8(&line1[50..52])?.parse::<i8>()? as i32),
            drag_term: line1[53..59].parse_decimal_point_assumed()?
                * 10.0_f64.powi(std::str::from_utf8(&line1[59..61])?.parse::<i8>()? as i32),
            ephemeris_type: std::str::from_utf8(&line1[62..63])?.trim_start().parse()?,
            element_set_number: std::str::from_utf8(&line1[64..68])?.trim_start().parse()?,
            inclination: std::str::from_utf8(&line2[8..16])?.trim_start().parse()?,
            right_ascension: std::str::from_utf8(&line2[17..25])?.trim_start().parse()?,
            eccentricity: line2[26..33].parse_decimal_point_assumed()?,
            argument_of_perigee: std::str::from_utf8(&line2[34..42])?.trim_start().parse()?,
            mean_anomaly: std::str::from_utf8(&line2[43..51])?.trim_start().parse()?,
            mean_motion: std::str::from_utf8(&line2[52..63])?.trim_start().parse()?,
            revolution_number: std::str::from_utf8(&line2[63..68])?.trim_start().parse()?,
        })
    }

    pub fn epoch(&self) -> f64 {
        // y₂₀₀₀ = (367 yᵤ - ⌊7 (yᵤ + ⌊(mᵤ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤ / 9⌋ + dᵤ - 730531) / 365.25
        //         + (3600 hᵤ + 60 minᵤ + sᵤ - 43200) / (24 × 60 × 60 × 365.25)
        //         + nsᵤ / (24 × 60 × 60 × 365.25 × 10⁹)
        (367 * self.datetime.year() as i32
            - (7 * (self.datetime.year() as i32 + (self.datetime.month() as i32 + 9) / 12)) / 4
            + 275 * self.datetime.month() as i32 / 9
            + self.datetime.day() as i32 - 730531) as f64 / 365.25
            + (self.datetime.num_seconds_from_midnight() as i32 - 43200) as f64 / (24.0 * 60.0 * 60.0 * 365.25)
            + (self.datetime.nanosecond() as f64) / (24.0 * 60.0 * 60.0 * 1e9 * 365.25)
    }

    pub fn epoch_afspc_compatibility_mode(&self) -> f64 {
        // y₂₀₀₀ = (367 yᵤ - ⌊7 (yᵤ + ⌊(mᵤ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤ / 9⌋ + dᵤ
        //         + 1721013.5
        //         + hᵤ / 24
        //         + minᵤ / (24 × 60)
        //         + sᵤ / (24 × 60 × 60)
        //         + nsᵤ / (24 × 60 × 60 × 10⁹)
        //         - 2451545)
        //         / 365.25
        (((367 * self.datetime.year() as u32
            - (7 * (self.datetime.year() as u32 + (self.datetime.month() + 9) / 12)) / 4
            + 275 * self.datetime.month() / 9
            + self.datetime.day()) as f64
            + 1721013.5
            + self.datetime.hour() as f64 / 24.0
            + self.datetime.minute() as f64 / (24.0 * 60.0)
            + self.datetime.second() as f64 / (24.0 * 60.0 * 60.0)
            + self.datetime.nanosecond() as f64 / (24.0 * 60.0 * 60.0 * 1e9))
            - 2451545.0)
            / 365.25
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_eq_f64(first: f64, second: f64) {
        if second == 0.0 {
            assert_eq!(first, 0.0);
        } else {
            assert!((first - second).abs() / second < f64::EPSILON);
        }
    }

    #[test]
    fn test_from_tle() -> Result<()> {
        let elements = Elements::from_tle(
            Some("ISS (ZARYA)".to_owned()),
            "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927".as_bytes(),
            "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537".as_bytes(),
        )?;
        match elements.object_name.as_ref() {
            Some(object_name) => assert_eq!(object_name, "ISS (ZARYA)"),
            None => panic!(),
        }
        assert_eq!(elements.norad_id, 25544);
        assert!(matches!(
            elements.classification,
            Classification::Unclassified
        ));
        assert_eq!(
            elements.international_designator.as_ref().unwrap(),
            "1998-067A"
        );
        assert_eq!(
            elements.datetime,
            chrono::NaiveDate::from_yo(2008, 264).and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight(44740, 104192001)
            )
        );
        assert_eq_f64(elements.epoch(), 8.720103559972621);
        assert_eq_f64(
            elements.epoch_afspc_compatibility_mode(),
            8.7201035599722125,
        );
        assert_eq_f64(elements.mean_motion_dot, -0.00002182);
        assert_eq_f64(elements.mean_motion_ddot, 0.0);
        assert_eq_f64(elements.drag_term, -0.11606e-4);
        assert_eq!(elements.ephemeris_type, 0);
        assert_eq!(elements.element_set_number, 292);
        assert_eq_f64(elements.inclination, 51.6416);
        assert_eq_f64(elements.right_ascension, 247.4627);
        assert_eq_f64(elements.eccentricity, 0.0006703);
        assert_eq_f64(elements.argument_of_perigee, 130.5360);
        assert_eq_f64(elements.mean_anomaly, 325.0288);
        assert_eq_f64(elements.mean_motion, 15.72125391);
        assert_eq!(elements.revolution_number, 56353);
        let elements = Elements::from_tle(
            None,
            "1 11801U          80230.29629788  .01431103  00000-0  14311-1 0    13".as_bytes(),
            "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848    13".as_bytes(),
        )?;
        assert!(elements.object_name.is_none());
        assert_eq!(elements.norad_id, 11801);
        assert!(matches!(
            elements.classification,
            Classification::Unclassified
        ));
        assert!(elements.international_designator.is_none());
        assert_eq!(
            elements.datetime,
            chrono::NaiveDate::from_yo(1980, 230).and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight(25600, 136832000)
            )
        );
        assert_eq_f64(elements.epoch(), -19.373589875756331);
        assert_eq_f64(
            elements.epoch_afspc_compatibility_mode(),
            -19.373589875756632,
        );
        assert_eq_f64(elements.mean_motion_dot, 0.01431103);
        assert_eq_f64(elements.mean_motion_ddot, 0.0);
        assert_eq_f64(elements.drag_term, 0.014311);
        assert_eq!(elements.ephemeris_type, 0);
        assert_eq!(elements.element_set_number, 1);
        assert_eq_f64(elements.inclination, 46.7916);
        assert_eq_f64(elements.right_ascension, 230.4354);
        assert_eq_f64(elements.eccentricity, 0.7318036);
        assert_eq_f64(elements.argument_of_perigee, 47.4722);
        assert_eq_f64(elements.mean_anomaly, 10.4117);
        assert_eq_f64(elements.mean_motion, 2.28537848);
        assert_eq!(elements.revolution_number, 1);
        Ok(())
    }

    #[test]
    fn test_from_omm() -> Result<()> {
        let elements: Elements = serde_json::from_str(
            r#"{
            "OBJECT_NAME": "ISS (ZARYA)",
            "OBJECT_ID": "1998-067A",
            "EPOCH": "2020-07-12T01:19:07.402656",
            "MEAN_MOTION": 15.49560532,
            "ECCENTRICITY": 0.0001771,
            "INCLINATION": 51.6435,
            "RA_OF_ASC_NODE": 225.4004,
            "ARG_OF_PERICENTER": 44.9625,
            "MEAN_ANOMALY": 5.1087,
            "EPHEMERIS_TYPE": 0,
            "CLASSIFICATION_TYPE": "U",
            "NORAD_CAT_ID": 25544,
            "ELEMENT_SET_NO": 999,
            "REV_AT_EPOCH": 23587,
            "BSTAR": 0.0049645,
            "MEAN_MOTION_DOT": 0.00289036,
            "MEAN_MOTION_DDOT": 0
        }"#,
        )?;
        match elements.object_name.as_ref() {
            Some(object_name) => assert_eq!(object_name, "ISS (ZARYA)"),
            None => panic!(),
        }
        assert_eq!(elements.norad_id, 25544);
        assert!(matches!(
            elements.classification,
            Classification::Unclassified
        ));
        assert_eq!(
            elements.international_designator.as_ref().unwrap(),
            "1998-067A"
        );
        assert_eq!(
            elements.datetime,
            chrono::NaiveDate::from_yo(2020, 194).and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight(4747, 402656000)
            )
        );
        assert_eq_f64(elements.epoch(), 20.527186712635181);
        assert_eq_f64(
            elements.epoch_afspc_compatibility_mode(),
            20.527186712635135,
        );
        assert_eq_f64(elements.mean_motion_dot, 0.00289036);
        assert_eq_f64(elements.mean_motion_ddot, 0.0);
        assert_eq_f64(elements.drag_term, 0.0049645);
        assert_eq!(elements.ephemeris_type, 0);
        assert_eq!(elements.element_set_number, 999);
        assert_eq_f64(elements.inclination, 51.6435);
        assert_eq_f64(elements.right_ascension, 225.4004);
        assert_eq_f64(elements.eccentricity, 0.0001771);
        assert_eq_f64(elements.argument_of_perigee, 44.9625);
        assert_eq_f64(elements.mean_anomaly, 5.1087);
        assert_eq_f64(elements.mean_motion, 15.49560532);
        assert_eq!(elements.revolution_number, 23587);
        Ok(())
    }
}
