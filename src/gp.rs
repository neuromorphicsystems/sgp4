use chrono::{Datelike, Timelike};
use serde::de::Deserialize;

/// Represents an SGP4 error
///
/// Errors can result from corrupted TLEs or OMMs, or if one of the orbital elements diverges during propagation.
#[derive(Debug, Clone)]
pub struct Error {
    message: String,
}

impl Error {
    /// Creates a new error from a string
    ///
    /// # Arguments
    ///
    /// * `message` - A description of the error
    ///
    /// # Examples
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// #     if false {
    /// Err(sgp4::Error::new("a useful message".to_owned()))
    /// #     } else {
    /// #         Ok::<(), sgp4::Error>(())
    /// #     }
    /// # }
    /// ```
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// #     if false {
    /// Err(sgp4::Error::new(format!("error code {}", 3)))
    /// #     } else {
    /// #         Ok::<(), sgp4::Error>(())
    /// #     }
    /// # }
    /// ```
    pub fn new(message: String) -> Error {
        Error { message }
    }
}

impl std::error::Error for Error {}

impl std::fmt::Display for Error {
    fn fmt(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(formatter, "{}", self.message)
    }
}

impl From<std::io::Error> for Error {
    fn from(error: std::io::Error) -> Self {
        Error::new(error.to_string())
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
        if let Some(stripped) = trimmed.strip_prefix('-') {
            Ok(format!("-.{}", stripped).parse::<f64>()?)
        } else if let Some(stripped) = trimmed.strip_prefix('+') {
            Ok(format!(".{}", stripped).parse::<f64>()?)
        } else {
            Ok(format!(".{}", trimmed).parse::<f64>()?)
        }
    }
}

/// The result type returned by SGP4 functions
pub type Result<T> = std::result::Result<T, Error>;

/// A satellite's elements classification
#[derive(serde::Serialize, serde::Deserialize)]
pub enum Classification {
    /// Declassfied objects or objects without a classification
    #[serde(rename = "U")]
    Unclassified,

    /// Would cause "serious damage" to national security if it were publicly available
    #[serde(rename = "C")]
    Classified,

    /// Would cause "damage" or be prejudicial to national security if publicly available
    #[serde(rename = "S")]
    Secret,
}

/// General perturbations orbital data parsed from a TLE or OMM
///
/// Elements can be retrieved from either a Two-Line Element Set (TLE) or an Orbit Mean-Elements Message (OMM).
/// See [https://celestrak.com/NORAD/documentation/gp-data-formats.php](https://celestrak.com/NORAD/documentation/gp-data-formats.php)
/// for more information on the difference between the two formats.
///
/// The fields' documentation is adapted from [https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/SSOP_Help/tle_def.html](https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/SSOP_Help/tle_def.html).
///
/// See [sgp4::Elements::from_tle](struct.Elements.html#method.from_tle) to parse a TLE.
///
/// `serde_json` can be used to parse a JSON OMM object (into a `sgp4::Elements`)
///  or a JSON list of OMM objects (into a `Vec<sgp4::Elements>`).
///
/// # Example
/// ```
/// # fn main() -> sgp4::Result<()> {
/// let elements: sgp4::Elements = serde_json::from_str(
///     r#"{
///         "OBJECT_NAME": "ISS (ZARYA)",
///         "OBJECT_ID": "1998-067A",
///         "EPOCH": "2020-07-12T01:19:07.402656",
///         "MEAN_MOTION": 15.49560532,
///         "ECCENTRICITY": 0.0001771,
///         "INCLINATION": 51.6435,
///         "RA_OF_ASC_NODE": 225.4004,
///         "ARG_OF_PERICENTER": 44.9625,
///         "MEAN_ANOMALY": 5.1087,
///         "EPHEMERIS_TYPE": 0,
///         "CLASSIFICATION_TYPE": "U",
///         "NORAD_CAT_ID": 25544,
///         "ELEMENT_SET_NO": 999,
///         "REV_AT_EPOCH": 23587,
///         "BSTAR": 0.0049645,
///         "MEAN_MOTION_DOT": 0.00289036,
///         "MEAN_MOTION_DDOT": 0
///     }"#,
/// )?;
/// #     Ok(())
/// # }
#[derive(serde::Serialize, serde::Deserialize)]
pub struct Elements {
    /// The name associated with the satellite
    #[serde(rename = "OBJECT_NAME")]
    pub object_name: Option<String>,

    /// The satellite's international designator
    ///
    /// It consists of the launch year, the launch number of that year and
    /// a letter code representing the sequential identifier of a piece in a launch.
    #[serde(rename = "OBJECT_ID")]
    pub international_designator: Option<String>,

    /// The catalog number USSPACECOM has designated for this object
    #[serde(rename = "NORAD_CAT_ID", deserialize_with = "u64_or_string")]
    pub norad_id: u64,

    /// The elements' classification
    #[serde(rename = "CLASSIFICATION_TYPE")]
    pub classification: Classification,

    /// The UTC timestamp of the elements
    #[serde(rename = "EPOCH")]
    pub datetime: chrono::naive::NaiveDateTime,

    /// Time derivative of the mean motion
    #[serde(rename = "MEAN_MOTION_DOT", deserialize_with = "f64_or_string")]
    pub mean_motion_dot: f64,

    /// Second time derivative of the mean motion
    #[serde(rename = "MEAN_MOTION_DDOT", deserialize_with = "f64_or_string")]
    pub mean_motion_ddot: f64,

    /// Radiation pressure coefficient in earth radii⁻¹
    #[serde(rename = "BSTAR", deserialize_with = "f64_or_string")]
    pub drag_term: f64,

    /// A running count of all 2 line element sets generated by USSPACECOM for this object
    #[serde(rename = "ELEMENT_SET_NO", deserialize_with = "u64_or_string")]
    pub element_set_number: u64,

    /// Angle between the equator and the orbit plane in deg
    #[serde(rename = "INCLINATION", deserialize_with = "f64_or_string")]
    pub inclination: f64,

    /// Angle between vernal equinox and the point where the orbit crosses the equatorial plane in deg
    #[serde(rename = "RA_OF_ASC_NODE", deserialize_with = "f64_or_string")]
    pub right_ascension: f64,

    /// The shape of the orbit
    #[serde(rename = "ECCENTRICITY", deserialize_with = "f64_or_string")]
    pub eccentricity: f64,

    /// Angle between the ascending node and the orbit's point of closest approach to the earth in deg
    #[serde(rename = "ARG_OF_PERICENTER", deserialize_with = "f64_or_string")]
    pub argument_of_perigee: f64,

    /// Angle of the satellite location measured from perigee in deg
    #[serde(rename = "MEAN_ANOMALY", deserialize_with = "f64_or_string")]
    pub mean_anomaly: f64,

    /// Mean number of orbits per day in day⁻¹ (Kozai convention)
    #[serde(rename = "MEAN_MOTION", deserialize_with = "f64_or_string")]
    pub mean_motion: f64,

    /// The orbit number at epoch
    #[serde(rename = "REV_AT_EPOCH", deserialize_with = "u64_or_string")]
    pub revolution_number: u64,

    /// NORAD internal use, always 0 in distributed data
    #[serde(rename = "EPHEMERIS_TYPE", deserialize_with = "u8_or_string")]
    pub ephemeris_type: u8,
}

fn u64_or_string<'de, D>(deserializer: D) -> std::result::Result<u64, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    match serde_json::value::Value::deserialize(deserializer)? {
        serde_json::value::Value::Number(number) => number
            .as_u64()
            .ok_or_else(|| serde::de::Error::custom("parsing the number as u64 failed")),
        serde_json::value::Value::String(string) => {
            string.parse().map_err(serde::de::Error::custom)
        }
        _ => Err(serde::de::Error::custom("expected a u64 or string")),
    }
}

fn u8_or_string<'de, D>(deserializer: D) -> std::result::Result<u8, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    match serde_json::value::Value::deserialize(deserializer)? {
        serde_json::value::Value::Number(number) => match number.as_u64() {
            Some(value) => Ok(value as u8),
            None => Err(serde::de::Error::custom("parsing the number as u64 failed")),
        },
        serde_json::value::Value::String(string) => {
            string.parse().map_err(serde::de::Error::custom)
        }
        _ => Err(serde::de::Error::custom("expected a u64 or string")),
    }
}

fn f64_or_string<'de, D>(deserializer: D) -> std::result::Result<f64, D::Error>
where
    D: serde::de::Deserializer<'de>,
{
    match serde_json::value::Value::deserialize(deserializer)? {
        serde_json::value::Value::Number(number) => number
            .as_f64()
            .ok_or_else(|| serde::de::Error::custom("parsing the number as f64 failed")),
        serde_json::value::Value::String(string) => {
            string.parse().map_err(serde::de::Error::custom)
        }
        _ => Err(serde::de::Error::custom("expected a f64 or string")),
    }
}

impl Elements {
    /// Parses a Two-Line Element Set (TLE) with an optionnal title
    ///
    /// # Arguments
    ///
    /// * `object_name` - The name of the satellite, usually given by a third line placed before the TLE
    /// * `line1` - The first line of the TLE composed of 69 ASCII characters
    /// * `line2` - The second line of the TLE composed of 69 ASCII characters
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() -> sgp4::Result<()> {
    /// let elements = sgp4::Elements::from_tle(
    ///     Some("ISS (ZARYA)".to_owned()),
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927".as_bytes(),
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537".as_bytes(),
    /// )?;
    /// #     Ok(())
    /// # }
    /// ```
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
        let norad_id = std::str::from_utf8(&line1[2..7])?
            .trim_start()
            .parse::<u64>()?;
        if norad_id
            != std::str::from_utf8(&line2[2..7])?
                .trim_start()
                .parse::<u64>()?
        {
            return Err(Error::new(
                "line 1 and 2 have different satellite numbers".to_owned(),
            ));
        }
        for line in &[line1, line2] {
            if (line[..68]
                .iter()
                .fold(0, |accumulator, character| match character {
                    b'-' => accumulator + 1,
                    character if (&b'0'..=&b'9').contains(&character) => {
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
            object_name,
            norad_id,
            classification: match line1[7] {
                b'U' => Classification::Unclassified,
                b'C' => Classification::Classified,
                b'S' => Classification::Secret,
                _ => return Err(Error::new("unknown classification".to_owned())),
            },
            international_designator: if line1[9..17].iter().all(|character| *character == b' ') {
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

    /// Returns the number of years since UTC 1 January 2000 12h00 (J2000)
    ///
    /// This is the recommended method to calculate the epoch
    pub fn epoch(&self) -> f64 {
        // y₂₀₀₀ = (367 yᵤ - ⌊7 (yᵤ + ⌊(mᵤ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤ / 9⌋ + dᵤ - 730531) / 365.25
        //         + (3600 hᵤ + 60 minᵤ + sᵤ - 43200) / (24 × 60 × 60 × 365.25)
        //         + nsᵤ / (24 × 60 × 60 × 365.25 × 10⁹)
        (367 * self.datetime.year() as i32
            - (7 * (self.datetime.year() as i32 + (self.datetime.month() as i32 + 9) / 12)) / 4
            + 275 * self.datetime.month() as i32 / 9
            + self.datetime.day() as i32
            - 730531) as f64
            / 365.25
            + (self.datetime.num_seconds_from_midnight() as i32 - 43200) as f64
                / (24.0 * 60.0 * 60.0 * 365.25)
            + (self.datetime.nanosecond() as f64) / (24.0 * 60.0 * 60.0 * 1e9 * 365.25)
    }

    /// Returns the number of years since UTC 1 January 2000 12h00 (J2000) using the AFSPC expression
    ///
    /// This function should be used if compatibility with the AFSPC implementation is needed
    pub fn epoch_afspc_compatibility_mode(&self) -> f64 {
        // y₂₀₀₀ = (367 yᵤ - ⌊7 (yᵤ + ⌊(mᵤ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤ / 9⌋ + dᵤ
        //         + 1721013.5
        //         + (((nsᵤ / 10⁹ + sᵤ) / 60 + minᵤ) / 60 + hᵤ) / 24
        //         - 2451545)
        //         / 365.25
        ((367 * self.datetime.year() as u32
            - (7 * (self.datetime.year() as u32 + (self.datetime.month() + 9) / 12)) / 4
            + 275 * self.datetime.month() / 9
            + self.datetime.day()) as f64
            + 1721013.5
            + (((self.datetime.nanosecond() as f64 / 1e9 + self.datetime.second() as f64) / 60.0
                + self.datetime.minute() as f64)
                / 60.0
                + self.datetime.hour() as f64)
                / 24.0
            - 2451545.0)
            / 365.25
    }
}

/// Parses a multi-line TL/2LE string into a list of `Elements`
///
/// Each pair of lines must represent a TLE, for example as in
/// [https://celestrak.com/NORAD/elements/gp.php?GROUP=stations&FORMAT=2le](https://celestrak.com/NORAD/elements/gp.php?GROUP=stations&FORMAT=2le).
///
/// # Arguments
///
/// * `tles` - A string containing multiple lines
pub fn parse_2les(tles: &str) -> Result<Vec<Elements>> {
    let mut line_buffer = "";
    let mut first = true;
    let mut elements_group = Vec::new();
    for line in tles.lines() {
        if first {
            line_buffer = line;
        } else {
            elements_group.push(Elements::from_tle(
                None,
                line_buffer.as_bytes(),
                line.as_bytes(),
            )?);
        }
        first = !first;
    }
    Ok(elements_group)
}

/// Parses a multi-line TL/3LE string into a list of `Elements`
///
/// Each triplet of lines must represent a TLE with an object name, for example as in
/// [https://celestrak.com/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle](https://celestrak.com/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle).
///
/// # Arguments
///
/// * `tles` - A string containing multiple lines
pub fn parse_3les(tles: &str) -> Result<Vec<Elements>> {
    let mut lines_buffer = ["", ""];
    let mut index = 0;
    let mut elements_group = Vec::new();
    for line in tles.lines() {
        match index {
            0 | 1 => {
                lines_buffer[index] = line;
                index += 1;
            }
            _ => {
                elements_group.push(Elements::from_tle(
                    Some(lines_buffer[0].to_owned()),
                    lines_buffer[1].as_bytes(),
                    line.as_bytes(),
                )?);
                index = 0;
            }
        }
    }
    Ok(elements_group)
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
    fn test_from_celestrak_omm() -> Result<()> {
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

    #[test]
    fn test_from_space_track_omm() -> Result<()> {
        let elements: Elements = serde_json::from_str(
            r#"{"CCSDS_OMM_VERS":"2.0",
                "COMMENT":"GENERATED VIA SPACE-TRACK.ORG API",
                "CREATION_DATE":"2020-12-13T17:26:09",
                "ORIGINATOR":"18 SPCS",
                "OBJECT_NAME":"ISS (ZARYA)",
                "OBJECT_ID":"1998-067A",
                "CENTER_NAME":"EARTH",
                "REF_FRAME":"TEME",
                "TIME_SYSTEM":"UTC",
                "MEAN_ELEMENT_THEORY":"SGP4",
                "EPOCH":"2020-12-13T16:36:04.502592",
                "MEAN_MOTION":"15.49181153",
                "ECCENTRICITY":"0.00017790",
                "INCLINATION":"51.6444",
                "RA_OF_ASC_NODE":"180.2777",
                "ARG_OF_PERICENTER":"128.5985",
                "MEAN_ANOMALY":"350.1361",
                "EPHEMERIS_TYPE":"0",
                "CLASSIFICATION_TYPE":"U",
                "NORAD_CAT_ID":"25544",
                "ELEMENT_SET_NO":"999",
                "REV_AT_EPOCH":"25984",
                "BSTAR":"0.00002412400000",
                "MEAN_MOTION_DOT":"0.00000888",
                "MEAN_MOTION_DDOT":"0.0000000000000",
                "SEMIMAJOR_AXIS":"6797.257",
                "PERIOD":"92.952",
                "APOAPSIS":"420.331",
                "PERIAPSIS":"417.912",
                "OBJECT_TYPE":"PAYLOAD",
                "RCS_SIZE":"LARGE",
                "COUNTRY_CODE":"ISS",
                "LAUNCH_DATE":"1998-11-20",
                "SITE":"TTMTR",
                "DECAY_DATE":null,
                "FILE":"2902442",
                "GP_ID":"167697146",
                "TLE_LINE0":"0 ISS (ZARYA)",
                "TLE_LINE1":"1 25544U 98067A   20348.69171878  .00000888  00000-0  24124-4 0  9995",
                "TLE_LINE2":"2 25544  51.6444 180.2777 0001779 128.5985 350.1361 15.49181153259845"
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
            chrono::NaiveDate::from_yo(2020, 348).and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight(59764, 502592000)
            )
        );
        assert_eq_f64(elements.epoch(), 20.95055912054757);
        assert_eq_f64(elements.epoch_afspc_compatibility_mode(), 20.95055912054749);
        assert_eq_f64(elements.mean_motion_dot, 0.00000888);
        assert_eq_f64(elements.mean_motion_ddot, 0.0);
        assert_eq_f64(elements.drag_term, 0.000024124);
        assert_eq!(elements.ephemeris_type, 0);
        assert_eq!(elements.element_set_number, 999);
        assert_eq_f64(elements.inclination, 51.6444);
        assert_eq_f64(elements.right_ascension, 180.2777);
        assert_eq_f64(elements.eccentricity, 0.0001779);
        assert_eq_f64(elements.argument_of_perigee, 128.5985);
        assert_eq_f64(elements.mean_anomaly, 350.1361);
        assert_eq_f64(elements.mean_motion, 15.49181153);
        assert_eq!(elements.revolution_number, 25984);
        Ok(())
    }

    #[test]
    fn test_from_celestrak_omms() -> Result<()> {
        let elements_group: Vec<Elements> = serde_json::from_str(
            r#"[{
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
            },{
                "OBJECT_NAME": "KESTREL EYE IIM (KE2M)",
                "OBJECT_ID": "1998-067NE",
                "EPOCH": "2020-07-12T01:38:52.903968",
                "MEAN_MOTION": 15.70564504,
                "ECCENTRICITY": 0.0002758,
                "INCLINATION": 51.6338,
                "RA_OF_ASC_NODE": 155.6245,
                "ARG_OF_PERICENTER": 166.8841,
                "MEAN_ANOMALY": 193.2228,
                "EPHEMERIS_TYPE": 0,
                "CLASSIFICATION_TYPE": "U",
                "NORAD_CAT_ID": 42982,
                "ELEMENT_SET_NO": 999,
                "REV_AT_EPOCH": 15494,
                "BSTAR": 7.2204e-5,
                "MEAN_MOTION_DOT": 8.489e-5,
                "MEAN_MOTION_DDOT": 0
            }]"#,
        )?;
        assert_eq!(elements_group.len(), 2);
        Ok(())
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
    fn test_parse_2les() -> Result<()> {
        let elements_group = parse_2les(
            "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992\n\
             2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008\n\
             1 42982U 98067NE  20194.06866787  .00008489  00000-0  72204-4 0  9997\n\
             2 42982  51.6338 155.6245 0002758 166.8841 193.2228 15.70564504154944\n",
        )?;
        assert_eq!(elements_group.len(), 2);
        Ok(())
    }

    #[test]
    fn test_parse_3les() -> Result<()> {
        let elements_group = parse_3les(
            "ISS (ZARYA)\n\
             1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992\n\
             2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008\n\
             KESTREL EYE IIM (KE2M)\n\
             1 42982U 98067NE  20194.06866787  .00008489  00000-0  72204-4 0  9997\n\
             2 42982  51.6338 155.6245 0002758 166.8841 193.2228 15.70564504154944\n",
        )?;
        assert_eq!(elements_group.len(), 2);
        Ok(())
    }
}
