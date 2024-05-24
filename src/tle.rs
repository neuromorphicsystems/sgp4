use chrono::{Datelike, Timelike};

#[cfg(feature = "alloc")]
use alloc::format;

#[cfg(feature = "alloc")]
use alloc::borrow::ToOwned;

#[cfg(not(feature = "std"))]
use num_traits::Float;

#[cfg(feature = "serde")]
use serde::de::Deserialize;

/// TLE error type
#[derive(Debug, Clone)]
pub enum ErrorWhat {
    BadChecksum,
    BadLength,
    BadFirstCharacter,
    ExpectedFloat,
    ExpectedFloatWithAssumedDecimalPoint,
    ExpectedInteger,
    ExpectedSpace,
    ExpectedString,
    FloatWithAssumedDecimalPointTooLong,
    NoradIdMismatch,
    UnknownClassification,
    FromYoOptFailed,
    FromNumSecondsFromMidnightFailed,
}

/// Input line where a parse error was found
#[derive(Debug, Clone)]
pub enum ErrorLine {
    Line1,
    Line2,
    Both,
}

/// Represents a TLE parse error
#[derive(Debug, Clone)]
pub struct Error {
    /// TLE error type
    what: ErrorWhat,

    /// TLE error line
    line: ErrorLine,

    /// Start character position of the line slice that caused the error
    start: usize,

    /// End character position of the line slice that caused the error
    end: usize,
}

impl core::fmt::Display for Error {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_fmt(format_args!("TLE parse error: {} {} between characters {} and {}",
            match self.what {
                ErrorWhat::BadChecksum => "Bad line checksum",
                ErrorWhat::BadLength => "Bad line length",
                ErrorWhat::BadFirstCharacter => "Bad first character",
                ErrorWhat::ExpectedFloat => "Parsing a float field failed",
                ErrorWhat::ExpectedFloatWithAssumedDecimalPoint => "Parsing a float field failed (special TLE float representation with assumed decimal point)",
                ErrorWhat::ExpectedInteger => "Parsing an integer field failed",
                ErrorWhat::ExpectedSpace => "Found a non-space character between fields",
                ErrorWhat::ExpectedString => "Parsing a string field failed",
                ErrorWhat::FloatWithAssumedDecimalPointTooLong => "Tried to parse a float (special TLE float representation with assumed decimal point) with more than 16 ASCII characters",
                ErrorWhat::NoradIdMismatch => "NORAD mismatch between TLE lines",
                ErrorWhat::UnknownClassification => "Unknown classification code",
                ErrorWhat::FromYoOptFailed => "Date generation failed due to an error in the year",
                ErrorWhat::FromNumSecondsFromMidnightFailed => "Date generation failed due to an error in the seconds from midnight",
            },
            match self.line {
                ErrorLine::Line1 => "on TLE line 1",
                ErrorLine::Line2 => "on TLE line 2",
                ErrorLine::Both => "(TLE lines mismatch)",
            },
            self.start,
            self.end,
        ))
    }
}

#[cfg(feature = "std")]
impl std::error::Error for Error {}

trait TrimStart {
    fn trim_ascii_start_polyfill(&self) -> &[u8];
}

impl TrimStart for [u8] {
    fn trim_ascii_start_polyfill(&self) -> &[u8] {
        let mut bytes = self;
        while let [first, rest @ ..] = bytes {
            if first.is_ascii_whitespace() {
                bytes = rest;
            } else {
                break;
            }
        }
        bytes
    }
}

trait FromU8: Sized {
    type Err;
    fn from_u8(s: &[u8]) -> core::result::Result<Self, Self::Err>;
}

impl FromU8 for i8 {
    type Err = core::num::ParseIntError;
    fn from_u8(s: &[u8]) -> core::result::Result<Self, Self::Err> {
        // unsafe: parse calls from_str_radix which converts back to bytes.
        unsafe { core::str::from_utf8_unchecked(s) }.parse()
    }
}

impl FromU8 for u8 {
    type Err = core::num::ParseIntError;
    fn from_u8(s: &[u8]) -> core::result::Result<Self, Self::Err> {
        // unsafe: parse calls from_str_radix which converts back to bytes.
        unsafe { core::str::from_utf8_unchecked(s) }.parse()
    }
}

impl FromU8 for u64 {
    type Err = core::num::ParseIntError;
    fn from_u8(s: &[u8]) -> core::result::Result<Self, Self::Err> {
        // unsafe: parse calls from_str_radix which converts back to bytes.
        unsafe { core::str::from_utf8_unchecked(s) }.parse()
    }
}
impl FromU8 for f64 {
    type Err = core::num::ParseFloatError;
    fn from_u8(s: &[u8]) -> core::result::Result<Self, Self::Err> {
        // unsafe: parse calls dec2flt which converts back to bytes.
        unsafe { core::str::from_utf8_unchecked(s) }.parse()
    }
}

trait Parse {
    fn parse<F: FromU8>(&self) -> core::result::Result<F, F::Err>;
}

impl Parse for [u8] {
    fn parse<F: FromU8>(&self) -> core::result::Result<F, F::Err> {
        FromU8::from_u8(self)
    }
}

trait DecimalPointAssumedRepresentation {
    fn parse_decimal_point_assumed(
        &self,
        line: ErrorLine,
        start: usize,
        end: usize,
    ) -> core::result::Result<f64, Error>;
}

impl DecimalPointAssumedRepresentation for [u8] {
    fn parse_decimal_point_assumed(
        &self,
        line: ErrorLine,
        start: usize,
        end: usize,
    ) -> core::result::Result<f64, Error> {
        let trimmed = self.trim_ascii_start_polyfill();
        if trimmed.is_empty() {
            return Err(Error {
                what: ErrorWhat::ExpectedFloatWithAssumedDecimalPoint,
                line,
                start,
                end,
            });
        }
        let mut raw_buffer = [0_u8; 16];
        let length;
        if trimmed[0] == b'-' {
            if trimmed.len() < 2 || trimmed.len() + 1 > raw_buffer.len() {
                return Err(Error {
                    what: ErrorWhat::FloatWithAssumedDecimalPointTooLong,
                    line,
                    start,
                    end,
                });
            }
            raw_buffer[0] = b'-';
            raw_buffer[1] = b'.';
            raw_buffer[2..trimmed.len() + 1].copy_from_slice(&trimmed[1..]);
            length = trimmed.len() + 1;
        } else if trimmed[0] == b'+' {
            if trimmed.len() < 2 || trimmed.len() > raw_buffer.len() {
                return Err(Error {
                    what: ErrorWhat::FloatWithAssumedDecimalPointTooLong,
                    line,
                    start,
                    end,
                });
            }
            raw_buffer[0] = b'.';
            raw_buffer[1..trimmed.len()].copy_from_slice(&trimmed[1..]);
            length = trimmed.len();
        } else {
            if trimmed.len() + 1 > raw_buffer.len() {
                return Err(Error {
                    what: ErrorWhat::FloatWithAssumedDecimalPointTooLong,
                    line,
                    start,
                    end,
                });
            }
            raw_buffer[0] = b'.';
            raw_buffer[1..trimmed.len() + 1].copy_from_slice(trimmed);
            length = trimmed.len() + 1;
        }
        raw_buffer[0..length].parse::<f64>().map_err(|_| Error {
            what: ErrorWhat::ExpectedFloatWithAssumedDecimalPoint,
            line,
            start,
            end,
        })
    }
}

/// A satellite's elements classification
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Classification {
    /// Declassfied objects or objects without a classification
    #[cfg_attr(feature = "serde", serde(rename = "U"))]
    Unclassified,

    /// Would cause "serious damage" to national security if it were publicly available
    #[cfg_attr(feature = "serde", serde(rename = "C"))]
    Classified,

    /// Would cause "damage" or be prejudicial to national security if publicly available
    #[cfg_attr(feature = "serde", serde(rename = "S"))]
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
/// # fn main() -> anyhow::Result<()> {
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Elements {
    /// The name associated with the satellite
    #[cfg_attr(
        all(feature = "alloc", feature = "serde"),
        serde(rename = "OBJECT_NAME")
    )]
    #[cfg(feature = "alloc")]
    #[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
    pub object_name: Option<alloc::string::String>,

    /// The satellite's international designator
    ///
    /// It consists of the launch year, the launch number of that year and
    /// a letter code representing the sequential identifier of a piece in a launch.
    #[cfg_attr(all(feature = "alloc", feature = "serde"), serde(rename = "OBJECT_ID"))]
    #[cfg(feature = "alloc")]
    #[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
    pub international_designator: Option<alloc::string::String>,

    /// The catalog number USSPACECOM has designated for this object
    #[cfg_attr(
        feature = "serde",
        serde(rename = "NORAD_CAT_ID", deserialize_with = "u64_or_string")
    )]
    pub norad_id: u64,

    /// The elements' classification
    #[cfg_attr(feature = "serde", serde(rename = "CLASSIFICATION_TYPE"))]
    pub classification: Classification,

    /// The UTC timestamp of the elements
    #[cfg_attr(feature = "serde", serde(rename = "EPOCH"))]
    pub datetime: chrono::NaiveDateTime,

    /// Time derivative of the mean motion
    #[cfg_attr(
        feature = "serde",
        serde(rename = "MEAN_MOTION_DOT", deserialize_with = "f64_or_string")
    )]
    pub mean_motion_dot: f64,

    /// Second time derivative of the mean motion
    #[cfg_attr(
        feature = "serde",
        serde(rename = "MEAN_MOTION_DDOT", deserialize_with = "f64_or_string")
    )]
    pub mean_motion_ddot: f64,

    /// Radiation pressure coefficient in earth radii⁻¹
    #[cfg_attr(
        feature = "serde",
        serde(rename = "BSTAR", deserialize_with = "f64_or_string")
    )]
    pub drag_term: f64,

    /// A running count of all 2 line element sets generated by USSPACECOM for this object
    #[cfg_attr(
        feature = "serde",
        serde(rename = "ELEMENT_SET_NO", deserialize_with = "u64_or_string")
    )]
    pub element_set_number: u64,

    /// Angle between the equator and the orbit plane in deg
    #[cfg_attr(
        feature = "serde",
        serde(rename = "INCLINATION", deserialize_with = "f64_or_string")
    )]
    pub inclination: f64,

    /// Angle between vernal equinox and the point where the orbit crosses the equatorial plane in deg
    #[cfg_attr(
        feature = "serde",
        serde(rename = "RA_OF_ASC_NODE", deserialize_with = "f64_or_string")
    )]
    pub right_ascension: f64,

    /// The shape of the orbit
    #[cfg_attr(
        feature = "serde",
        serde(rename = "ECCENTRICITY", deserialize_with = "f64_or_string")
    )]
    pub eccentricity: f64,

    /// Angle between the ascending node and the orbit's point of closest approach to the earth in deg
    #[cfg_attr(
        feature = "serde",
        serde(rename = "ARG_OF_PERICENTER", deserialize_with = "f64_or_string")
    )]
    pub argument_of_perigee: f64,

    /// Angle of the satellite location measured from perigee in deg
    #[cfg_attr(
        feature = "serde",
        serde(rename = "MEAN_ANOMALY", deserialize_with = "f64_or_string")
    )]
    pub mean_anomaly: f64,

    /// Mean number of orbits per day in day⁻¹ (Kozai convention)
    #[cfg_attr(
        feature = "serde",
        serde(rename = "MEAN_MOTION", deserialize_with = "f64_or_string")
    )]
    pub mean_motion: f64,

    /// The orbit number at epoch
    #[cfg_attr(
        feature = "serde",
        serde(rename = "REV_AT_EPOCH", deserialize_with = "u64_or_string")
    )]
    pub revolution_number: u64,

    /// NORAD internal use, always 0 in distributed data
    #[cfg_attr(
        feature = "serde",
        serde(rename = "EPHEMERIS_TYPE", deserialize_with = "u8_or_string")
    )]
    pub ephemeris_type: u8,
}

#[cfg(feature = "serde")]
fn u64_or_string<'de, D>(deserializer: D) -> core::result::Result<u64, D::Error>
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

#[cfg(feature = "serde")]
fn u8_or_string<'de, D>(deserializer: D) -> core::result::Result<u8, D::Error>
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

#[cfg(feature = "serde")]
fn f64_or_string<'de, D>(deserializer: D) -> core::result::Result<f64, D::Error>
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

/// Returns the number of years since UTC 1 January 2000 12h00 (J2000)
///
/// This is the recommended method to calculate the epoch
pub fn julian_years_since_j2000(datetime: &chrono::NaiveDateTime) -> f64 {
    // y₂₀₀₀ = (367 yᵤ - ⌊7 (yᵤ + ⌊(mᵤ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤ / 9⌋ + dᵤ - 730531) / 365.25
    //         + (3600 hᵤ + 60 minᵤ + sᵤ - 43200) / (24 × 60 × 60 × 365.25)
    //         + nsᵤ / (24 × 60 × 60 × 365.25 × 10⁹)
    (367 * datetime.year() - (7 * (datetime.year() + (datetime.month() as i32 + 9) / 12)) / 4
        + 275 * datetime.month() as i32 / 9
        + datetime.day() as i32
        - 730531) as f64
        / 365.25
        + (datetime.num_seconds_from_midnight() as i32 - 43200) as f64
            / (24.0 * 60.0 * 60.0 * 365.25)
        + (datetime.nanosecond() as f64) / (24.0 * 60.0 * 60.0 * 1e9 * 365.25)
}

/// Returns the number of years since UTC 1 January 2000 12h00 (J2000) using the AFSPC expression
///
/// This function should be used if compatibility with the AFSPC implementation is needed
pub fn julian_years_since_j2000_afspc_compatibility_mode(datetime: &chrono::NaiveDateTime) -> f64 {
    // y₂₀₀₀ = (367 yᵤ - ⌊7 (yᵤ + ⌊(mᵤ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤ / 9⌋ + dᵤ
    //         + 1721013.5
    //         + (((nsᵤ / 10⁹ + sᵤ) / 60 + minᵤ) / 60 + hᵤ) / 24
    //         - 2451545)
    //         / 365.25
    ((367 * datetime.year() as u32
        - (7 * (datetime.year() as u32 + (datetime.month() + 9) / 12)) / 4
        + 275 * datetime.month() / 9
        + datetime.day()) as f64
        + 1721013.5
        + (((datetime.nanosecond() as f64 / 1e9 + datetime.second() as f64) / 60.0
            + datetime.minute() as f64)
            / 60.0
            + datetime.hour() as f64)
            / 24.0
        - 2451545.0)
        / 365.25
}

/// Minutes ellapsed since the elements' epoch
///
/// This number can be negative since SGP4 can propagate back in time.
#[derive(Debug, Clone, Copy)]
pub struct MinutesSinceEpoch(pub f64);

/// Nanoseconds overflow while converting from datetime to minutes since epoch
///
/// 2⁶⁴ nanoseconds correspond to about 585 years.
/// Overflows are almost certainly caused by data corruption or code bugs.
#[derive(Debug, Clone)]
pub struct DatetimeToMinutesSinceEpochError {
    from: chrono::NaiveDateTime,
    to: chrono::NaiveDateTime,
}

impl core::fmt::Display for DatetimeToMinutesSinceEpochError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        formatter.write_fmt(format_args!(
            "Nanoseconds overflow when calculating {} - {}",
            self.to, self.from
        ))
    }
}

#[cfg(feature = "std")]
impl std::error::Error for DatetimeToMinutesSinceEpochError {}

/// Nanoseconds overflow while converting from minutes since epoch to datetime
///
/// 2⁶⁴ nanoseconds correspond to about 585 years.
/// Overflows are almost certainly caused by data corruption or code bugs.
#[derive(Debug, Clone)]
pub enum MinutesSinceEpochToDatetimeError {
    MinutesToNanoseconds(f64),
    Add {
        datetime: chrono::NaiveDateTime,
        duration: chrono::Duration,
    },
}

impl core::fmt::Display for MinutesSinceEpochToDatetimeError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MinutesSinceEpochToDatetimeError::MinutesToNanoseconds(nanoseconds) => formatter
                .write_fmt(format_args!(
                    "Nanoseconds overflow when calculating {nanoseconds} * 60e9"
                )),
            MinutesSinceEpochToDatetimeError::Add { datetime, duration } => formatter.write_fmt(
                format_args!("Nanoseconds overflow when calculating {datetime} + {duration}"),
            ),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for MinutesSinceEpochToDatetimeError {}

impl Elements {
    fn from_lines(line1: &[u8], line2: &[u8]) -> core::result::Result<Elements, Error> {
        if line1.len() != 69 {
            return Err(Error {
                what: ErrorWhat::BadLength,
                line: ErrorLine::Line1,
                start: 0,
                end: line1.len(),
            });
        }
        if line2.len() != 69 {
            return Err(Error {
                what: ErrorWhat::BadLength,
                line: ErrorLine::Line2,
                start: 0,
                end: line2.len(),
            });
        }
        if line1[0] != b'1' {
            return Err(Error {
                what: ErrorWhat::BadFirstCharacter,
                line: ErrorLine::Line1,
                start: 0,
                end: 1,
            });
        }
        if line2[0] != b'2' {
            return Err(Error {
                what: ErrorWhat::BadFirstCharacter,
                line: ErrorLine::Line2,
                start: 0,
                end: 1,
            });
        }
        for index in [1, 8, 17, 32, 43, 52, 61, 63].iter() {
            if line1[*index] != b' ' {
                return Err(Error {
                    what: ErrorWhat::ExpectedSpace,
                    line: ErrorLine::Line1,
                    start: *index,
                    end: *index + 1,
                });
            }
        }
        for index in [1, 7, 16, 25, 33, 42, 51].iter() {
            if line2[*index] != b' ' {
                return Err(Error {
                    what: ErrorWhat::ExpectedSpace,
                    line: ErrorLine::Line2,
                    start: *index,
                    end: *index + 1,
                });
            }
        }
        let norad_id = line1[2..7]
            .trim_ascii_start_polyfill()
            .parse::<u64>()
            .map_err(|_| Error {
                what: ErrorWhat::ExpectedInteger,
                line: ErrorLine::Line1,
                start: 2,
                end: 7,
            })?;
        if norad_id
            != line2[2..7]
                .trim_ascii_start_polyfill()
                .parse::<u64>()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedInteger,
                    line: ErrorLine::Line2,
                    start: 2,
                    end: 7,
                })?
        {
            return Err(Error {
                what: ErrorWhat::NoradIdMismatch,
                line: ErrorLine::Both,
                start: 2,
                end: 7,
            });
        }
        for (line, content) in [(ErrorLine::Line1, &line1), (ErrorLine::Line2, &line2)] {
            if (content[..68]
                .iter()
                .fold(0, |accumulator, character| match character {
                    b'-' => accumulator + 1,
                    character if (&b'0'..=&b'9').contains(&character) => {
                        accumulator + (character - b'0') as u16
                    }
                    _ => accumulator,
                })
                % 10) as u8
                != content[68] - b'0'
            {
                return Err(Error {
                    what: ErrorWhat::BadChecksum,
                    line,
                    start: 68,
                    end: 69,
                });
            }
        }
        Ok(Elements {
            #[cfg(feature = "alloc")]
            object_name: None,
            norad_id,
            classification: match line1[7] {
                b'U' => Classification::Unclassified,
                b'C' => Classification::Classified,
                b'S' => Classification::Secret,
                _ => {
                    return Err(Error {
                        what: ErrorWhat::UnknownClassification,
                        line: ErrorLine::Line1,
                        start: 7,
                        end: 8,
                    })
                }
            },
            #[cfg(feature = "alloc")]
            international_designator: if line1[9..17]
                .iter()
                .all(|character| character.is_ascii_whitespace())
            {
                None
            } else {
                Some(format!(
                    "{}-{}",
                    match line1[9..11].parse::<u8>().map_err(|_| Error {
                        what: ErrorWhat::ExpectedInteger,
                        line: ErrorLine::Line1,
                        start: 9,
                        end: 11,
                    })? {
                        launch_year if launch_year < 57 => 2000 + launch_year as u16,
                        launch_year => 1900 + launch_year as u16,
                    },
                    core::str::from_utf8(&line1[11..17])
                        .map_err(|_| Error {
                            what: ErrorWhat::ExpectedString,
                            line: ErrorLine::Line1,
                            start: 11,
                            end: 17,
                        })?
                        .trim()
                ))
            },
            datetime: {
                let day = line1[20..32]
                    .trim_ascii_start_polyfill()
                    .parse::<f64>()
                    .map_err(|_| Error {
                        what: ErrorWhat::ExpectedFloat,
                        line: ErrorLine::Line1,
                        start: 20,
                        end: 32,
                    })?;
                let seconds = day.fract() * (24.0 * 60.0 * 60.0);
                chrono::NaiveDate::from_yo_opt(
                    match line1[18..20].parse::<u8>().map_err(|_| Error {
                        what: ErrorWhat::ExpectedFloat,
                        line: ErrorLine::Line1,
                        start: 18,
                        end: 20,
                    })? {
                        year if year < 57 => year as i32 + 2000,
                        year => year as i32 + 1900,
                    },
                    day as u32,
                )
                .ok_or(Error {
                    what: ErrorWhat::FromYoOptFailed,
                    line: ErrorLine::Line1,
                    start: 18,
                    end: 20,
                })?
                .and_time(
                    chrono::NaiveTime::from_num_seconds_from_midnight_opt(
                        seconds as u32,
                        (seconds.fract() * 1e9).round() as u32,
                    )
                    .ok_or(Error {
                        what: ErrorWhat::FromNumSecondsFromMidnightFailed,
                        line: ErrorLine::Line1,
                        start: 20,
                        end: 32,
                    })?,
                )
            },
            mean_motion_dot: line1[33..43]
                .trim_ascii_start_polyfill()
                .parse()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedFloat,
                    line: ErrorLine::Line1,
                    start: 33,
                    end: 43,
                })?,
            mean_motion_ddot: line1[44..50].parse_decimal_point_assumed(
                ErrorLine::Line1,
                44,
                50,
            )? * 10.0_f64.powi(line1[50..52].parse::<i8>().map_err(|_| Error {
                what: ErrorWhat::ExpectedFloat,
                line: ErrorLine::Line1,
                start: 50,
                end: 52,
            })? as i32),
            drag_term: line1[53..59].parse_decimal_point_assumed(ErrorLine::Line1, 53, 59)?
                * 10.0_f64.powi(line1[59..61].parse::<i8>().map_err(|_| Error {
                    what: ErrorWhat::ExpectedFloat,
                    line: ErrorLine::Line1,
                    start: 59,
                    end: 61,
                })? as i32),
            ephemeris_type: line1[62..63]
                .trim_ascii_start_polyfill()
                .parse()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedInteger,
                    line: ErrorLine::Line1,
                    start: 62,
                    end: 63,
                })?,
            element_set_number: line1[64..68].trim_ascii_start_polyfill().parse().map_err(
                |_| Error {
                    what: ErrorWhat::ExpectedInteger,
                    line: ErrorLine::Line1,
                    start: 64,
                    end: 68,
                },
            )?,
            inclination: line2[8..16]
                .trim_ascii_start_polyfill()
                .parse()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedFloat,
                    line: ErrorLine::Line2,
                    start: 8,
                    end: 16,
                })?,
            right_ascension: line2[17..25]
                .trim_ascii_start_polyfill()
                .parse()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedFloat,
                    line: ErrorLine::Line2,
                    start: 17,
                    end: 25,
                })?,
            eccentricity: line2[26..33].parse_decimal_point_assumed(ErrorLine::Line2, 26, 33)?,
            argument_of_perigee: line2[34..42].trim_ascii_start_polyfill().parse().map_err(
                |_| Error {
                    what: ErrorWhat::ExpectedFloat,
                    line: ErrorLine::Line2,
                    start: 34,
                    end: 42,
                },
            )?,
            mean_anomaly: line2[43..51]
                .trim_ascii_start_polyfill()
                .parse()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedFloat,
                    line: ErrorLine::Line2,
                    start: 43,
                    end: 51,
                })?,
            mean_motion: line2[52..63]
                .trim_ascii_start_polyfill()
                .parse()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedFloat,
                    line: ErrorLine::Line2,
                    start: 52,
                    end: 63,
                })?,
            revolution_number: line2[63..68]
                .trim_ascii_start_polyfill()
                .parse()
                .map_err(|_| Error {
                    what: ErrorWhat::ExpectedInteger,
                    line: ErrorLine::Line2,
                    start: 63,
                    end: 68,
                })?,
        })
    }

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
    /// # fn main() -> anyhow::Result<()> {
    /// let elements = sgp4::Elements::from_tle(
    ///     Some("ISS (ZARYA)".to_owned()),
    ///     "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927".as_bytes(),
    ///     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537".as_bytes(),
    /// )?;
    /// #     Ok(())
    /// # }
    /// ```
    #[cfg(feature = "alloc")]
    pub fn from_tle(
        object_name: Option<alloc::string::String>,
        line1: &[u8],
        line2: &[u8],
    ) -> core::result::Result<Elements, Error> {
        let mut result = Self::from_lines(line1, line2)?;
        result.object_name = object_name;
        Ok(result)
    }

    #[cfg(not(feature = "alloc"))]
    pub fn from_tle(line1: &[u8], line2: &[u8]) -> core::result::Result<Elements, Error> {
        Self::from_lines(line1, line2)
    }

    /// Returns the number of years since UTC 1 January 2000 12h00 (J2000)
    ///
    /// This is the recommended method to calculate the epoch
    pub fn epoch(&self) -> f64 {
        julian_years_since_j2000(&self.datetime)
    }

    /// Returns the number of years since UTC 1 January 2000 12h00 (J2000) using the AFSPC expression
    ///
    /// This function should be used if compatibility with the AFSPC implementation is needed
    pub fn epoch_afspc_compatibility_mode(&self) -> f64 {
        julian_years_since_j2000_afspc_compatibility_mode(&self.datetime)
    }

    /// Returns the time difference in minutes between the given datetime and the elements' epoch
    ///
    /// This method does not take leap seconds into account
    pub fn datetime_to_minutes_since_epoch(
        &self,
        datetime: &chrono::NaiveDateTime,
    ) -> core::result::Result<MinutesSinceEpoch, DatetimeToMinutesSinceEpochError> {
        (*datetime - self.datetime)
            .num_nanoseconds()
            .ok_or(DatetimeToMinutesSinceEpochError {
                from: self.datetime,
                to: *datetime,
            })
            .map(|nanoseconds| MinutesSinceEpoch(nanoseconds as f64 / 60e9))
    }

    /// Builds a datetime from a number of minutes since epoch
    ///
    /// This method does not take leap seconds into account
    pub fn minutes_since_epoch_to_datetime(
        &self,
        minutes_since_epoch: &MinutesSinceEpoch,
    ) -> core::result::Result<chrono::NaiveDateTime, MinutesSinceEpochToDatetimeError> {
        let nanoseconds = minutes_since_epoch.0 * 60e9;
        if nanoseconds > i64::MAX as f64 || nanoseconds < i64::MIN as f64 {
            Err(MinutesSinceEpochToDatetimeError::MinutesToNanoseconds(
                minutes_since_epoch.0,
            ))
        } else {
            let duration = chrono::Duration::nanoseconds(nanoseconds.round() as i64);
            self.datetime.checked_add_signed(duration).ok_or(
                MinutesSinceEpochToDatetimeError::Add {
                    datetime: self.datetime,
                    duration,
                },
            )
        }
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
#[cfg(feature = "alloc")]
#[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
pub fn parse_2les(tles: &str) -> core::result::Result<alloc::vec::Vec<Elements>, Error> {
    let mut line_buffer = "";
    let mut first = true;
    let mut elements_vec = alloc::vec::Vec::new();
    for line in tles.lines() {
        if first {
            line_buffer = line;
        } else {
            elements_vec.push(Elements::from_tle(
                None,
                line_buffer.as_bytes(),
                line.as_bytes(),
            )?);
        }
        first = !first;
    }
    Ok(elements_vec)
}

/// Parses a multi-line TL/3LE string into a list of `Elements`
///
/// Each triplet of lines must represent a TLE with an object name, for example as in
/// [https://celestrak.com/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle](https://celestrak.com/NORAD/elements/gp.php?GROUP=stations&FORMAT=tle).
///
/// # Arguments
///
/// * `tles` - A string containing multiple lines
#[cfg(feature = "alloc")]
#[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
pub fn parse_3les(tles: &str) -> core::result::Result<alloc::vec::Vec<Elements>, Error> {
    let mut lines_buffer = ["", ""];
    let mut index = 0;
    let mut elements_vec = alloc::vec::Vec::new();
    for line in tles.lines() {
        match index {
            0 | 1 => {
                lines_buffer[index] = line;
                index += 1;
            }
            _ => {
                elements_vec.push(Elements::from_tle(
                    Some(lines_buffer[0].to_owned()),
                    lines_buffer[1].as_bytes(),
                    line.as_bytes(),
                )?);
                index = 0;
            }
        }
    }
    Ok(elements_vec)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "serde")]
    fn assert_eq_f64(first: f64, second: f64) {
        if second == 0.0 {
            assert_eq!(first, 0.0);
        } else {
            assert!((first - second).abs() / second < f64::EPSILON);
        }
    }

    #[cfg(feature = "serde")]
    #[test]
    fn test_from_celestrak_omm() -> anyhow::Result<()> {
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
        )
        .map_err(|error| anyhow::anyhow!("{error}"))?;
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
            chrono::NaiveDate::from_yo_opt(2020, 194).unwrap().and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight_opt(4747, 402656000).unwrap()
            )
        );
        assert_eq_f64(elements.epoch(), 20.527_186_712_635_18);
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
    #[cfg(feature = "alloc")]
    fn test_from_space_track_omm() -> anyhow::Result<()> {
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
        )
        .map_err(|error| anyhow::anyhow!("{error}"))?;
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
            chrono::NaiveDate::from_yo_opt(2020, 348).unwrap().and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight_opt(59764, 502592000).unwrap()
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
    #[cfg(feature = "alloc")]
    fn test_from_celestrak_omms() -> anyhow::Result<()> {
        let elements_vec: [Elements; 2] = serde_json::from_str(
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
        )
        .map_err(|error| anyhow::anyhow!("{error}"))?;
        assert_eq!(elements_vec.len(), 2);
        Ok(())
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn test_from_tle() -> core::result::Result<(), Error> {
        let elements = Elements::from_tle(
            Some("ISS (ZARYA)".into()),
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
            chrono::NaiveDate::from_yo_opt(2008, 264).unwrap().and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight_opt(44740, 104192001).unwrap()
            )
        );
        assert_eq_f64(elements.epoch(), 8.720103559972621);
        assert_eq_f64(
            elements.epoch_afspc_compatibility_mode(),
            8.720_103_559_972_213,
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
            chrono::NaiveDate::from_yo_opt(1980, 230).unwrap().and_time(
                chrono::NaiveTime::from_num_seconds_from_midnight_opt(25600, 136832000).unwrap()
            )
        );
        assert_eq_f64(elements.epoch(), -19.373_589_875_756_33);
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
    #[cfg(feature = "alloc")]
    fn test_parse_2les() -> core::result::Result<(), Error> {
        let elements_vec = parse_2les(
            "1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992\n\
             2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008\n\
             1 42982U 98067NE  20194.06866787  .00008489  00000-0  72204-4 0  9997\n\
             2 42982  51.6338 155.6245 0002758 166.8841 193.2228 15.70564504154944\n",
        )?;
        assert_eq!(elements_vec.len(), 2);
        Ok(())
    }

    #[test]
    #[cfg(feature = "alloc")]
    fn test_parse_3les() -> core::result::Result<(), Error> {
        let elements_vec = parse_3les(
            "ISS (ZARYA)\n\
             1 25544U 98067A   20194.88612269 -.00002218  00000-0 -31515-4 0  9992\n\
             2 25544  51.6461 221.2784 0001413  89.1723 280.4612 15.49507896236008\n\
             KESTREL EYE IIM (KE2M)\n\
             1 42982U 98067NE  20194.06866787  .00008489  00000-0  72204-4 0  9997\n\
             2 42982  51.6338 155.6245 0002758 166.8841 193.2228 15.70564504154944\n",
        )?;
        assert_eq!(elements_vec.len(), 2);
        Ok(())
    }
}
