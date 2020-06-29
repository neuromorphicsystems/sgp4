const MONTHS_ENDS: [u16; 12] = [32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366];
const LEAP_YEAR_MONTHS_ENDS: [u16; 12] = [32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367];

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

pub enum Classification {
    Unclassified,
    Classified,
    Secret,
}

// https://spaceflight.nasa.gov/realdata/sightings/SSapplications/Post/JavaSSOP/SSOP_Help/tle_def.html
// epoch: years since UTC 1 January 2000, 0h00
// drag_term: the radiation pressure coefficient B*, in earth radii⁻¹
// inclination: the angle between the equator and the orbit plane iₒ, in degrees
// right_ascension: the angle between vernal equinox and
//     the point where the orbit crosses the equatorial plane, in degrees
// eccentricity: the shape of the orbit eₒ
// argument_of_perigee: the angle between the ascending node and
//     the orbit's point of closest approach to the earth ωₒ, in degrees
// mean_anomaly: the angle of the satellite location measured from perigee Mₒ, in degrees
// mean_motion: mean number of orbits per day nₒ, in days⁻¹
pub struct Tle {
    pub satellite_name: Option<String>,
    pub satellite_number: u32,
    pub classification: Classification,
    pub launch_year: Option<u16>,
    pub launch_number: Option<u16>,
    pub launch_piece: Option<String>,
    pub year: u16,
    pub day: f64,
    pub mean_motion_first_derivative: f64,
    pub mean_motion_second_derivative: f64,
    pub drag_term: f64,
    pub element_set_number: u16,
    pub inclination: f64,
    pub right_ascension: f64,
    pub eccentricity: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
    pub mean_motion: f64,
    pub revolution_number: u32,
}

pub fn parse(title: Option<String>, line1: &[u8], line2: &[u8]) -> Result<Tle> {
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
    let satellite_number = std::str::from_utf8(&line1[2..7])?.parse::<u32>()?;
    if satellite_number != std::str::from_utf8(&line2[2..7])?.parse()? {
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
    Ok(Tle {
        satellite_name: title,
        satellite_number: satellite_number,
        classification: match line1[7] {
            b'U' => Classification::Unclassified,
            b'C' => Classification::Classified,
            b'S' => Classification::Secret,
            _ => return Err(Error::new("unknown classification".to_owned())),
        },
        launch_year: if line1[9..11].iter().all(|character| *character == ' ' as u8) {
            None
        } else {
            Some(match std::str::from_utf8(&line1[9..11])?.parse::<u8>()? {
                launch_year if launch_year < 57 => 2000 + launch_year as u16,
                launch_year => 1900 + launch_year as u16,
            })
        },
        launch_number: if line1[11..14]
            .iter()
            .all(|character| *character == ' ' as u8)
        {
            None
        } else {
            Some(std::str::from_utf8(&line1[11..14])?.parse()?)
        },
        launch_piece: if line1[14..15]
            .iter()
            .all(|character| *character == ' ' as u8)
        {
            None
        } else {
            Some(std::str::from_utf8(&line1[14..15])?.trim_end().to_owned())
        },
        year: match std::str::from_utf8(&line1[18..20])?.parse::<u8>()? {
            year if year < 57 => year as u16 + 2000,
            year => year as u16 + 1900,
        },
        day: std::str::from_utf8(&line1[20..32])?
            .trim_start()
            .parse::<f64>()?,
        mean_motion_first_derivative: std::str::from_utf8(&line1[33..43])?.trim_start().parse()?,
        mean_motion_second_derivative: line1[44..50].parse_decimal_point_assumed()?
            * 10.0_f64.powi(std::str::from_utf8(&line1[50..52])?.parse::<i8>()? as i32),
        drag_term: line1[53..59].parse_decimal_point_assumed()?
            * 10.0_f64.powi(std::str::from_utf8(&line1[59..61])?.parse::<i8>()? as i32),
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

impl Tle {
    pub fn international_designator(&self) -> Option<String> {
        if self.launch_year.is_none() || self.launch_number.is_none() || self.launch_piece.is_none()
        {
            None
        } else {
            Some(format!(
                "{}-{:03}{}",
                self.launch_year.unwrap(),
                self.launch_number.unwrap(),
                self.launch_piece.as_ref().unwrap()
            ))
        }
    }

    pub fn epoch(&self) -> f64 {
        let months_ends = if self.year % 4 == 0 {
            &LEAP_YEAR_MONTHS_ENDS
        } else {
            &MONTHS_ENDS
        };
        let month = months_ends
            .iter()
            .position(|month_end| *month_end > self.day as u16)
            .unwrap_or(11) as u32
            + 1;
        ((367 * self.year as u32 - (7 * (self.year as u32 + (month + 9) / 12)) / 4
            + 275 * month / 9) as f64
            - 730531.5
            + (self.day
                - if month == 1 {
                    0
                } else {
                    months_ends[(month - 2) as usize] - 1
                } as f64))
            / 365.25
    }

    pub fn epoch_afspc_compatibility_mode(&self) -> f64 {
        let months_ends = if self.year % 4 == 0 {
            &LEAP_YEAR_MONTHS_ENDS
        } else {
            &MONTHS_ENDS
        };
        let month = months_ends
            .iter()
            .position(|month_end| *month_end > self.day as u16)
            .unwrap_or(11) as u32
            + 1;
        ((367 * self.year as u32 - (7 * (self.year as u32 + (month + 9) / 12)) / 4
            + 275 * month / 9) as f64
            + 1721013.5
            + (self.day
                - if month == 1 {
                    0
                } else {
                    months_ends[(month - 2) as usize] - 1
                } as f64)
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
    fn test_parse() -> Result<()> {
        let tle = parse(
            Some("ISS (ZARYA)".to_owned()),
            "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927".as_bytes(),
            "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537".as_bytes(),
        )?;
        match tle.satellite_name.as_ref() {
            Some(satellite_name) => assert_eq!(satellite_name, "ISS (ZARYA)"),
            None => panic!(),
        }
        assert_eq!(tle.satellite_number, 25544);
        assert!(matches!(tle.classification, Classification::Unclassified));
        assert_eq!(tle.launch_year.unwrap(), 1998);
        assert_eq!(tle.launch_number.unwrap(), 67);
        assert_eq!(tle.launch_piece.as_ref().unwrap(), "A");
        assert_eq!(tle.international_designator().unwrap(), "1998-067A");
        assert_eq!(tle.year, 2008);
        assert_eq_f64(tle.day, 264.51782528);
        assert_eq_f64(tle.mean_motion_first_derivative, -0.00002182);
        assert_eq_f64(tle.mean_motion_second_derivative, 0.0);
        assert_eq_f64(tle.drag_term, -0.11606e-4);
        assert_eq!(tle.element_set_number, 292);
        assert_eq_f64(tle.inclination, 51.6416);
        assert_eq_f64(tle.right_ascension, 247.4627);
        assert_eq_f64(tle.eccentricity, 0.0006703);
        assert_eq_f64(tle.argument_of_perigee, 130.5360);
        assert_eq_f64(tle.mean_anomaly, 325.0288);
        assert_eq_f64(tle.mean_motion, 15.72125391);
        assert_eq!(tle.revolution_number, 56353);
        let tle = parse(
            None,
            "1 11801U          80230.29629788  .01431103  00000-0  14311-1 0    13".as_bytes(),
            "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848    13".as_bytes(),
        )?;
        assert!(tle.satellite_name.is_none());
        assert_eq!(tle.satellite_number, 11801);
        assert!(matches!(tle.classification, Classification::Unclassified));
        assert!(tle.launch_year.is_none());
        assert!(tle.launch_number.is_none());
        assert!(tle.launch_piece.is_none());
        assert!(tle.international_designator().is_none());
        assert_eq!(tle.year, 1980);
        assert_eq_f64(tle.day, 230.29629788);
        assert_eq_f64(tle.mean_motion_first_derivative, 0.01431103);
        assert_eq_f64(tle.mean_motion_second_derivative, 0.0);
        assert_eq_f64(tle.drag_term, 0.014311);
        assert_eq!(tle.element_set_number, 1);
        assert_eq_f64(tle.inclination, 46.7916);
        assert_eq_f64(tle.right_ascension, 230.4354);
        assert_eq_f64(tle.eccentricity, 0.7318036);
        assert_eq_f64(tle.argument_of_perigee, 47.4722);
        assert_eq_f64(tle.mean_anomaly, 10.4117);
        assert_eq_f64(tle.mean_motion, 2.28537848);
        assert_eq!(tle.revolution_number, 1);
        Ok(())
    }
}
