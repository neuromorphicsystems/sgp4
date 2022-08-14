```toml
[dependencies]
sgp4 = "0.6"
```

The SGP4 algorithm, ported to Rust from the reference Celestrak implementation [[1]](#1).

The code was entirely refactored to leverage Rust's algebraic data types and highlight the relationship between the implementation and the reference mathematical equations [[2]](#2).

The numerical predictions are almost identical to those of the Celestrak implementation. The observed differences (less than 2 × 10⁻⁷ km for the position and 10⁻⁹ km.s⁻¹ for the velocity three and a half years after epoch) are well below the accuracy of the algorithm.

We drew inspiration from the incomplete https://github.com/natronics/rust-sgp4 to write mathematical expressions using UTF-8 characters.

- [Documentation](#documentation)
- [Benchmark](#benchmark)
- [Variables and mathematical expressions](#variables-and-mathematical-expressions)
    - [Variables](#variables)
        - [Initialization variables](#initialization-variables)
        - [Propagation variables](#propagation-variables)
        - [Third-body initialization variables](#third-body-initialization-variables)
        - [Third-body propagation variables](#third-body-propagation-variables)
    - [Mathematical expressions](#mathematical-expressions)
        - [UT1 to Julian conversion](#ut1-to-julian-conversion)
        - [Common initialization](#common-initialization)
        - [Near earth initialization](#near-earth-initialization)
        - [High altitude near earth initialization](#high-altitude-near-earth-initialization)
        - [Elliptic high altitude near earth initialization](#elliptic-high-altitude-near-earth-initialization)
        - [Deep space initialization](#deep-space-initialization)
        - [Third body perturbations](#third-body-perturbations)
        - [Resonant deep space initialization](#resonant-deep-space-initialization)
        - [Geosynchronous deep space initialization](#geosynchronous-deep-space-initialization)
        - [Molniya deep space initialization](#molniya-deep-space-initialization)
        - [Common propagation](#common-propagation)
        - [Near earth propagation](#near-earth-propagation)
        - [High altitude near earth propagation](#high-altitude-near-earth-propagation)
        - [Deep space propagation](#deep-space-propagation)
        - [Third body propagation](#third-body-propagation)
        - [Resonant deep space propagation](#resonant-deep-space-propagation)
        - [Lyddane deep space propagation](#lyddane-deep-space-propagation)
- [References](#references)

## Documentation

The code documentation is hosted at [https://docs.rs/sgp4/0.6.0/sgp4/](https://docs.rs/sgp4/0.6.0/sgp4/).

Examples can be found in this repository's _examples_ directory:

-   _examples/celestrak.rs_ retrieves the most recent Galileo OMMs from Celestrak and propagates them
-   _examples/omm.rs_ parses and propagates a JSON-encoded OMM
-   _examples/space-track.rs_ retrieves the 20 most recent launches OMMs from Space-Track and propagates them
-   _examples/tle.rs_ parses and propagates a TLE
-   _examples/tle_afspc.rs_ parses and propagates a TLE using the AFSPC compatibility mode
-   _examples/advanced.rs_ leverages the advanced API to (marginally) accelerate the propagation of deep space resonant satellites

To run an example (here _examples/celestrak.rs_), use:

```sh
cargo run --example celestrak
```

To run the Space-Track example, you must first assign your Space-Track.org credentials to the fields `identity` and `password` (see lines 3 and 4 in _examples/space-track.rs_).

## Benchmark

The benchmark code is available at https://github.com/neuromorphicsystems/sgp4-benchmark. It compares two SGP4 implementations in different configurations:

-   `cpp`: the Celestrak implementation [[1]](#1) in improved mode
-   `cpp-afspc`: the Celestrak implementation [[1]](#1) in AFSPC compatibility mode
-   `cpp-fastmath`: the Celestrak implementation [[1]](#1) in improved mode with the `fast-math` compiler flag
-   `cpp-afspc-fastmath`: the Celestrak implementation [[1]](#1) in AFSPC compatibility mode with the `fast-math` compiler flag
-   `rust`: our Rust implementation in default mode
-   `rust-afspc`: our Rust implementation in AFSPC compatibility mode

This benchmark must not be confused with the code in this repository's _bench_ directory. The latter considers only a small subset of the Celestrak catalogue (the tests recommended in [[1]](#1)) and does not measure the original C++ implementation.

The present results were obtained using a machine with the following configuration:

-   **CPU** - Intel Core i7-8700 @ 3.20GHz
-   **RAM** - Kingston DDR4 @ 2.667 GHz
-   **OS** - Ubuntu 16.04
-   **Compilers** - Rust 1.44.1 and gcc 9.3.0

Accuracy measures the maximum propagation error of each implementation with respect to the reference implementation (`cpp-afspc`) over the full Celestrak catalogue (1 minute timestep over 24 hours).

| implementation       | maximum position error | maximum speed error |
| -------------------- | ---------------------- | ------------------- |
| `cpp-afspc`          | (reference)            | (reference)         |
| `cpp`                | 1.05 km                | 1.30 × 10⁻³ km.s⁻¹  |
| `cpp-fastmath`       | 1.05 km                | 1.30 × 10⁻³ km.s⁻¹  |
| `cpp-afspc-fastmath` | 4.21 × 10⁻⁸ km         | 7.51 × 10⁻¹² km.s⁻¹ |
| `rust`               | 1.05 km                | 1.30 × 10⁻³ km.s⁻¹  |
| `rust-afspc`         | 4.19 × 10⁻⁸ km         | 7.46 × 10⁻¹² km.s⁻¹ |

The Rust and C++ fast-math errors have the same order of magnitude. In both cases, they can be attributed to mathematically identical expressions implemented with different floating-point operations.

Speed measures the time it takes to propagate every satellite in the Celestrak catalogue (1 minute timestep over 24 hours) using a single thread. 100 values are sampled per implementation.

| implementation       | minimum | Q1     | median | Q3     | maximum | relative difference |
| -------------------- | ------- | ------ | ------ | ------ | ------- | ------------------- |
| `cpp-afspc`          | 8.95 s  | 9.02 s | 9.03 s | 9.06 s | 9.18 s  | (reference)         |
| `cpp`                | 8.95 s  | 9.01 s | 9.04 s | 9.06 s | 9.25 s  | + 0 %               |
| `cpp-fastmath`       | 7.67 s  | 7.74 s | 7.77 s | 7.79 s | 7.90 s  | - 14 %              |
| `cpp-afspc-fastmath` | 7.70 s  | 7.74 s | 7.76 s | 7.79 s | 7.86 s  | - 14 %              |
| `rust`               | 8.36 s  | 8.41 s | 8.43 s | 8.45 s | 8.53 s  | - 7 %               |
| `rust-afspc`         | 8.36 s  | 8.41 s | 8.43 s | 8.46 s | 8.59 s  | - 7 %               |

Rust fast-math support is a work in progress (see https://github.com/rust-lang/rust/issues/21690). Similarly to C++, it should have a very small impact on accuracy while providing a substantial speed gain.

## Variables and mathematical expressions

### Variables

Each variable is used to store the result of one and only one expression. Most variables are immutable, with the exception of the variable `(E + ω)ᵢ` used to solve Kepler's equation and the state variables `tᵢ`, `nᵢ` and `λᵢ` used to integrate the resonance effects of Earth gravity.

The following tables list the variables used in the code and their associated mathematical symbol. Where possible, we used symbols from [[2]](#2). Partial expressions without a name in [[2]](#2) follow the convention `kₙ, n ∈ ℕ` if they are shared between initialization and propagation, and `pₙ, n ∈ ℕ` if they are local to initialization or propagation.

1. [Initialization variables](#initialization-variables)
2. [Propagation variables](#propagation-variables)
3. [Third-body initialization variables](#third-body-initialization-variables)
4. [Third-body propagation variables](#third-body-propagation-variables)

---

#### Initialization variables

The following variables depend solely on epoch elements.

| variable                          | symbol   | description                                                                                                                                                                                  |
| :-------------------------------- | :------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Elements::datetime.year()`       | `yᵤ`     | Gregorian calendar year                                                                                                                                                                      |
| `Elements::datetime.month()`      | `mᵤ`     | Gregorian calendar month in the range `[1, 12]`                                                                                                                                              |
| `Elements::datetime.day()`        | `dᵤ`     | Gregorian calendar day in the range `[1, 31]`                                                                                                                                                |
| `Elements::datetime.hour()`       | `hᵤ`     | Hours since midnight in the range `[0, 23]`                                                                                                                                                  |
| `Elements::datetime.minute()`     | `minᵤ`   | Minutes since the hour in the range `[0, 59]`                                                                                                                                                |
| `Elements::datetime.second()`     | `sᵤ`     | Seconds since the minute in the range `[0, 59]`                                                                                                                                              |
| `Elements::datetime.nanosecond()` | `nsᵤ`    | Nanoseconds since the second in the range `[0, 10⁹[`                                                                                                                                         |
| `epoch`                           | `y₂₀₀₀`  | Julian years since UTC 1 January 2000 12h00 (J2000)                                                                                                                                          |
| `d1900`                           | `d₁₉₀₀`  | Julian days since UTC 1 January 1900 12h00 (J1900)                                                                                                                                           |
| `d1970`                           | `d₁₉₇₀`  | Julian days since UTC 1 January 1970 12h00 (J1970)                                                                                                                                           |
| `c2000`                           | `c₂₀₀₀`  | Julian centuries since UTC 1 January 2000 12h00 (J2000)                                                                                                                                      |
| `geopotential.ae`                 | `aₑ`     | equatorial radius of the earth in km                                                                                                                                                         |
| `geopotential.ke`                 | `kₑ`     | square root of the earth's gravitational parameter in earth radii³ min⁻²                                                                                                                     |
| `geopotential.j2`                 | `J₂`     | un-normalised second zonal harmonic                                                                                                                                                          |
| `geopotential.j3`                 | `J₃`     | un-normalised third zonal harmonic                                                                                                                                                           |
| `geopotential.j4`                 | `J₄`     | un-normalised fourth zonal harmonic                                                                                                                                                          |
| `kozai_mean_motion`               | `n₀`     | mean number of orbits per day (Kozai convention) at epoch in rad.min⁻¹                                                                                                                       |
| `a1`                              | `a₁`     | semi-major axis at epoch (Kozai convention)                                                                                                                                                  |
| `p0`                              | `p₀`     | partial expression of `𝛿₀` and `𝛿₁`                                                                                                                                                          |
| `d1`                              | `𝛿₁`     | used in the Kozai to Brouwer conversion                                                                                                                                                      |
| `d0`                              | `𝛿₀`     | used in the Kozai to Brouwer conversion                                                                                                                                                      |
| `B*`                              | `B*`     | radiation pressure coefficient in earth radii⁻¹                                                                                                                                              |
| `orbit_0.inclination`             | `I₀`     | angle between the equator and the orbit plane at epoch in rad                                                                                                                                |
| `orbit_0.right_ascension`         | `Ω₀`     | angle between vernal equinox and the point where the orbit crosses the equatorial plane at epoch in rad                                                                                      |
| `orbit_0.eccentricity`            | `e₀`     | shape of the orbit at epoch                                                                                                                                                                  |
| `orbit_0.argument_of_perigee`     | `ω₀`     | angle between the ascending node and the orbit's point of closest approach to the earth at epoch in rad                                                                                      |
| `orbit_0.mean_anomaly`            | `M₀`     | angle of the satellite location measured from perigee at epoch in rad                                                                                                                        |
| `orbit_0.mean_motion`             | `n₀"`    | mean number of orbits per day (Brouwer convention) at epoch in rad.min⁻¹                                                                                                                     |
| `p1`                              | `p₁`     | cosine of the inclination at epoch used in multiple expressions during initialization (`θ` in [[2]](#2), renamed to avoid confusion with the sidereal time)                                  |
| `p2`                              | `p₂`     | partial expression of multiple initialization expressions                                                                                                                                    |
| `a0`                              | `a₀"`    | semi-major axis at epoch (Brouwer convention)                                                                                                                                                |
| `p3`                              | `p₃`     | perigee in earth radii                                                                                                                                                                       |
| `p4`                              | `p₄`     | height of perigee in km                                                                                                                                                                      |
| `p5`                              | `p₅`     | partial expression of `s`                                                                                                                                                                    |
| `s`                               | `s`      | altitude parameter of the atmospheric drag expression                                                                                                                                        |
| `p6`                              | `p₆`     | partial expression of the atmospheric drag                                                                                                                                                   |
| `xi`                              | `ξ`      | partial expression of multiple initialization expressions                                                                                                                                    |
| `p7`                              | `p₇`     | partial expression of multiple initialization expressions                                                                                                                                    |
| `eta`                             | `η`      | partial expression of multiple initialization expressions and of the argument of perigee and mean anomaly in eccentric high altitude near earth propagation                                  |
| `p8`                              | `p₈`     | partial expression of multiple initialization expressions                                                                                                                                    |
| `p9`                              | `p₉`     | partial expression of multiple initialization expressions                                                                                                                                    |
| `c1`                              | `C₁`     | partial expression of multiple initializationa and propagation expressions                                                                                                                   |
| `p10`                             | `p₁₀`    | partial expression of multiple initialization expressions                                                                                                                                    |
| `b0`                              | `β₀`     | partial expression of multiple initialization expressions                                                                                                                                    |
| `p11`                             | `p₁₁`    | partial expression of multiple initialization expressions                                                                                                                                    |
| `p12`                             | `p₁₂`    | partial expression of multiple initialization expressions                                                                                                                                    |
| `p13`                             | `p₁₃`    | partial expression of multiple initialization expressions                                                                                                                                    |
| `p14`                             | `p₁₄`    | partial expression of multiple initialization expressions                                                                                                                                    |
| `p15`                             | `p₁₅`    | partial expression of multiple initialization expressions                                                                                                                                    |
| `k14`                             | `k₁₄`    | first order coefficient of the argument of perigee before adding solar and lunar perturbations                                                                                               |
| `c4`                              | `C₄`     | partial expression of multiple initializationa and propagation expressions (differs from the `C₄` expression in [[2]](#2) by a factor B\*)                                                   |
| `right_ascension_dot`             | `Ω̇`      | first order coefficient of the right ascension                                                                                                                                               |
| `argument_of_perigee_dot`         | `ω̇`      | first order coefficient of the argument of perigee                                                                                                                                           |
| `mean_anomaly_dot`                | `Ṁ`      | first order coefficient of the mean anomaly                                                                                                                                                  |
| `k0`                              | `k₀`     | second order coefficient of the right ascension before adding perturbations                                                                                                                  |
| `k1`                              | `k₁`     | partial expression of the second order coefficient of the mean anomaly                                                                                                                       |
| `k2`                              | `k₂`     | partial expression of `aᵧₙ` in near earth propagation                                                                                                                                        |
| `k3`                              | `k₃`     | partial expression of `rₖ`, `ṙₖ` and `rḟₖ` in near earth propagation                                                                                                                         |
| `k4`                              | `k₄`     | partial expression of `uₖ` in near earth propagation                                                                                                                                         |
| `k5`                              | `k₅`     | partial expression of the initial Kepler variable `p₃₈` in near earth propagation                                                                                                            |
| `k6`                              | `k₆`     | partial expression of multiple initialization expressions and of `rₖ` and `rḟₖ` in near earth propagation                                                                                    |
| `d2`                              | `D₂`     | partial expression of multiple near earth initialization expressions and of the semi-major axis in near earth propagation                                                                    |
| `p16`                             | `p₁₆`    | partial expression of multiple near earth initialization expressions                                                                                                                         |
| `d3`                              | `D₃`     | partial expression of multiple near earth initialization expressions and of the semi-major axis in near earth propagation                                                                    |
| `d4`                              | `D₄`     | partial expression of multiple near earth initialization expressions and of the semi-major axis in near earth propagation                                                                    |
| `c5`                              | `C₅`     | partial expression of multiple initializationa and propagation expressions (differs from the `C₅` expression in [[2]](#2) by a factor B\*)                                                   |
| `k7`                              | `k₇`     | sine of the mean anomaly at epoch                                                                                                                                                            |
| `k8`                              | `k₈`     | partial expression of the mean anomaly third order coefficient in high altitude near earth propagation                                                                                       |
| `k9`                              | `k₉`     | partial expression of the mean anomaly fourth order coefficient in high altitude near earth propagation                                                                                      |
| `k10`                             | `k₁₀`    | partial expression of the mean anomaly fifth order coefficient in high altitude near earth propagation                                                                                       |
| `k11`                             | `k₁₁`    | partial expression of the argument of perigee and mean anomaly in eccentric high altitude near earth propagation                                                                             |
| `k12`                             | `k₁₂`    | partial expression of the argument of perigee and mean anomaly in eccentric high altitude near earth propagation                                                                             |
| `k13`                             | `k₁₃`    | partial expression of the argument of perigee and mean anomaly in eccentric high altitude near earth propagation                                                                             |
| `lunar_right_ascension_epsilon`   | `Ωₗₑ`    | lunar right ascension of the ascending node                                                                                                                                                  |
| `lunar_right_ascension_sine`      | `sin Ωₗ` | sine of the lunar right ascension of the ascending node referred to the equator                                                                                                              |
| `lunar_right_ascension_cosine`    | `cos Ωₗ` | cosine of the lunar right ascension of the ascending node referred to the equator                                                                                                            |
| `lunar_argument_of_perigee`       | `ωₗ`     | lunar argument of perigee                                                                                                                                                                    |
| `sidereal_time_0`                 | `θ₀`     | Greenwich sidereal time at epoch                                                                                                                                                             |
| `lambda_0`                        | `λ₀`     | Earth gravity resonance variable at epoch                                                                                                                                                    |
| `lambda_dot_0`                    | `λ̇₀`     | time derivative of the Earth gravity resonance variable at epoch                                                                                                                             |
| `p17`                             | `p₁₇`    | partial expression of `𝛿ᵣ₁`, `𝛿ᵣ₂` and `𝛿ᵣ₃`                                                                                                                                                 |
| `dr1`                             | `𝛿ᵣ₁`    | first Earth gravity resonance coefficient for geosynchronous orbits (`𝛿₁` in [[2]](#2), renamed to avoid confusion with `𝛿₁` used in the Kozai to Brouwer conversion)                        |
| `dr2`                             | `𝛿ᵣ₂`    | second Earth gravity resonance coefficient for geosynchronous orbits (`𝛿₂` in [[2]](#2), renamed to match `𝛿ᵣ₁`)                                                                             |
| `dr3`                             | `𝛿ᵣ₃`    | third Earth gravity resonance coefficient for geosynchronous orbits (`𝛿₃` in [[2]](#2), renamed to match `𝛿ᵣ₁`)                                                                              |
| `p18`                             | `p₁₈`    | partial expression of `D₂₂₀₋₁` and `D₂₂₁₁`                                                                                                                                                   |
| `p19`                             | `p₁₉`    | partial expression of `D₃₂₁₀` and `D₃₂₂₂`                                                                                                                                                    |
| `p20`                             | `p₂₀`    | partial expression of `D₄₄₁₀` and `D₄₄₂₂`                                                                                                                                                    |
| `p21`                             | `p₂₁`    | partial expression of `D₅₂₂₀`, `D₅₂₃₂`, `D₅₄₂₁` and `D₅₄₃₃`                                                                                                                                  |
| `f220`                            | `F₂₂₀`   | partial expression of `D₂₂₀₋₁` and `D₄₄₁₀`                                                                                                                                                   |
| `g211`                            | `G₂₁₁`   | partial expression of `D₂₂₁₁`                                                                                                                                                                |
| `g310`                            | `G₃₁₀`   | partial expression of `D₃₂₁₀`                                                                                                                                                                |
| `g322`                            | `G₃₂₂`   | partial expression of `D₃₂₂₂`                                                                                                                                                                |
| `g410`                            | `G₄₁₀`   | partial expession of `D₄₄₁₀`                                                                                                                                                                 |
| `g422`                            | `G₄₂₂`   | partial expession of `D₄₄₂₂`                                                                                                                                                                 |
| `g520`                            | `G₅₂₀`   | partial expression of `D₅₂₂₀`                                                                                                                                                                |
| `g532`                            | `G₅₃₂`   | partial expression of `D₅₂₃₂`                                                                                                                                                                |
| `g521`                            | `G₅₂₁`   | partial expression of `D₅₄₂₁`                                                                                                                                                                |
| `g533`                            | `G₅₃₃`   | partial expression of `D₅₄₃₃`                                                                                                                                                                |
| `d220₋1`                          | `D₂₂₀₋₁` | gravity resonance coefficient for Molniya orbits (the `Dₗₘₚₖ` expression in [[2]](#2) is missing a factor `l - 2p + k` from the original equation in [[4]](#4) with `k = -1` instead of `1`) |
| `d2211`                           | `D₂₂₁₁`  | gravity resonance coefficient for Molniya orbits (the `Dₗₘₚₖ` expression in [[2]](#2) is missing a factor `l - 2p + k` from the original equation in [[4]](#4))                              |
| `d3210`                           | `D₃₂₁₀`  | see `D₂₂₁₁`                                                                                                                                                                                  |
| `d3222`                           | `D₃₂₂₂`  | see `D₂₂₁₁`                                                                                                                                                                                  |
| `d4410`                           | `D₄₄₁₀`  | see `D₂₂₁₁`                                                                                                                                                                                  |
| `d4422`                           | `D₄₄₂₂`  | see `D₂₂₁₁`                                                                                                                                                                                  |
| `d5220`                           | `D₅₂₂₀`  | see `D₂₂₁₁`                                                                                                                                                                                  |
| `d5232`                           | `D₅₂₃₂`  | see `D₂₂₁₁`                                                                                                                                                                                  |
| `d5421`                           | `D₅₄₂₁`  | see `D₂₂₁₁`                                                                                                                                                                                  |
| `d5433`                           | `D₅₄₃₃`  | see `D₂₂₁₁`                                                                                                                                                                                  |

#### Propagation variables

The following expressions depend on the propagation time `t`.

| variable                      | symbol      | description                                                                                                               |
| :---------------------------- | :---------- | :------------------------------------------------------------------------------------------------------------------------ |
| `t`                           | `t`         | minutes elapsed since epoch (can be negative)                                                                             |
| `p22`                         | `p₂₂`       | right ascension of the ascending node with neither Earth gravity resonance nor Sun and Moon contributions                 |
| `p23`                         | `p₂₃`       | argument of perigee with neither high altitude drag effects, Earth gravity resonance nor Sun and Moon contributions       |
| `orbit.inclination`           | `I`         | inclination at epoch plus `t` without the short-period effects of Earth gravity                                           |
| `orbit.right_ascension`       | `Ω`         | right ascension of the ascending node at epoch plus `t` without the short-period effects of Earth gravity                 |
| `orbit.eccentricity`          | `e`         | eccentricity at epoch plus `t` without the short-period effects of Earth gravity                                          |
| `orbit.argument_of_perigee`   | `ω`         | argument of perigee at epoch plus `t` without the short-period effects of Earth gravity                                   |
| `orbit.mean_anomaly`          | `M`         | mean anomaly at epoch plus `t` without the short-period effects of Earth gravity                                          |
| `orbit.mean_motion`           | `n`         | mean motion at epoch plus `t` without the short-period effects of Earth gravity                                           |
| `a`                           | `a`         | semi-major axis                                                                                                           |
| `p32`                         | `p₃₂`       | partial expression of `aᵧₙ`                                                                                               |
| `p33`                         | `p₃₃`       | partial expression of `rₖ`, `ṙₖ` and `rḟₖ`                                                                                |
| `p34`                         | `p₃₄`       | partial expression of `uₖ`                                                                                                |
| `p35`                         | `p₃₅`       | partial expression of the initial Kepler variable `p₃₈`                                                                   |
| `p36`                         | `p₃₆`       | partial expression of `rₖ` and `rḟₖ`                                                                                      |
| `p37`                         | `p₃₇`       | partial expression of `aᵧₙ` and the initial Kepler variable `p₃₈`                                                         |
| `axn`                         | `aₓₙ`       | normalized linear eccentricity projected on the line of nodes                                                             |
| `ayn`                         | `aᵧₙ`       | normalized linear eccentricity projected on the normal to the line of nodes                                               |
| `p38`                         | `p₃₈`       | initial Kepler variable (`U` in [[2]](#2), renamed to avoid confusion with the true anomaly plus argument of perigee `u`) |
| `ew`                          | `(E + ω)ᵢ`  | Kepler variable used in an iterative process to estimate the eccentric anomaly `E`                                        |
| `delta `                      | `Δ(E + ω)ᵢ` | correction to the Kepler variable at iteration `i`                                                                        |
| `p39`                         | `p₃₉`       | eccentricity at epoch plus `t`                                                                                            |
| `pl`                          | `pₗ`        | semi-latus rectum                                                                                                         |
| `p40`                         | `p₄₀`       | normalized linear eccentricity projected on the semi-minor axis                                                           |
| `r`                           | `r`         | radius (distance to the focus) without the short-period effects of Earth gravity                                          |
| `r_dot`                       | `ṙ`         | radius time derivative without the short-period effects of Earth gravity                                                  |
| `b`                           | `β`         | semi-minor axis over semi-major axis                                                                                      |
| `p41`                         | `p₄₁`       | partial expression of `p₄₂` and `p₄₃`                                                                                     |
| `p42`                         | `p₄₂`       | sine of `u`                                                                                                               |
| `p43`                         | `p₄₃`       | cosine of `u`                                                                                                             |
| `u`                           | `u`         | true anomaly plus argument of perigee without the short-period effects of Earth gravity                                   |
| `p44`                         | `p₄₄`       | `sin(2 u)`, partial expression of `uₖ`, `Ωₖ` and `ṙₖ`                                                                     |
| `p45`                         | `p₄₅`       | `cos(2 u)`, partial expression of `rₖ`, `Iₖ` and `rḟₖ`                                                                    |
| `p46`                         | `p₄₆`       | partial expression of `rₖ`, `uₖ`, `Iₖ` and `Ωₖ`                                                                           |
| `rk`                          | `rₖ`        | radius (distance to the focus)                                                                                            |
| `uk`                          | `uₖ`        | true anomaly plus argument of perigee                                                                                     |
| `inclination_k`               | `Iₖ`        | inclination at epoch plus `t`                                                                                             |
| `right_ascension_k`           | `Ωₖ`        | right ascension at epoch plus `t`                                                                                         |
| `rk_dot`                      | `ṙₖ`        | radius time derivative                                                                                                    |
| `rfk_dot`                     | `rḟₖ`       | radius times the true anomaly derivative                                                                                  |
| `u0`                          | `u₀`        | x component of the position unit vector                                                                                   |
| `u1`                          | `u₁`        | y component of the position unit vector                                                                                   |
| `u2`                          | `u₂`        | z component of the position unit vector                                                                                   |
| `prediction.position[0]`      | `r₀`        | x component of the position vector in km (True Equator, Mean Equinox (TEME) of epoch reference frame)                     |
| `prediction.position[1]`      | `r₁`        | y component of the position vector in km (True Equator, Mean Equinox (TEME) of epoch reference frame)                     |
| `prediction.position[2]`      | `r₂`        | z component of the position vector in km (True Equator, Mean Equinox (TEME) of epoch reference frame)                     |
| `prediction.velocity[0]`      | `ṙ₀`        | x component of the velocity vector in km.s⁻¹ (True Equator, Mean Equinox (TEME) of epoch reference frame)                 |
| `prediction.velocity[1]`      | `ṙ₁`        | y component of the velocity vector in km.s⁻¹ (True Equator, Mean Equinox (TEME) of epoch reference frame)                 |
| `prediction.velocity[2]`      | `ṙ₂`        | z component of the velocity vector in km.s⁻¹ (True Equator, Mean Equinox (TEME) of epoch reference frame)                 |
| `p24`                         | `p₂₄`       | mean anomaly without drag contributions in near earth propagation                                                         |
| `p25`                         | `p₂₅`       | partial expression of `ω` and `M` in near earth propagation                                                               |
| `p26`                         | `p₂₆`       | mean anomaly with elliptic correction and without drag contributions in near earth propagation                            |
| `p27`                         | `p₂₇`       | non-clamped eccentricity in near earth propagation                                                                        |
| `p28`                         | `p₂₈`       | semi-major axis with resonance correction in deep space propagation                                                       |
| `p29`                         | `p₂₉`       | mean anomaly with resonance correction in deep space propagation                                                          |
| `p31`                         | `p₃₁`       | non-clamped eccentricity in deep space propagation                                                                        |
| `sidereal_time`               | `θ`         | sidereal time at epoch plus `t`                                                                                           |
| `delta_t`                     | `Δt`        | time step used in the integration of resonance effects of Earth gravity in min (either `720` or `-720`)                   |
| `lambda_dot`                  | `λ̇ᵢ`        | resonance effects of Earth gravity variable's time derivative at epoch plus `i Δt`                                        |
| `ni_dot`                      | `ṅᵢ`        | mean motion time derivative at epoch plus `i Δt`                                                                          |
| `ni_ddot`                     | `n̈ᵢ`        | mean motion second time derivative at epoch plus `i Δt`                                                                   |
| `ResonanceState::t`           | `tᵢ`        | resonance effects of Earth gravity integrator time (`i Δt`)                                                               |
| `ResonanceState::mean_motion` | `nᵢ`        | mean motion time derivative at epoch plus `Δt i`                                                                          |
| `ResonanceState::lambda`      | `λᵢ`        | resonance effects of Earth gravity variable at epoch plus `i Δt`                                                          |
| `p30`                         | `p₃₀`       | non-normalised `Ω` in Lyddane deep space propagation                                                                      |

#### Third-body initialization variables

The contribution of the Sun and the Moon to the orbital elements are calculated with a unique set of expressions. _src/third_body.rs_ provides a generic implementation of these expressions. Variables specific to the third body (either the Sun or the Moon) are annotated with `x`. In every other file, these variables are annotated with `s` if they correspond to solar perturbations, and `l` if they correspond to lunar perturbations.

The `aₓₙ`, `Xₓₙ`, `Zₓₙ` (`n ∈ ℕ`), `Fₓ₂` and `Fₓ₃` variables correspond to the `aₙ`, `Xₙ`, `Zₙ`, `F₂` and `F₃` variables in [[2]](#2). The added `x` highlights the dependence on the perturbing third body.

The following variables depend solely on epoch elements.

| variable                                | symbol         | description                                                                                                                                                        |
| :-------------------------------------- | :------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `third_body_inclination_sine`           | `sin Iₓ`       | sine of the inclination of the Sun (`sin Iₛ`) or the Moon (`sin Iₗ`)                                                                                               |
| `third_body_inclination_cosine`         | `cos Iₓ`       | cosine of the inclination of the Sun (`cos Iₛ`) or the Moon (`cos Iₗ`)                                                                                             |
| `delta_right_ascension_sine`            | `sin(Ω₀ - Ωₓ)` | sine of the difference between the right ascension of the ascending node of the satellite at epoch and the Sun's (`sin(Ω₀ - Ωₛ)`) or the Moon's (`sin(Ω₀ - Ωₗ)`)   |
| `delta_right_ascension_cosine`          | `cos(Ω₀ - Ωₓ)` | cosine of the difference between the right ascension of the ascending node of the satellite at epoch and the Sun's (`cos(Ω₀ - Ωₛ)`) or the Moon's (`cos(Ω₀ - Ωₗ)`) |
| `third_body_argument_of_perigee_sine`   | `sin ωₓ`       | sine of the argument of perigee of the Sun (`sin ωₛ`) or the Moon (`sin ωₗ`)                                                                                       |
| `third_body_argument_of_perigee_cosine` | `cos ωₓ`       | cosine of the argument of perigee of the Sun (`sin ωₛ`) or the Moon (`cos ωₗ`)                                                                                     |
| `third_body_mean_anomaly_0`             | `Mₓ₀`          | mean anomaly at epoch of the Sun (`Mₛ₀`) or the Moon (`Mₗ₀`)                                                                                                       |
| `ax1`                                   | `aₓ₁`          | partial expression of multiple `Xₓₙ` and `Zₓₙ` expressions                                                                                                         |
| `ax3`                                   | `aₓ₃`          | partial expression of multiple `Xₓₙ` and `Zₓₙ` expressions                                                                                                         |
| `ax7`                                   | `aₓ₇`          | partial expression of multiple `aₓ₂` and `aₓ₅`                                                                                                                     |
| `ax8`                                   | `aₓ₈`          | partial expression of multiple `aₓ₂` and `aₓ₅`                                                                                                                     |
| `ax9`                                   | `aₓ₉`          | partial expression of multiple `aₓ₄` and `aₓ₆`                                                                                                                     |
| `ax10`                                  | `aₓ₁₀`         | partial expression of multiple `aₓ₄` and `aₓ₆`                                                                                                                     |
| `ax2`                                   | `aₓ₂`          | partial expression of multiple `Xₓₙ` and `Zₓₙ` expressions                                                                                                         |
| `ax4`                                   | `aₓ₄`          | partial expression of multiple `Xₓₙ` and `Zₓₙ` expressions                                                                                                         |
| `ax5`                                   | `aₓ₅`          | partial expression of multiple `Xₓₙ` and `Zₓₙ` expressions                                                                                                         |
| `ax6`                                   | `aₓ₆`          | partial expression of multiple `Xₓₙ` and `Zₓₙ` expressions                                                                                                         |
| `xx1`                                   | `Xₓ₁`          | partial expression of multiple `Zₓₙ` expressions, `kₓ₀`, `kₓ₁` and `ėₓ`                                                                                            |
| `xx2`                                   | `Xₓ₂`          | partial expression of multiple `Zₓₙ` expressions, `kₓ₀`, `kₓ₁` and `ėₓ`                                                                                            |
| `xx3`                                   | `Xₓ₃`          | partial expression of multiple `Zₓₙ` expressions, `kₓ₀`, `kₓ₁` and `ėₓ`                                                                                            |
| `xx4`                                   | `Xₓ₄`          | partial expression of multiple `Zₓₙ` expressions, `kₓ₀`, `kₓ₁` and `ėₓ`                                                                                            |
| `xx5`                                   | `Xₓ₅`          | partial expression of multiple `Zₓₙ` expressions                                                                                                                   |
| `xx6`                                   | `Xₓ₆`          | partial expression of multiple `Zₓₙ` expressions                                                                                                                   |
| `xx7`                                   | `Xₓ₇`          | partial expression of multiple `Zₓₙ` expressions                                                                                                                   |
| `xx8`                                   | `Xₓ₈`          | partial expression of multiple `Zₓₙ` expressions                                                                                                                   |
| `zx31`                                  | `Zₓ₃₁`         | partial expression of `Zₓ₃`, `kₓ₈` and `ω̇ₓ`                                                                                                                        |
| `zx32`                                  | `Zₓ₃₂`         | partial expression of `Zₓ₂`, `kₓ₇` and `ω̇ₓ`                                                                                                                        |
| `zx33`                                  | `Zₓ₃₃`         | partial expression of `Zₓ₃`, `kₓ₈` and `ω̇ₓ`                                                                                                                        |
| `zx11`                                  | `Zₓ₁₁`         | partial expression of `kₓ₃` and `İₓ`                                                                                                                               |
| `zx13`                                  | `Zₓ₁₃`         | partial expression of `kₓ₃` and `İₓ`                                                                                                                               |
| `zx21`                                  | `Zₓ₂₁`         | partial expression of `kₓ₁₁` and `Ω̇ₓ`                                                                                                                              |
| `zx23`                                  | `Zₓ₂₃`         | partial expression of `kₓ₁₁` and `Ω̇ₓ`                                                                                                                              |
| `zx1`                                   | `Zₓ₁`          | partial expression of `kₓ₅` and `Ṁₓ`                                                                                                                               |
| `zx3`                                   | `Zₓ₃`          | partial expression of `kₓ₅` and `Ṁₓ`                                                                                                                               |
| `px0`                                   | `pₓ₀`          | partial expression of multiple `kₓₙ` expressions and `Ṁₓ`                                                                                                          |
| `px1`                                   | `pₓ₁`          | partial expression of multiple `kₓₙ` expressions and `İₓ`                                                                                                          |
| `px2`                                   | `pₓ₂`          | partial expression of multiple `kₓₙ` expressions and `ω̇ₓ`                                                                                                          |
| `px3`                                   | `pₓ₃`          | partial expression of multiple `kₓₙ` expressions and `ėₓ`                                                                                                          |
| `kx0`                                   | `kₓ₀`          | `Fₓ₂` coefficient of `δeₓ`                                                                                                                                         |
| `kx1`                                   | `kₓ₁`          | `Fₓ₃` coefficient of `δeₓ`                                                                                                                                         |
| `kx2`                                   | `kₓ₂`          | `Fₓ₂` coefficient of `δIₓ`                                                                                                                                         |
| `kx3`                                   | `kₓ₃`          | `Fₓ₃` coefficient of `δIₓ`                                                                                                                                         |
| `kx4`                                   | `kₓ₄`          | `Fₓ₂` coefficient of `δMₓ`                                                                                                                                         |
| `kx5`                                   | `kₓ₅`          | `Fₓ₃` coefficient of `δMₓ`                                                                                                                                         |
| `kx6`                                   | `kₓ₆`          | `sin fₓ` coefficient of `δMₓ`                                                                                                                                      |
| `kx7`                                   | `kₓ₇`          | `Fₓ₂` coefficient of `pₓ₄`                                                                                                                                         |
| `kx8`                                   | `kₓ₈`          | `Fₓ₃` coefficient of `pₓ₄`                                                                                                                                         |
| `kx9`                                   | `kₓ₉`          | `sin fₓ` coefficient of `pₓ₄`                                                                                                                                      |
| `kx10`                                  | `kₓ₁₀`         | `Fₓ₂` coefficient of `pₓ₅`                                                                                                                                         |
| `kx11`                                  | `kₓ₁₁`         | `Fₓ₃` coefficient of `pₓ₅`                                                                                                                                         |
| `third_body_dots.inclination`           | `İₓ`           | secular contribution of the Sun (`İₛ`) or the Moon (`İₗ`) to the inclination                                                                                       |
| `third_body_right_ascension_dot`        | `Ω̇ₓ`           | secular contribution of the Sun (`Ω̇ₛ`) or the Moon (`Ω̇ₗ`) to the right ascension of the ascending node                                                             |
| `third_body_dots.eccentricity`          | `ėₓ`           | secular contribution of the Sun (`ėₛ`) or the Moon (`ėₗ`) to the eccentricity                                                                                      |
| `third_body_dots.agument_of_perigee`    | `ω̇ₓ`           | secular contribution of the Sun (`ω̇ₛ`) or the Moon (`ω̇ₗ`) to the argument of perigee                                                                               |
| `third_body_dots.mean_anomaly`          | `Ṁₓ`           | secular contribution of the Sun (`Ṁₛ`) or the Moon (`Ṁₗ`) to the mean anomaly                                                                                      |

#### Third-body propagation variables

The following variables depend on the propagation time `t`.

| variable                        | symbol | description                                                                                                                                                                 |
| :------------------------------ | :----- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `third_body_mean_anomaly`       | `Mₓ`   | mean anomaly of the Sun (`Mₛ`) or the Moon (`Mₗ`)                                                                                                                           |
| `fx`                            | `fₓ`   | third body true anomaly                                                                                                                                                     |
| `fx2`                           | `Fₓ₂`  | partial expression of the third body long-period periodic contribution                                                                                                      |
| `fx3`                           | `Fₓ₃`  | partial expression of the third body long-period periodic contribution                                                                                                      |
| `third_body_delta_eccentricity` | `δeₓ`  | long-period periodic contribution of the Sun (`δeₛ`) or the Moon (`δeₗ`) to the eccentricity                                                                                |
| `third_body_delta_inclination`  | `δIₓ`  | long-period periodic contribution of the Sun (`δIₛ`) or the Moon (`δIₗ`) to the inclination                                                                                 |
| `third_body_delta_mean_mootion` | `δMₓ`  | long-period periodic contribution of the Sun (`δMₛ`) or the Moon (`δMₗ`) to the mean motion                                                                                 |
| `px4`                           | `pₓ₄`  | partial expression of the long-period periodic contribution of the Sun (`pₛ₄`) or the Moon (`pₗ₄`) to the right ascension of the ascending node and the argument of perigee |
| `px5`                           | `pₓ₅`  | partial expression of the long-period periodic contribution of the Sun (`pₛ₅`) or the Moon (`pₗ₅`) to the right ascension of the ascending node                             |

### Mathematical expressions

#### UT1 to Julian conversion

The epoch (Julian years since UTC 1 January 2000 12h00) can be calculated with either the AFSPC formula:

```
y₂₀₀₀ = (367 yᵤ - ⌊7 (yᵤ + ⌊(mᵤ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤ / 9⌋ + dᵤ
        + 1721013.5
        + (((nsᵤ / 10⁹ + sᵤ) / 60 + minᵤ) / 60 + hᵤ) / 24
        - 2451545)
        / 365.25
```

or the more accurate version of the same formula:

```
y₂₀₀₀ = (367 yᵤₜ₁ - ⌊7 (yᵤₜ₁ + ⌊(mᵤₜ₁ + 9) / 12⌋) / 4⌋ + 275 ⌊mᵤₜ₁ / 9⌋ + dᵤₜ₁ - 730531) / 365.25
        + (3600 hᵤₜ₁ + 60 minᵤₜ₁ + sᵤₜ₁ - 43200) / (24 × 60 × 60 × 365.25)
        + nsᵤₜ₁ / (24 × 60 × 60 × 365.25 × 10⁹)
```

#### Common initialization

```
a₁ = (kₑ / n₀)²ᐟ³

     3      3 cos²I₀
p₀ = - J₂ -----------
     4    (1 − e₀²)³ᐟ²

𝛿₁ = p₂ / a₁²

𝛿₀ = p₂ / (a₁ (1 - ¹/₃ 𝛿₁ - 𝛿₁² - ¹³⁴/₈₁ 𝛿₁³))²

n₀" = n₀ / (1 + 𝛿₀)

p₁ = cos I₀

p₂ = 1 − e₀²

k₆ = 3 p₁² - 1

a₀" = (kₑ / n₀")²ᐟ³

p₃ = a₀" (1 - e₀)

p₄ = aₑ (p₃ - 1)

p₅ = │ 20      if p₄ < 98
     │ p₄ - 78 if 98 ≤ p₄ < 156
     │ 78      otherwise

s = p₅ / aₑ + 1

p₆ = ((120 - p₅) / aₑ)⁴

ξ = 1 / (a₀" - s)

p₇ = p₆ ξ⁴

η = a₀" e₀ ξ

p₈ = |1 - η²|

p₉ = p₇ / p₈⁷ᐟ²

C₁ = B* p₉ n₀" (a₀" (1 + ³/₂ η² + e₀ η (4 + η²))
     + ³/₈ J₂ ξ k₆ (8 + 3 η² (8 + η²)) / p₈)

p₁₀ = (a₀" p₂)⁻²

β₀ = p₂¹ᐟ²

p₁₁ = ³/₂ J₂ p₁₀ n₀"

p₁₂ = ¹/₂ p₁₁ J₂ p₁₀

p₁₃ = - ¹⁵/₃₂ J₄ p₁₀² n₀"

p₁₄ = - p₁₁ p₁ + (¹/₂ p₁₂ (4 - 19 p₁²) + 2 p₁₃ (3 - 7 p₁²)) p₁

k₁₄ = - ¹/₂ p₁₁ (1 - 5 p₁²) + ¹/₁₆ p₁₂ (7 - 114 p₁² + 395 p₁⁴)

p₁₅ = n₀" + ¹/₂ p₁₁ β₀ k₆ + ¹/₁₆ p₁₂ β₀ (13 - 78 p₁² + 137 p₁⁴)

C₄ = 2 B* n₀" p₉ a₀" p₂ (
     η (2 + ¹/₂ η²)
     + e₀ (¹/₂ + 2 η²)
     - J₂ ξ / (a p₈) (-3 k₆ (1 - 2 e₀ η + η² (³/₂ - ¹/₂ e₀ η))
     + ³/₄ (1 - p₁²) (2 η² - e₀ η (1 + η²)) cos 2 ω₀)

k₀ = - ⁷/₂ p₂ p₁₁ p₁ C₁

k₁ = ³/₂ C₁

Ω̇ = │ p₁₄            if n₀" > 2π / 225
    │ p₁₄ + (Ω̇ₛ + Ω̇ₗ) otherwise

ω̇ = │ k₁₄            if n₀" > 2π / 225
    │ k₁₄ + (ω̇ₛ + ω̇ₗ) otherwise

Ṁ = │ p₁₅            if n₀" > 2π / 225
    │ p₁₅ + (Ṁₛ + Ṁₗ) otherwise
```

#### Near earth initialization

Defined only if `n₀" > 2π / 225` (near earth).

```
       1 J₃
k₂ = - - -- sin I₀
       2 J₂

k₃ = 1 - p₁²

k₄ = 7 p₁² - 1

     │   1 J₃        3 + 5 p₁
k₅ = │ - - -- sin I₀ --------    if |1 + p₁| > 1.5 × 10⁻¹²
     │   4 J₂         1 + p₁
     │   1 J₃         3 + 5 p₁
     │ - - -- sin I₀ ----------- otherwise
     │   4 J₂        1.5 × 10⁻¹²
```

#### High altitude near earth initialization

Defined only if `n₀" > 2π / 225` (near earth) and `p₃ ≥ 220 / (aₑ + 1)` (high altitude).

```
D₂ = 4 a₀" ξ C₁²

p₁₆ = D₂ ξ C₁ / 3

D₃ = (17 a + s) p₁₆

D₄ = ¹/₂ p₁₆ a₀" ξ (221 a₀" + 31 s) C₁

C₅ = 2 B* p₉ a₀" p₂ (1 + 2.75 (η² + η e₀) + e₀ η³)

k₁₁ = (1 + η cos M₀)³

k₇ = sin M₀

k₈ = D₂ + 2 C₁²

k₉ = ¹/₄ (3 D₃ + C₁ (12 D₂ + 10 C₁²))

k₁₀ = ¹/₅ (3 D₄ + 12 C₁ D₃ + 6 D₂² + 15 C₁² (2 D₂ + C₁²))
```

#### Elliptic high altitude near earth initialization

Defined only if `n₀" > 2π / 225` (near earth), `p₃ ≥ 220 / (aₑ + 1)` (high altitude) and `e₀ > 10⁻⁴` (elliptic).

```
                    J₃ p₇ ξ  n₀" sin I₀
k₁₂ = - 2 B* cos ω₀ -- ----------------
                    J₂        e₀

        2 p₇ B*
k₁₃ = - - -----
        3 e₀ η
```

#### Deep space initialization

Defined only if `n₀" ≤ 2π / 225` (deep space).

```
e₁₉₀₀ = 365.25 (t₀ + 100)

sin Iₛ = 0.39785416

cos Iₛ = 0.91744867

sin(Ω₀ - Ωₛ) = sin Ω₀

cos(Ω₀ - Ωₛ) = cos Ω₀

sin ωₛ = -0.98088458

cos ωₛ = 0.1945905

Mₛ₀ = (6.2565837 + 0.017201977 e₁₉₀₀) rem 2π

Ωₗₑ = 4.523602 - 9.2422029 × 10⁻⁴ e₁₉₀₀ rem 2π

cos Iₗ = 0.91375164 - 0.03568096 Ωₗₑ

sin Iₗ = (1 - cos²Iₗ)¹ᐟ²

sin Ωₗ = 0.089683511 sin Ωₗₑ / sin Iₗ

cos Ωₗ = (1 - sin²Ωₗ)¹ᐟ²

ωₗ = 5.8351514 + 0.001944368 e₁₉₀₀
                    0.39785416 sin Ωₗₑ / sin Iₗ
     + tan⁻¹ ------------------------------------------ - Ωₗₑ
             cos Ωₗ cos Ωₗₑ + 0.91744867 sin Ωₗ sin Ωₗₑ

sin(Ω₀ - Ωₗ) = sin Ω₀ cos Ωₗ - cos Ω₀ sin Ωₗ

cos(Ω₀ - Ωₗ) = cos Ωₗ cos Ω₀ + sin Ωₗ sin Ω₀

Mₗ₀ = (-1.1151842 + 0.228027132 e₁₉₀₀) rem 2π
```

#### Third body perturbations

Defined only if `n₀" ≤ 2π / 225` (deep space).

The following variables are evaluated for two third bodies, the Sun (solar perturbations `s`) and the Moon (lunar perturbations `l`). Variables specific to the third body are annotated with `x`. In other sections, `x` is either `s` or `l`.

```
aₓ₁ = cos ωₓ cos(Ω₀ - Ωₓ) + sin ωₓ cos Iₓ sin(Ω₀ - Ωₓ)

aₓ₃ = - sin ωₓ cos(Ω₀ - Ωₓ) + cos ωₓ cos Iₓ sin(Ω₀ - Ωₓ)

aₓ₇ = - cos ωₓ sin(Ω₀ - Ωₓ) + sin ωₓ cos Iₓ cos(Ω₀ - Ωₓ)

aₓ₈ = sin ωₓ sin Iₓ

aₓ₉ = sin ωₓ sin(Ω₀ - Ωₓ) + cos ωₓ cos Iₓ cos(Ω₀ - Ωₓ)

aₓ₁₀ = cos ωₓ sin Iₓ

aₓ₂ = aₓ₇ cos i₀ + aₓ₈ sin I₀

aₓ₄ = aₓ₉ cos i₀ + aₓ₁₀ sin I₀

aₓ₅ = - aₓ₇ sin I₀ + aₓ₈ cos I₀

aₓ₆ = - aₓ₉ sin I₀ + aₓ₁₀ cos I₀

Xₓ₁ = aₓ₁ cos ω₀ + aₓ₂ sin ω₀

Xₓ₂ = aₓ₃ cos ω₀ + aₓ₄ sin ω₀

Xₓ₃ = - aₓ₁ sin ω₀ + aₓ₂ cos ω₀

Xₓ₄ = - aₓ₃ sin ω₀ + aₓ₄ cos ω₀

Xₓ₅ = aₓ₅ sin ω₀

Xₓ₆ = aₓ₆ sin ω₀

Xₓ₇ = aₓ₅ cos ω₀

Xₓ₈ = aₓ₆ cos ω₀

Zₓ₃₁ = 12 Xₓ₁² - 3 Xₓ₃²

Zₓ₃₂ = 24 Xₓ₁ Xₓ₂ - 6 Xₓ₃ Xₓ₄

Zₓ₃₃ = 12 Xₓ₂² - 3 Xₓ₄²

Zₓ₁₁ = - 6 aₓ₁ aₓ₅ + e₀² (- 24 Xₓ₁ Xₓ₇ - 6 Xₓ₃ Xₓ₅)

Zₓ₁₃ = - 6 aₓ₃ aₓ₆ + e₀² (-24 Xₓ₂ Xₓ₈ - 6 Xₓ₄ Xₓ₆)

Zₓ₂₁ = 6 aₓ₂ aₓ₅ + e₀² (24.0 Xₓ₁ Xₓ₅ - 6 Xₓ₃ Xₓ₇)

Zₓ₂₃ = 6 aₓ₄ aₓ₆ + e₀² (24 Xₓ₂ Xₓ₆ - 6 Xₓ₄ Xₓ₈)

Zₓ₁ = 2 (3 (aₓ₁² + aₓ₂²) + Zₓ₃₁ e₀²) + p₁ Zₓ₃₁

Zₓ₃ = 2 (3 (aₓ₃² + aₓ₄²) + Zₓ₃₃ e₀²) + p₁ Zₓ₃₃

pₓ₀ = Cₓ / n₀"

        1 pₓ₀
pₓ₁ = - - ---
        2 β₀

pₓ₂ = pₓ₀ β₀

pₓ₃ = - 15 e₀ pₓ₂

Ω̇ₓ = │ 0                               if I₀ < 5.2359877 × 10⁻²
     │                                 or I₀ > π - 5.2359877 × 10⁻²
     │ - nₓ pₓ₁ (Zₓ₂₁ + Zₓ₂₃) / sin I₀ otherwise

kₓ₀ = 2 pₓ₃ (Xₓ₂ Xₓ₃ + Xₓ₁ Xₓ₄)

kₓ₁ = 2 pₓ₃ (Xₓ₂ Xₓ₄ - Xₓ₁ Xₓ₃)

kₓ₂ = 2 pₓ₁ (- 6 (aₓ₁ aₓ₆ + aₓ₃ aₓ₅) + e₀² (- 24 (Xₓ₂ Xₓ₇ + Xₓ₁ Xₓ₈) - 6 (Xₓ₃ Xₓ₆ + Xₓ₄ Xₓ₅)))

kₓ₃ = 2 pₓ₁ (Zₓ₁₃ - Zₓ₁₁)

kₓ₄ = - 2 pₓ₀ (2 (6 (aₓ₁ aₓ₃ + aₓ₂ aₓ₄) + Zₓ₃₂ e₀²) + p₁ Zₓ₃₂)

kₓ₅ = - 2 pₓ₀ (Zₓ₃ - Zₓ₁)

kₓ₆ = - 2 pₓ₀ (- 21 - 9 e₀²) eₓ

kₓ₇ = 2 pₓ₂ Zₓ₃₂

kₓ₈ = 2 pₓ₂ (Zₓ₃₃ - Zₓ₃₁)

kₓ₉ = - 18 pₓ₂ eₓ

kₓ₁₀ = - 2 pₓ₁ (6 (aₓ₄ aₓ₅ + aₓ₂ aₓ₆) + e₀² (24 (Xₓ₂ Xₓ₅ + Xₓ₁ Xₓ₆) - 6 (Xₓ₄ Xₓ₇ + Xₓ₃ Xₓ₈)))

kₓ₁₁ = - 2 pₓ₁ (Zₓ₂₃ - Zₓ₂₁)

İₓ = pₓ₁ nₓ (Zₓ₁₁ + Zₓ₁₃)

ėₓ = pₓ₃ nₓ (Xₓ₁ Xₓ₃ + Xₓ₂ Xₓ₄)

ω̇ₓ = pₓ₂ nₓ (Zₓ₃₁ + Zₓ₃₃ - 6) - cos I₀ Ω̇ₓ

Ṁₓ = - nₓ pₓ₀ (Zₓ₁ + Zₓ₃ - 14 - 6 e₀²)
```

#### Resonant deep space initialization

Defined only if `n₀" ≤ 2π / 225` (deep space) and either:

-   `0.0034906585 < n₀" < 0.0052359877` (geosynchronous)
-   `8.26 × 10⁻³ ≤ n₀" ≤ 9.24 × 10⁻³` and `e₀ ≥ 0.5` (Molniya)

The sidereal time `θ₀` at epoch can be calculated with either the IAU formula:

```
c₂₀₀₀ = y₂₀₀₀ / 100

θ₀ = ¹/₂₄₀ (π / 180) (- 6.2 × 10⁻⁶ c₂₀₀₀³ + 0.093104 c₂₀₀₀²
     + (876600 × 3600 + 8640184.812866) c₂₀₀₀ + 67310.54841) mod 2π
```

or the AFSPC formula:

```
d₁₉₇₀ = 365.25 (y₂₀₀₀ + 30) + 1

θ₀ = 1.7321343856509374 + 1.72027916940703639 × 10⁻² ⌊d₁₉₇₀ + 10⁻⁸⌋
     + (1.72027916940703639 × 10⁻² + 2π) (d₁₉₇₀ - ⌊d₁₉₇₀ + 10⁻⁸⌋)
     + 5.07551419432269442 × 10⁻¹⁵ d₁₉₇₀² mod 2π
```

```
λ₀ = │ M₀ + Ω₀ + ω₀ − θ₀ rem 2π if geosynchronous
     │ M₀ + 2 Ω₀ − 2 θ₀ rem 2π  otherwise

λ̇₀ = │ p₁₅ + (k₁₄ + p₁₄) − θ̇ + (Ṁₛ + Ṁₗ) + (ω̇ₛ + ω̇ₗ) + (Ω̇ₛ + Ω̇ₗ) - n₀" if geosynchronous
     │ p₁₅ + (Ṁₛ + Ṁₗ) + 2 (p₁₄ + (Ω̇ₛ + Ω̇ₗ) - θ̇) - n₀"                otherwise
```

#### Geosynchronous deep space initialization

Defined only if `n₀" ≤ 2π / 225` (deep space) and `0.0034906585 < n₀" < 0.0052359877` (geosynchronous orbit).

```
p₁₇ = 3 (n / a₀")²

𝛿ᵣ₁ = p₁₇ (¹⁵/₁₆ sin²I₀ (1 + 3 p₁) - ³/₄ (1 + p₁))
          (1 + 2 e₀²) 2.1460748 × 10⁻⁶ / a₀"²

𝛿ᵣ₂ = 2 p₁₇ (³/₄ (1 + p₁)²)
     (1 + e₀² (- ⁵/₂ + ¹³/₁₆ e₀²)) 1.7891679 × 10⁻⁶

𝛿ᵣ₃ = 3 p₁₇ (¹⁵/₈ (1 + p₁)³) (1 + e₀² (- 6 + 6.60937 e₀²))
      2.2123015 × 10⁻⁷ / a₀"²
```

#### Molniya deep space initialization

Defined only if `n₀" ≤ 2π / 225` (deep space) and `8.26 × 10⁻³ ≤ n₀" ≤ 9.24 × 10⁻³` and `e₀ ≥ 0.5` (Molniya).

```
p₁₈ = 3 n₀"² / a₀"²

p₁₉ = p₁₈ / a₀"

p₂₀ = p₁₉ / a₀"

p₂₁ = p₂₀ / a₀"

F₂₂₀ = ³/₄ (1 + 2 p₁ + p₁²)

G₂₁₁ = │ 3.616 - 13.247 e₀ + 16.29 e₀²                     if e₀ ≤ 0.65
       │ - 72.099 + 331.819 e₀ - 508.738 e₀² + 266.724 e₀³ otherwise

G₃₁₀ = │ - 19.302 + 117.39 e₀ - 228.419 e₀² + 156.591 e₀³      if e₀ ≤ 0.65
       │ - 346.844 + 1582.851 e₀ - 2415.925 e₀² + 1246.113 e₀³ otherwise

G₃₂₂ = │ - 18.9068 + 109.7927 e₀ - 214.6334 e₀² + 146.5816 e₀³ if e₀ ≤ 0.65
       │ - 342.585 + 1554.908 e₀ - 2366.899 e₀² + 1215.972 e₀³ otherwise

G₄₁₀ = │ - 41.122 + 242.694 e₀ - 471.094 e₀² + 313.953 e₀³      if e₀ ≤ 0.65
       │ - 1052.797 + 4758.686 e₀ - 7193.992 e₀² + 3651.957 e₀³ otherwise

G₄₂₂ = │ - 146.407 + 841.88 e₀ - 1629.014 e₀² + 1083.435 e₀³   if e₀ ≤ 0.65
       │ - 3581.69 + 16178.11 e₀ - 24462.77 e₀² + 12422.52 e₀³ otherwise

G₅₂₀ = │ - 532.114 + 3017.977 e₀ - 5740.032 e₀² + 3708.276 e₀³ if e₀ ≤ 0.65
       │ 1464.74 - 4664.75 e₀ + 3763.64 e₀²                    if 0.65 < e₀ < 0.715
       │ - 5149.66 + 29936.92 e₀ - 54087.36 e₀² + 31324.56 e₀³ otherwise

G₅₃₂ = │ - 853.666 + 4690.25 e₀ - 8624.77 e₀² + 5341.4 e₀³         if e₀ < 0.7
       │ - 40023.88 + 170470.89 e₀ - 242699.48 e₀² + 115605.82 e₀³ otherwise

G₅₂₁ = │ - 822.71072 + 4568.6173 e₀ - 8491.4146 e₀² + 5337.524 e₀³  if e₀ < 0.7
       │ - 51752.104 + 218913.95 e₀ - 309468.16 e₀² + 146349.42 e₀³ otherwise

G₅₃₃ = │ - 919.2277 + 4988.61 e₀ - 9064.77 e₀² + 5542.21 e₀³      if e₀ < 0.7
       │ - 37995.78 + 161616.52 e₀ - 229838.2 e₀² + 109377.94 e₀³ otherwise

D₂₂₀₋₁ = p₁₈ 1.7891679 × 10⁻⁶ F₂₂₀ (- 0.306 - 0.44 (e₀ - 0.64))

D₂₂₁₁ = p₁₈ 1.7891679 × 10⁻⁶ (³/₂ sin²I₀) G₂₁₁

D₃₂₁₀ = p₁₉ 3.7393792 × 10⁻⁷ (¹⁵/₈ sin I₀ (1 - 2 p₁ - 3 p₁²)) G₃₁₀

D₃₂₂₂ = p₁₉ 3.7393792 × 10⁻⁷ (- ¹⁵/₈ sin I₀ (1 + 2 p₁ - 3 p₁²)) G₃₂₂

D₄₄₁₀ = 2 p₂₀ 7.3636953 × 10⁻⁹ (35 sin²I₀ F₂₂₀) G₄₁₀

D₄₄₂₂ = 2 p₂₀ 7.3636953 × 10⁻⁹ (³¹⁵/₈ sin⁴I₀) G₄₂₂

D₅₂₂₀ = p₂₁ 1.1428639 × 10⁻⁷ (³¹⁵/₃₂ sin I₀
        (sin²I₀ (1 - 2 p₁ - 5 p₁²)
        + 0.33333333 (- 2 + 4 p₁ + 6 p₁²))) G₅₂₀

D₅₂₃₂ = p₂₁ 1.1428639 × 10⁻⁷ (sin I₀
        (4.92187512 sin²I₀ (- 2 - 4 p₁ + 10 p₁²)
        + 6.56250012 (1 + p₁ - 3 p₁²))) G₅₃₂

D₅₄₂₁ = 2 p₂₁ 2.1765803 × 10⁻⁹ (⁹⁴⁵/₃₂ sin I₀
        (2 - 8 p₁ + p₁² (- 12 + 8 p₁ + 10 p₁²))) G₅₂₁

D₅₄₃₃ = 2 p₂₁ 2.1765803 × 10⁻⁹ (⁹⁴⁵/₃₂ sin I₀
        (- 2 - 8 p₁ + p₁² (12 + 8 p₁ - 10 p₁²))) G₅₃₃
```

#### Common propagation

The following values depend on the propagation time `t` (minutes since epoch).

Named conditions have the following meaning:

-   `near earth`: `n₀" ≤ 2π / 225`
-   `low altitude near earth`: `near earth` and `p₃ < 220 / (aₑ + 1)`
-   `high altitude near earth`: `near earth` and `p₃ ≥ 220 / (aₑ + 1)`
-   `elliptic high altitude near earth`: `high altitude near earth` and `e₀ > 10⁻⁴`
-   `non-elliptic near earth`: `low altitude near earth` or `high altitude near earth` and `e₀ ≤ 10⁻⁴`
-   `deep space`: `n₀" > 2π / 225`
-   `non-Lyddane deep space`: `deep space` and `I ≥ 0.2`
-   `Lyddane deep space`: `deep space` and `I < 0.2`
-   `AFSPC Lyddane deep space`: `Lyddane deep space` and use the same expression as the original AFSPC implementation, with an `ω` discontinuity at `p₂₂ = 0`

```
p₂₂ = Ω₀ + Ω̇ t + k₀ t²

p₂₃ = ω₀ + ω̇ t

I = │ I₀                    if near earth
    │ I₀ + İ t + (δIₛ + δIₗ) otherwise

Ω = │ p₂₂                      if near earth
    │ p₂₂ + (pₛ₅ + pₗ₅) / sin I if non-Lyddane deep space
    │ p₃₀ + 2π                 if Lyddane deep space and p₃₀ + π < p₂₂ rem 2π
    │ p₃₀ - 2π                 if Lyddane deep space and p₃₀ - π > p₂₂ rem 2π
    │ p₃₀                      otherwise

e = │ 10⁻⁶              if near earth and p₂₇ < 10⁻⁶
    │ p₂₇               if near earth and p₂₇ ≥ 10⁻⁶
    │ 10⁻⁶ + (δeₛ + δeₗ) if deep space and p₃₁ < 10⁻⁶
    │ p₃₁ + (δeₛ + δeₗ)  otherwise

ω = │ p₂₃ - p₂₅                                   if elliptic high altitude near earth
    │ p₂₃                                         if non-elliptic near earth
    │ p₂₃ + (pₛ₄ + pₗ₄) - cos I (pₛ₅ + pₗ₅) / sin I if non-Lyddane deep space
    │ p₂₃ + (pₛ₄ + pₗ₄) + cos I ((p₂₂ rem 2π) - Ω)
    │ - (δIₛ + δIₗ) (p₂₂ mod 2π) sin I             if AFSPC Lyddane deep space
    │ p₂₃ + (pₛ₄ + pₗ₄) + cos I ((p₂₂ rem 2π) - Ω)
    │ - (δIₛ + δIₗ) (p₂₂ rem 2π) sin I             otherwise

M = │ p₂₆ + n₀" (k₁ t² + k₈ t³ + t⁴ (k₉ + t k₁₀) if high altitude near earth
    │ p₂₄ + n₀" k₁ t²                            if low altitude near earth
    │ p₂₉ + (δMₛ + δMₗ) + n₀" k₁ t²               otherwise


a = │ a₀" (1 - C₁ t - D₂ t² - D₃ t³ - D₄ t⁴)² if high altitude near earth
    │ a₀" (1 - C₁ t)²                         if low altitude near earth
    │ p₂₈ (1 - C₁ t)²                         otherwise

n = kₑ / a³ᐟ²

p₃₂ = │ k₂           if near earth
      │   1 J₃
      │ - - -- sin I othewise
      │   2 J₂

p₃₃ = │ k₃        if near earth
      │ 1 - cos²I othewise

p₃₄ = │ k₄          if near earth
      │ 7 cos²I - 1 otherwise

p₃₅ = │ k₅                       if near earth
      │   1 J₃       3 + 5 cos I
      │ - - -- sin I ----------- if deep space and |1 + cos I| > 1.5 × 10⁻¹²
      │   4 J₂        1 + cos I
      │   1 J₃       3 + 5 cos I
      │ - - -- sin I ----------- otherwise
      │   4 J₂       1.5 × 10⁻¹²

p₃₆ = │ k₆          if near earth
      │ 3 cos²I - 1 otherwise

p₃₇ = 1 / (a (1 - e²))

aₓₙ = e cos ω

aᵧₙ = e sin ω + p₃₇ p₃₂

p₃₈ = M + ω + p₃₇ p₃₅ aₓₙ rem 2π

(E + ω)₀ = p₃₈

            p₃₈ - aᵧₙ cos (E + ω)ᵢ + aₓₙ sin (E + ω)ᵢ - (E + ω)ᵢ
Δ(E + ω)ᵢ = ---------------------------------------------------
                  1 - cos (E + ω)ᵢ aₓₙ - sin (E + ω)ᵢ aᵧₙ

(E + ω)ᵢ₊₁ = (E + ω)ᵢ + Δ(E + ω)ᵢ|[-0.95, 0.95]

E + ω = │ (E + ω)₁₀ if ∀ j ∈ [0, 9], Δ(E + ω)ⱼ ≥ 10⁻¹²
        │ (E + ω)ⱼ  otherwise, with j the smallest integer | Δ(E + ω)ⱼ < 10⁻¹²

p₃₉ = aₓₙ² + aᵧₙ²

pₗ = a (1 - p₃₉)

p₄₀ = aₓₙ sin(E + ω) - aᵧₙ cos(E + ω)

r = a (1 - (aₓₙ cos(E + ω) + aᵧₙ sin(E + ω)))

ṙ = a¹ᐟ² p₄₀ / r

β = (1 - p₃₉)¹ᐟ²

p₄₁ = p₄₀ / (1 + β)

p₄₂ = a / r (sin(E + ω) - aᵧₙ - aₓₙ p₄₁)

p₄₃ = a / r (cos(E + ω) - aₓₙ + aᵧₙ p₄₁)

          p₄₂
u = tan⁻¹ ---
          p₄₃

p₄₄ = 2 p₄₃ p₄₂

p₄₅ = 1 - 2 p₄₂²

p₄₆ = (¹/₂ J₂ / pₗ) / pₗ

rₖ = r (1 - ³/₂ p₄₆ β p₃₆) + ¹/₂ (¹/₂ J₂ / pₗ) p₃₃ p₄₅

uₖ = u - ¹/₄ p₄₆ p₃₄ p₄₄

Ωₖ = Ω + ³/₂ p₄₆ cos I p₄₄

Iₖ = I + ³/₂ p₄₆ cos I sin I p₄₅

ṙₖ = ṙ + n (¹/₂ J₂ / pₗ) p₃₃ / kₑ

rḟₖ = pₗ¹ᐟ² / r + n (¹/₂ J₂ / pₗ) (p₃₃ p₄₅ + ³/₂ p₃₆) / kₑ

u₀ = - sin Ωₖ cos Iₖ sin uₖ + cos Ωₖ cos uₖ

u₁ = cos Ωₖ cos Iₖ sin uₖ + sin Ωₖ cos uₖ

u₂ = sin Iₖ sin uₖ

r₀ = rₖ u₀ aₑ

r₁ = rₖ u₁ aₑ

r₂ = rₖ u₂ aₑ

ṙ₀ = (ṙₖ u₀ + rḟₖ (- sin Ωₖ cos Iₖ cos uₖ - cos Ωₖ sin uₖ)) aₑ kₑ / 60

ṙ₁ = (ṙₖ u₁ + rḟₖ (cos Ωₖ cos Iₖ cos uₖ - sin Ωₖ sin uₖ)) aₑ kₑ / 60

ṙ₂ = (ṙₖ u₂ + rḟₖ (sin Iₖ cos uₖ)) aₑ kₑ / 60
```

#### Near earth propagation

Defined only if `n₀" > 2π / 225` (near earth).

```
p₂₄ = M₀ + Ṁ t

p₂₇ = | e₀ - (C₄ t + C₅ (sin p₂₆ - k₇)) if high altitude
      | e₀ - C₄ t                       otherwise
```

#### High altitude near earth propagation

Defined only if `n₀" > 2π / 225` (near earth) and `p₃ ≥ 220 / (aₑ + 1)` (high altitude).

`elliptic` means `e₀ > 10⁻⁴`.

```
p₂₅ = k₁₃ ((1 + η cos p₂₄)³ - k₁₁) + k₁₂ t

p₂₆ = │ p₂₄ + p₂₅ if elliptic
      │ p₂₄       otherwise
```

#### Deep space propagation

Defined only if `n₀" ≤ 2π / 225` (deep space).

```
p₂₈ = │ (kₑ / (nⱼ + ṅⱼ (t - tⱼ) + ¹/₂ n̈ⱼ (t - tⱼ)²))²ᐟ³ if geosynchronous or Molniya
      │ a₀"                                            otherwise

p₂₉ = │ λⱼ + λ̇ⱼ (t - tⱼ) + ¹/₂ ṅᵢ (t - tⱼ)² - p₂₂ - p₂₃ + θ if geosynchronous
      │ λⱼ + λ̇ⱼ (t - tⱼ) + ¹/₂ ṅᵢ (t - tⱼ)² - 2 p₂₂ + 2 θ   if Molniya
      │ M₀ + Ṁ t                                            otherwise

j is │ the largest positive integer | tⱼ ≤ t  if t > 0
     │ the smallest negative integer | tⱼ ≥ t if t < 0
     │ 0                                      otherwise

p₃₁ = e₀ + ė t - C₄ t
```

#### Third body propagation

Defined only if `n₀" ≤ 2π / 225` (deep space).

The following variables are evaluated for two third bodies, the Sun (solar perturbations `s`) and the Moon (lunar perturbations `l`). Variables specific to the third body are annotated with `x`. In other sections, `x` is either `s` or `l`.

```
Mₓ = Mₓ₀ + nₓ t

fₓ = Mₓ + 2 eₓ sin Mₓ

Fₓ₂ = ¹/₂ sin²fₓ - ¹/₄

Fₓ₃ = - ¹/₂ sin fₓ cos fₓ

δeₓ = kₓ₀ Fₓ₂ + kₓ₁ Fₓ₃

δIₓ = kₓ₂ Fₓ₂ + kₓ₃ Fₓ₃

δMₓ = kₓ₄ Fₓ₂ + kₓ₅ Fₓ₃ + kₓ₆ sin fₓ

pₓ₄ = kₓ₇ Fₓ₂ + kₓ₈ Fₓ₃ + kₓ₉ sin fₓ

pₓ₅ = kₓ₁₀ Fₓ₂ + kₓ₁₁ Fₓ₃
```

#### Resonant deep space propagation

Defined only if `n₀" ≤ 2π / 225` (deep space) and either:

-   `0.0034906585 < n₀" < 0.0052359877` (geosynchronous)
-   `8.26 × 10⁻³ ≤ n₀" ≤ 9.24 × 10⁻³` and `e₀ ≥ 0.5` (Molniya)

```
θ = θ₀ + 4.37526908801129966 × 10⁻³ t rem 2π

Δt = │ |Δt|  if t > 0
     │ -|Δt| if t < 0
     │ 0     otherwise

λ̇ᵢ = nᵢ + λ̇₀

ṅᵢ = │ 𝛿ᵣ₁ sin(λᵢ - λ₃₁) + 𝛿ᵣ₂ sin(2 (λᵢ - λ₂₂)) + 𝛿ᵣ₃ sin(3 (λᵢ - λ₃₃)) if geosynchronous
     │ Σ₍ₗₘₚₖ₎ Dₗₘₚₖ sin((l - 2 p) ωᵢ + m / 2 λᵢ - Gₗₘ)                    otherwise

n̈ᵢ = │ (𝛿ᵣ₁ cos(λᵢ - λ₃₁) + 𝛿ᵣ₂ cos(2 (λᵢ - λ₂₂)) + 𝛿ᵣ₃ cos(3 (λᵢ - λ₃₃))) λ̇ᵢ if geosynchronous
     │ (Σ₍ₗₘₚₖ₎ m / 2 Dₗₘₚₖ cos((l - 2 p) ωᵢ + m / 2 λᵢ - Gₗₘ)) λ̇ᵢ               otherwise

(l, m, p, k) ∈ {(2, 2, 0, -1), (2, 2, 1, 1), (3, 2, 1, 0),
    (3, 2, 2, 2), (4, 4, 1, 0), (4, 4, 2, 2), (5, 2, 2, 0),
    (5, 2, 3, 2), (5, 4, 2, 1), (5, 4, 3, 3)}

tᵢ₊₁ = tᵢ + Δt

nᵢ₊₁ = nᵢ + ṅᵢ Δt + n̈ᵢ (Δt² / 2)

λᵢ₊₁ = λᵢ + λ̇ᵢ Δt + ṅᵢ (Δt² / 2)
```

#### Lyddane deep space propagation

Defined only if `n₀" ≤ 2π / 225` (deep space) and `I < 0.2` (Lyddane).

```
            sin I sin p₂₂ + (pₛ₅ + pₗ₅) cos p₂₂ + (δIₛ + δIₗ) cos I sin p₂₂
p₃₀ = tan⁻¹ -------------------------------------------------------------
            sin I cos p₂₂ - (pₛ₅ + pₗ₅) sin p₂₂ + (δIₛ + δIₗ) cos I cos p₂₂
```

## References

<a id="1">[1]</a> David A. Vallado, Paul Crawford, R. S. Hujsak and T. S. Kelso, "Revisiting Spacetrack Report #3", presented at the AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO, 2006 August 21–24, https://doi.org/10.2514/6.2006-6753

<a id="2">[2]</a> F. R. Hoots, P. W. Schumacher Jr. and R. A. Glover, "History of Analytical Orbit Modeling in the U. S. Space Surveillance System", Journal of Guidance, Control, and Dynamics, 27(2), 174–185, 2004, https://doi.org/10.2514/1.9161

<a id="3">[3]</a> F. R. Hoots and R. L. Roehrich, "Spacetrack Report No. 3: Models for propagation of NORAD element sets", U.S. Air Force Aerospace Defense Command, Colorado Springs, CO, 1980, https://www.celestrak.com/NORAD/documentation/

<a id="4">[4]</a> R. S. Hujsak, "A Restricted Four Body Solution for Resonating Satellites Without Drag", Project SPACETRACK, Rept. 1, U.S. Air Force Aerospace Defense Command, Colorado Springs, CO, Nov. 1979, https://doi.org/10.21236/ada081263
