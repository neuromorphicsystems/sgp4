# SGP4

The SGP4 algorithm, ported to Rust from the reference Celestrak implementation [[1]](#1).

The code was entirely refactored to leverage Rust's algebraic data types and highlight the relationship between the reference mathematical equations and the implementation [[2]](#2).

The numerical predictions are almost identical to those of the Celestrak implementation. The observed differences (less than 2 × 10⁻⁷ km for the position and 10⁻⁹ km.s⁻¹ for the velocity three and a half years after epoch) are well below the accuracy of the algorithm.

We drew inspiration from the incomplete https://github.com/natronics/rust-sgp4 to format mathematical expressions.

## Install

## Usage

## Variables and expressions

The variables used in the implementation are immutable and uniquely defined, making it easier to retrieve their mathematical expressions. The only exceptions are the four state variables used by the alogrithm's integrators.

The following table lists all the variables used in the code and their associated mathematical symbol. Where possible, we used symbols from [[2]](#2). Sub-expressions without a name in [[2]](#2) follow the convention `kₙ, n ∈ ℕ` if they are shared between initialization and propagation, and `pₙ, n ∈ ℕ` if they are local to initialization or propagation.

| variable | symbol | description |
|:--------:|:------:|:-----------:|

### Mathematical expressions

1. [Common initialization](#common-initialization)
2. [Near earth initialization](#near-earth-initialization)
3. [High altitude near earth initialization](#high-altitude-near-earth-initialization)
4. [Elliptic high altitude near earth initialization](#elliptic-high-altitude-near-earth-initialization)
5. [Deep space initialization](#deep-space-initialization)
6. [Third body perturbations](#third-body-perturbations)
7. [Resonant deep space initialization](#resonant-deep-space-initialization)
8. [Geosynchronous deep space initialization](#geosynchronous-deep-space-initialization)
9. [Molniya deep space initialization](#molniya-deep-space-initialization)
10. [Common propagation](#common-propagation)
11. [Near earth propagation](#near-earth-propagation)
12. [High altitude near earth propagation](#high-altitude-near-earth-propagation)
13. [Deep space propagation](#deep-space-propagation)
14. [Third body propagation](#third-body-propagation)
15. [Resonant deep space propagation](#resonant-deep-space-propagation)
16. [Lyddane deep space propagation](#lyddane-deep-space-propagation)

---

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

perigee = aₑ (p₃ - 1)

p₄ = │ 20           if perigee < 98
     │ perigee - 78 if 98 ≤ perigee < 156
     │ 78           otherwise

s = p₄ / aₑ + 1

p₅ = ((120 - p₄) / aₑ)⁴

ξ = 1 / (a₀" - s)

p₆ = p₅ ξ⁴

η = a₀" e₀ ξ

p₇ = |1 - η²|

p₈ = p₆ / p₇⁷ᐟ²

C₁ = B* p₈ n₀" (a₀" (1 + ³/₂ η² + e₀ η (4 + η²))
     + ³/₈ J₂ ξ k₆ (8 + 3 η² (8 + η²)) / p₇)

p₉ = (a₀" p₂)⁻²

β₀ = p₂¹ᐟ²

p₁₀ = ³/₂ J₂ p₉ n₀"

p₁₁ = ¹/₂ p₁₀ J₂ p₉

p₁₂ = - ¹⁵/₃₂ J₄ p₉² n₀"

p₁₃ = - p₁₀ p₁ + (¹/₂ p₁₁ (4 - 19 p₁²) + 2 p₁₂ (3 - 7 p₁²)) p₁

k₁₄ = - ¹/₂ p₁₀ (1 - 5 p₁²) + ¹/₁₆ p₁₁ (7 - 114 p₁² + 395 p₁⁴)

p₁₄ = n₀" + ¹/₂ p₁₀ β₀ k₆ + ¹/₁₆ p₁₁ β₀ (13 - 78 p₁² + 137 p₁⁴)

C₄ = 2 B* n₀" p₈ a₀" p₂ (
     η (2 + ¹/₂ η²)
     + e₀ (¹/₂ + 2 η²)
     - J₂ ξ / (a p₇) (-3 k₆ (1 - 2 e₀ η + η² (³/₂ - ¹/₂ e₀ η))
     + ³/₄ (1 - p₁²) (2 η² - e₀ η (1 + η²)) cos 2 ω₀)

k₀ = - ⁷/₂ p₂ p₁₀ p₁ C₁

k₁ = ³/₂ C₁

Ω̇ = │ p₁₃            if n₀" > 2π / 255
    │ p₁₃ + (Ω̇ₛ + Ω̇ₗ) otherwise

ω̇ = │ k₁₄            if n₀" > 2π / 255
    │ k₁₄ + (ω̇ₛ + ω̇ₗ) otherwise

Ṁ = │ p₁₄            if n₀" > 2π / 255
    │ p₁₄ + (Ṁₛ + Ṁₗ) otherwise
```

#### Near earth initialization
Defined only if `n₀" > 2π / 255` (near earth).
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
Defined only if `n₀" > 2π / 255` (near earth) and `p₃ ≥ 220 / (aₑ + 1)` (high altitude).
```
D₂ = 4 a₀" ξ C₁²

p₁₅ = D₂ ξ C₁ / 3

D₃ = (17 a + s) p₁₅

D₄ = ¹/₂ p₁₅ a₀" ξ (221 a₀" + 31 s) C₁

C₅ = 2 B* p₈ a₀" p₂ (1 + 2.75 (η² + η e₀) + e₀ η³)

k₇ = (1 + η cos M₀)³

k₈ = sin M₀

k₉ = D₂ + 2 C₁²

k₁₀ = ¹/₄ (3 D₃ + C₁ (12 D₂ + 10 C₁²))

k₁₁ = ¹/₅ (3 D₄ + 12 C₁ D₃ + 6 D₂² + 15 C₁² (2 D₂ + C₁²))
```

#### Elliptic high altitude near earth initialization
Defined only if `n₀" > 2π / 255` (near earth), `p₃ ≥ 220 / (aₑ + 1)` (high altitude) and `e₀ > 10⁻⁴` (elliptic).
```
                    J₃ p₆ ξ  n₀" sin I₀
k₁₂ = - 2 B* cos ω₀ -- ----------------
                    J₂        e₀

        2 p₆ B*
k₁₃ = - - -----
        3 e₀ η
```

#### Deep space initialization
Defined only if `n₀" ≤ 2π / 255` (deep space).
```
t₁₉₀₀ = 365.25 (t₀ + 100)

sin Iₛ = 0.39785416

cos Iₛ = 0.91744867

sin(Ω₀ - Ωₛ) = sin Ω₀

cos(Ω₀ - Ωₛ) = cos Ω₀

sin ωₛ = -0.98088458

cos ωₛ = 0.1945905

Mₛ₀ = (6.2565837 + 0.017201977 t₁₉₀₀) rem 2π

Ωₗₑ = 4.523602 - 9.2422029 × 10⁻⁴ t₁₉₀₀ rem 2π

cos Iₗ = 0.91375164 - 0.03568096 Ωₗₑ

sin Iₗ = (1 - cos²Iₗ)¹ᐟ²

sin Ωₗ = 0.089683511 sin Ωₗₑ / sin Iₗ

cos Ωₗ = (1 - sin²Ωₗ)¹ᐟ²

ωₗ = 5.8351514 + 0.001944368 t₁₉₀₀
                    0.39785416 sin Ωₗₑ / sin Iₗ
     + tan⁻¹ ------------------------------------------ - Ωₗₑ
             cos Ωₗ cos Ωₗₑ + 0.91744867 sin Ωₗ sin Ωₗₑ

sin(Ω₀ - Ωₗ) = sin Ω₀ cos Ωₗ - cos Ω₀ sin Ωₗ

cos(Ω₀ - Ωₗ) = cos Ωₗ cos Ω₀ + sin Ωₗ sin Ω₀

Mₗ₀ = (-1.1151842 + 0.228027132 t₁₉₀₀) rem 2π
```

#### Third body perturbations
Defined only if `n₀" ≤ 2π / 255` (deep space).

The following variables are evaluated for two third bodies, the sun (solar perturbations `s`) and the moon (lunar perturbations `l`). Variables specific to the third body are annotated with `x`. In other sections, `x` is either `s` or `l`.
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

Zₓ₁₂ = - 6 (aₓ₁ aₓ₆ + aₓ₃ aₓ₅) + e₀² (- 24 (Xₓ₂ Xₓ₇ + Xₓ₁ Xₓ₈) - 6 (Xₓ₃ Xₓ₆ + Xₓ₄ Xₓ₅))

Zₓ₁₃ = - 6 aₓ₃ aₓ₆ + e₀² (-24 Xₓ₂ Xₓ₈ - 6 Xₓ₄ Xₓ₆)

Zₓ₂₁ = 6 aₓ₂ aₓ₅ + e₀² (24.0 Xₓ₁ Xₓ₅ - 6 Xₓ₃ Xₓ₇)

Zₓ₂₂ = 6 (aₓ₄ aₓ₅ + aₓ₂ aₓ₆) + e₀² (24 (Xₓ₂ Xₓ₅ + Xₓ₁ Xₓ₆) - 6 (Xₓ₄ Xₓ₇ + Xₓ₃ Xₓ₈))

Zₓ₂₃ = 6 aₓ₄ aₓ₆ + e₀² (24 Xₓ₂ Xₓ₆ - 6 Xₓ₄ Xₓ₈)

Zₓ₁ = 2 (3 (aₓ₁² + aₓ₂²) + Zₓ₃₁ e₀²) + p₁ Zₓ₃₁

Zₓ₂ = 2 (6 (aₓ₁ aₓ₃ + aₓ₂ aₓ₄) + Zₓ₃₂ e₀²) + p₁ Zₓ₃₂

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

kₓ₂ = 2 pₓ₁ Zₓ₁₂

kₓ₃ = 2 pₓ₁ (Zₓ₁₃ - Zₓ₁₁)

kₓ₄ = - 2 pₓ₀ Zₓ₂

kₓ₅ = - 2 pₓ₀ (Zₓ₃ - Zₓ₁)

kₓ₆ = - 2 pₓ₀ (- 21 - 9 e₀²) eₓ

kₓ₇ = 2 pₓ₂ Zₓ₃₂

kₓ₈ = 2 pₓ₂ (Zₓ₃₃ - Zₓ₃₁)

kₓ₉ = - 18 pₓ₂ eₓ

kₓ₁₀ = - 2 pₓ₁ Zₓ₂₂

kₓ₁₁ = - 2 pₓ₁ (Zₓ₂₃ - Zₓ₂₁)

İₓ = pₓ₁ nₓ (Zₓ₁₁ + Zₓ₁₃)

ėₓ = pₓ₃ nₓ (Xₓ₁ Xₓ₃ + Xₓ₂ Xₓ₄)

ω̇ₓ = pₓ₂ nₓ (Zₓ₃₁ + Zₓ₃₃ - 6) - cos I₀ Ω̇ₓ

Ṁₓ = - nₓ pₓ₀ (Zₓ₁ + Zₓ₃ - 14 - 6 e₀²)
```

#### Resonant deep space initialization
Defined only if `n₀" ≤ 2π / 255` (deep space) and either:
- `0.0034906585 < n₀" < 0.0052359877` (geosynchronous)
- `8.26 × 10⁻³ ≤ n₀" ≤ 9.24 × 10⁻³` and `e₀ ≥ 0.5` (Molniya)

The sidereal time `θ₀` at epoch can be calculated with either the AFSPC formula:
```
t₁₉₇₀ = 365.25 (t₀ + 30)

θ₀ = 1.7321343856509374 + 1.72027916940703639 × 10⁻² ⌊t₁₉₇₀ + 10⁻⁸⌋
     + (1.72027916940703639 × 10⁻² + 2π) (t₁₉₇₀ - ⌊t₁₉₇₀ + 10⁻⁸⌋)
     + 5.07551419432269442 × 10⁻¹⁵ t₁₉₇₀² mod 2π
```
or the IAU formula:
```
t₂₀₀₀ = t₀ / 100

θ₀ = ¹/₂₄₀ (π / 180) (- 6.2 × 10⁻⁶ t₂₀₀₀³ + 0.093104 t₂₀₀₀²
     + (876600 × 3600 + 8640184.812866) t₂₀₀₀ + 67310.54841) mod 2π
```

```
λ₀ = │ M₀ + Ω₀ + ω₀ − θ₀ rem 2π if geosynchronous
     │ M₀ + 2 Ω₀ − 2 θ₀ rem 2π  otherwise

λ̇₀ = │ p₁₄ + (k₁₄ + p₁₃) − θ̇ + (Ṁₛ + Ṁₗ) + (ω̇ₛ + ω̇ₗ) + (Ω̇ₛ + Ω̇ₗ) - n₀" if geosynchronous
     │ p₁₄ + (Ṁₛ + Ṁₗ) + 2 (p₁₃ + (Ω̇ₛ + Ω̇ₗ) - θ̇) - n₀"                otherwise
```

#### Geosynchronous deep space initialization
Defined only if `n₀" ≤ 2π / 255` (deep space) and `0.0034906585 < n₀" < 0.0052359877` (geosynchronous orbit).
```
p₁₆ = 3 (n / a₀")²

𝛿ᵣ₁ = p₁₆ (¹⁵/₁₆ sin²I₀ (1 + 3 p₁) - ³/₄ (1 + p₁))
          (1 + 2 e₀²) 2.1460748 × 10⁻⁶ / a₀"²

𝛿ᵣ₂ = 2 p₁₆ (³/₄ (1 + p₁)²)
     (1 + e₀² (- ⁵/₂ + ¹³/₁₆ e₀²)) 1.7891679 × 10⁻⁶

𝛿ᵣ₃ = 3 p₁₆ (¹⁵/₈ (1 + p₁)³) (1 + e₀² (- 6 + 6.60937 e₀²))
      2.2123015 × 10⁻⁷ / a₀"²
```

#### Molniya deep space initialization
Defined only if `n₀" ≤ 2π / 255` (deep space) and `8.26 × 10⁻³ ≤ n₀" ≤ 9.24 × 10⁻³` and `e₀ ≥ 0.5` (Molniya).
```
p₁₇ = 3 n₀"² / a₀"²

p₁₈ = p₁₇ / a₀"

p₁₉ = p₁₈ / a₀"

p₂₀ = p₁₉ / a₀"

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

D₂₂₀₋₁ = p₁₇ 1.7891679 × 10⁻⁶ F₂₂₀ (- 0.306 - 0.44 (e₀ - 0.64))

D₂₂₁₁ = p₁₇ 1.7891679 × 10⁻⁶ (³/₂ sin²I₀) G₂₁₁

D₃₂₁₀ = p₁₈ 3.7393792 × 10⁻⁷ (¹⁵/₈ sin I₀ (1 - 2 p₁ - 3 p₁²)) G₃₁₀

D₃₂₂₂ = p₁₈ 3.7393792 × 10⁻⁷ (- ¹⁵/₈ sin I₀ (1 + 2 p₁ - 3 p₁²)) G₃₂₂

D₄₄₁₀ = 2 p₁₉ 7.3636953 × 10⁻⁹ (35 sin²I₀ F₂₂₀) G₄₁₀

D₄₄₂₂ = 2 p₁₉ 7.3636953 × 10⁻⁹ (³¹⁵/₈ sin⁴I₀) G₄₂₂

D₅₂₂₀ = p₂₀ 1.1428639 × 10⁻⁷ (³¹⁵/₃₂ sin I₀
        (sin²I₀ (1 - 2 p₁ - 5 p₁²)
        + 0.33333333 (- 2 + 4 p₁ + 6 p₁²))) G₅₂₀

D₅₂₃₂ = p₂₀ 1.1428639 × 10⁻⁷ (sin I₀
        (4.92187512 sin²I₀ (- 2 - 4 p₁ + 10 p₁²)
        + 6.56250012 (1 + p₁ - 3 p₁²))) G₅₃₂

D₅₄₂₁ = 2 p₂₀ 2.1765803 × 10⁻⁹ (⁹⁴⁵/₃₂ sin I₀
        (2 - 8 p₁ + p₁² (- 12 + 8 p₁ + 10 p₁²))) G₅₂₁

D₅₄₃₃ = 2 p₂₀ 2.1765803 × 10⁻⁹ (⁹⁴⁵/₃₂ sin I₀
        (- 2 - 8 p₁ + p₁² (12 + 8 p₁ - 10 p₁²))) G₅₃₃
```

#### Common propagation
The following values depend on the propagation time `t` (minutes since epoch).
```
p₂₁ = Ω₀ + Ω̇ t + k₀ t²

p₂₂ = ω₀ + ω̇ t

I = │ I₀                    if near earth
    │ I₀ + İ t + (δIₛ + δIₗ) otherwise

Ω = │ p₂₁                      if near earth
    │ p₂₁ + (pₛ₅ + pₗ₅) / sin I if deep space and I ≥ 0.2
    │ p₂₈ + 2π                 if deep space, I < 0.2 and p₂₈ + π < p₂₁ rem 2π
    │ p₂₈ - 2π                 if deep space, I < 0.2 and p₂₈ - π > p₂₁ rem 2π
    │ p₂₈                      otherwise

e = │ 10⁻⁶              if near earth and p₂₅ < 10⁻⁶
    │ p₂₅               if near earth
    │ 10⁻⁶ + (δeₛ + δeₗ) if deep space and p₂₉ < 10⁻⁶
    │ p₂₉ + (δeₛ + δeₗ)  otherwise

ω = │ p₂₂ - p₂₄                                   if elliptic high altitude near earth
    │ p₂₂                                         if near earth
    │ p₂₂ + (pₛ₄ + pₗ₄) - cos I (pₛ₅ + pₗ₅) / sin I if deep space and I ≥ 0.2
    │ p₂₂ + (pₛ₄ + pₗ₄) + cos I ((p₂₁ rem 2π) - Ω)
    │ - (δIₛ + δIₗ) (p₂₁ mod 2π) sin I             if deep space, I < 0.2
    │                                             and AFSPC compatibility mode
    │ p₂₂ + (pₛ₄ + pₗ₄) + cos I ((p₂₁ rem 2π) - Ω)
    │ - (δIₛ + δIₗ) (p₂₁ rem 2π) sin I             otherwise

M = │ p₂₃ + p₂₄        if elliptic high altitude near earth
    │ p₂₃              if near earth
    │ p₂₇ + (δMₛ + δMₗ) otherwise

a = │ a₀" (1 - C₁ t - D₂ t² - D₃ t³ - D₄ t⁴)² if high altitude near earth
    │ a₀" (1 - C₁ t)²                         if near earth
    │ p₂₆ (1 - C₁ t)²                         otherwise

n = kₑ / a³ᐟ²

𝕃 = │ M + n₀" (k₁ t² + k₉ t³ + t⁴ (k₁₀ + t k₁₁) if high altitude near earth
    │ p₂₃ + n₀" k₁ t²                           if near earth
    │ M + n₀" k₁ t²                             otherwise

p₃₀ = │ k₂           if near earth
      │   1 J₃
      │ - - -- sin I othewise
      │   2 J₂

p₃₁ = │ k₃        if near earth
      │ 1 - cos²I othewise

p₃₂ = │ k₄          if near earth
      │ 7 cos²I - 1 otherwise

p₃₃ = │ k₅                       if near earth
      │   1 J₃       3 + 5 cos I
      │ - - -- sin I ----------- if deep space and |1 + cos I| > 1.5 × 10⁻¹²
      │   4 J₂        1 + cos I
      │   1 J₃       3 + 5 cos I
      │ - - -- sin I ----------- otherwise
      │   4 J₂       1.5 × 10⁻¹²

p₃₄ = │ k₆          if near earth
      │ 3 cos²I - 1 otherwise

p₃₅ = 1 / (a (1 - e²))

aₓₙ = e cos ω

aᵧₙ = e sin ω + p₃₅ p₃₀

p₃₆ = 𝕃 + ω + p₃₅ p₃₃ aₓₙ rem 2π

(E + ω)₀ = p₃₆

            p₃₆ - aᵧₙ cos (E + ω)ᵢ + aₓₙ sin (E + ω)ᵢ - (E + ω)ᵢ
Δ(E + ω)ᵢ = ---------------------------------------------------
                  1 - cos (E + ω)ᵢ aₓₙ - sin (E + ω)ᵢ aᵧₙ

(E + ω)ᵢ₊₁ = (E + ω)ᵢ + Δ(E + ω)ᵢ|[-0.95, 0.95]

E + ω = │ (E + ω)₁₀ if ∀ j ∈ [0, 9], Δ(E + ω)ⱼ ≥ 10⁻¹²
        │ (E + ω)ⱼ  otherwise, with j the smallest integer | Δ(E + ω)ⱼ < 10⁻¹²

p₃₇ = aₓₙ² + aᵧₙ²

pₗ = a (1 - p₃₇)

p₃₈ = aₓₙ cos(E + ω) + aᵧₙ sin(E + ω)

p₃₉ = aₓₙ sin(E + ω) - aᵧₙ cos(E + ω)

r = a (1 - p₃₈)

ṙ = a¹ᐟ² p₃₉ / r

β = (1 - p₃₇)¹ᐟ²

p₄₀ = p₃₉ / (1 + β)

p₄₁ = a / r (sin(E + ω) - aᵧₙ - aₓₙ p₄₀)

p₄₂ = a / r (cos(E + ω) - aₓₙ + aᵧₙ p₄₀)

          p₄₁
u = tan⁻¹ ---
          p₄₂

p₄₃ = 2 p₄₂ p₄₁

p₄₄ = 1 - 2 p₄₁²

p₄₅ = (¹/₂ J₂ / pₗ) / pₗ

rₖ = r (1 - ³/₂ p₄₅ β p₃₄) + ¹/₂ (¹/₂ J₂ / pₗ) p₃₁ p₄₄

uₖ = u - ¹/₄ p₄₅ p₃₂ p₄₃

Ωₖ = Ω + ³/₂ p₄₅ cos I p₄₃

Iₖ = I + ³/₂ p₄₅ cos I sin I p₄₄

ṙₖ = ṙ + n (¹/₂ J₂ / pₗ) p₃₁ / kₑ

rḟₖ = pₗ¹ᐟ² / r + n (¹/₂ J₂ / pₗ) (p₃₁ p₄₄ + ³/₂ p₃₄) / kₑ

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
Defined only if `n₀" > 2π / 255` (near earth).
```
p₂₃ = M₀ + Ṁ t

p₂₅ = | e₀ - (C₄ t + C₅ (sin M - k₈)) if high altitude
      | e₀ - C₄ t                     otherwise
```

#### High altitude near earth propagation
Defined only if `n₀" > 2π / 255` (near earth) and `p₃ ≥ 220 / (aₑ + 1)` (high altitude).
```
p₂₄ = k₁₃ ((1 + η cos p₂₃)³ - k₇) + k₁₂ t
```

#### Deep space propagation
Defined only if `n₀" ≤ 2π / 255` (deep space).
```
p₂₆ = │ (kₑ / (nⱼ + ṅⱼ (t - tⱼ) + ¹/₂ n̈ⱼ (t - tⱼ)²))²ᐟ³ if geosynchronous or Molniya
      │ a₀"                                            otherwise

p₂₇ = │ λⱼ + λ̇ⱼ (t - tⱼ) + ¹/₂ ṅᵢ (t - tⱼ)² - p₂₁ - p₂₂ + θ if geosynchronous
      │ λⱼ + λ̇ⱼ (t - tⱼ) + ¹/₂ ṅᵢ (t - tⱼ)² - 2 p₂₁ + 2 θ   if Molniya
      │ M₀ + Ṁ t                                            otherwise

j is │ the largest positive integer | tⱼ ≤ t  if t > 0
     │ the smallest negative integer | tⱼ ≥ t if t < 0
     │ 0                                      otherwise

p₂₉ = e₀ + ė t - C₄ t
```

#### Third body propagation
Defined only if `n₀" ≤ 2π / 255` (deep space).

The following variables are evaluated for two third bodies, the sun (solar perturbations `s`) and the moon (lunar perturbations `l`). Variables specific to the third body are annotated with `x`. In other sections, `x` is either `s` or `l`.
```
Mₓ = Mₓ₀ + nₓ t

fₓ = Mₓ + 2 eₓ sin Mₓ

fₓ₂ = ¹/₂ sin²fₓ - ¹/₄

fₓ₃ = - ¹/₂ sin fₓ cos fₓ

δeₓ = kₓ₀ fₓ₂ + kₓ₁ fₓ₃

δIₓ = kₓ₂ fₓ₂ + kₓ₃ fₓ₃

δMₓ = kₓ₄ fₓ₂ + kₓ₅ fₓ₃ + kₓ₆ sin fₓ

pₓ₄ = kₓ₇ fₓ₂ + kₓ₈ fₓ₃ + kₓ₉ sin fₓ

pₓ₅ = kₓ₁₀ fₓ₂ + kₓ₁₁ fₓ₃
```

#### Resonant deep space propagation
Defined only if `n₀" ≤ 2π / 255` (deep space) and either:
- `0.0034906585 < n₀" < 0.0052359877` (geosynchronous)
- `8.26 × 10⁻³ ≤ n₀" ≤ 9.24 × 10⁻³` and `e₀ ≥ 0.5` (Molniya)
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
Defined only if `n₀" ≤ 2π / 255` (deep space) and `I < 0.2` (Lyddane).
```
            sin I sin p₂₁ + (pₛ₅ + pₗ₅) cos p₂₁ + (δIₛ + δIₗ) cos I sin p₂₁
p₂₈ = tan⁻¹ -------------------------------------------------------------
            sin I cos p₂₁ - (pₛ₅ + pₗ₅) sin p₂₁ + (δIₛ + δIₗ) cos I cos p₂₁
```

## References

<a id="1">[1]</a> David A. Vallado, Paul Crawford, R. S. Hujsak, and T.S. Kelso, "Revisiting Spacetrack Report #3", presented at the AIAA/AAS Astrodynamics Specialist Conference, Keystone, CO, 2006 August 21–24, https://celestrak.com/publications/AIAA/2006-6753/

<a id="2">[2]</a> F. R. Hoots, P. W. Schumacher Jr.  & R. A. Glover, "History of Analytical Orbit Modeling in the U. S. Space Surveillance System", Journal of Guidance, Control, and Dynamics, 27(2), 174–185, https://doi.org/10.2514/1.9161/

<a id="3">[3]</a> R. S. Hujsak, "A Restricted Four Body Solution for Resonating Satellites Without Drag", Project SPACETRACK, Rept. 1, U.S. Air Force Aerospace Defense Command, Colorado Springs, CO, Nov. 1979, https://doi.org/10.21236/ada081263/
