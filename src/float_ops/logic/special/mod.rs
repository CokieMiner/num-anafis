//! Special mathematical functions — generic over float and integer types.
//!
//! Each sub-module works via the `SpecFloat` and `SpecInt` traits so that
//! algorithms are backend-switchable with precision-specific constants
//! behind compilation walls.
//!
//! # Bibliography
//!
//! References cited throughout the `special` sub-modules follow this key:
//! - **[Abramowitz64]** Abramowitz, M. & Stegun, I.A. (1964). *Handbook of Mathematical Functions*. National Bureau of Standards.
//! - **[Borwein87]** Borwein, J. M., & Borwein, P. B. (1987). *Pi and the AGM*. Wiley.
//! - **[Borwein00]** Borwein, J. M., Bradley, D. M., & Crandall, R. E. (2000). "Computational strategies for the Riemann zeta function." *J. Comput. Appl. Math.*, 121(1-2), 247-296.
//! - **[Carlson95]** Carlson, B. C. (1995). "Numerical computation of real or complex elliptic integrals." *Numerical Algorithms*, 10(1), 13-26.
//! - **[Cephes]** Moshier, S. L. (1989). *Methods and Programs for Mathematical Functions*. Ellis Horwood Limited. (Cephes Mathematical Library)
//! - **[Clenshaw55]** Clenshaw, C. W. (1955). "A note on the summation of Chebyshev series." *Math. Tables Aids Comput.*, 9(51), 118-120.
//! - **[Corless96]** Corless, R. M., et al. (1996). "On the Lambert W function." *Adv. Comput. Math.*, 5(1), 329-359.
//! - **[DLMF]** NIST Digital Library of Mathematical Functions. `<https://dlmf.nist.gov/>`, Release 1.1.12 of 2023-12-15. F. W. J. Olver et al., eds.
//! - **[GKP94]** Graham, R. L., Knuth, D. E., & Patashnik, O. (1994). *Concrete Mathematics* (2nd ed.). Addison-Wesley.
//! - **[Lanczos64]** Lanczos, C. (1964). "A precision approximation of the gamma function." *J. SIAM Numer. Anal. Ser. B*, 1(1), 86-96.
//! - **[Miller52]** Miller, J. C. P. (1952). "A method for the determination of converging factors..." *Proc. Camb. Philos. Soc.*, 48(2), 243-254.
use core::fmt::Debug;
use core::ops::{Add, Div, Mul, Neg, Rem, Sub};

mod bessel_ik;
mod bessel_jy;
mod beta;
mod complex;
mod elliptic;
mod erf;
mod gamma;
mod helpers;
mod hermite;
mod impls;
mod lambert_w;
mod legendre;
mod polygamma;
mod zeta;
mod zeta_deriv;

pub use bessel_ik::{besseli, besselk};
pub use bessel_jy::{besselj, bessely};
pub use beta::beta;
pub use elliptic::{elliptic_e, elliptic_k};
pub use erf::{erf, erfc};
pub use gamma::{gamma, lgamma};
pub use hermite::hermite;
pub use lambert_w::{lambertw0, lambertwm1};
pub use legendre::{assoc_legendre, spherical_harmonic};
pub use polygamma::{digamma, polygamma_n, tetragamma, trigamma};
pub use zeta::zeta;
pub use zeta_deriv::zeta_deriv;

// ============================================================================
// SpecInt — integer abstraction
// ============================================================================

pub trait SpecInt:
    Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Rem<Output = Self>
    + Debug
    + 'static
{
    fn zero() -> Self;
    fn one() -> Self;
    fn abs(self) -> Self;
    fn is_zero(self) -> bool;
    fn is_negative(self) -> bool;
    fn to_usize(self) -> usize;
    fn from_usize(v: usize) -> Self;
}

// ============================================================================
// SpecFloat — float abstraction
// ============================================================================

pub trait SpecFloat:
    Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + 'static
{
    type Int: SpecInt;

    fn eps() -> Self;
    const ERF_TERMS: usize;
    const DIGAMMA_SHIFT: usize;
    const LAMBERT_MAX_ITERATIONS: usize;
    const ZETA_BORWEIN_N: usize;
    const ZETA_EM_TERMS: usize;

    fn zero() -> Self;
    fn one() -> Self;
    fn neg_one() -> Self {
        -Self::one()
    }
    fn two() -> Self {
        Self::one() + Self::one()
    }
    fn half() -> Self {
        Self::one() / Self::two()
    }
    fn pi() -> Self;
    fn e() -> Self;
    fn frac_2_pi() -> Self;
    fn infinity() -> Self;
    fn neg_infinity() -> Self {
        -Self::infinity()
    }
    fn nan() -> Self;
    fn max_value() -> Self;

    // Cody-Waite range reduction constants for π/4 and 3π/4.
    fn pio4_hi() -> Self;
    fn pio4_lo() -> Self;
    fn pio34_hi() -> Self;
    fn pio34_lo() -> Self;

    fn from_int<I: SpecInt>(v: I) -> Self {
        let a = v.abs();
        let f = Self::from_usize(a.to_usize());
        if v.is_negative() { -f } else { f }
    }
    fn from_usize(v: usize) -> Self;

    fn abs(self) -> Self;
    fn signum(self) -> Self;
    fn sqrt(self) -> Self;
    fn sin(self) -> Self;
    fn cos(self) -> Self;
    fn atan2(self, other: Self) -> Self;
    fn max(self, other: Self) -> Self;
    fn mul_add(self, a: Self, b: Self) -> Self;
    fn exp(self) -> Self;
    fn ln(self) -> Self;
    fn fract(self) -> Self;
    fn round(self) -> Self;
    fn floor(self) -> Self;
    fn powf(self, exp: Self) -> Self;
    fn pow_int<I: SpecInt>(self, n: I) -> Self;

    fn is_nan(self) -> bool;
    fn is_infinite(self) -> bool;
    fn is_sign_negative(self) -> bool;
    fn is_sign_positive(self) -> bool {
        !self.is_sign_negative() && !self.is_nan()
    }

    // --- Coefficient arrays for special functions ---
    // Each backend stores these as const arrays in native precision.

    fn lanczos_coeffs() -> &'static [Self];
    fn bernoulli_pairs() -> &'static [(Self, Self)];

    fn stieltjes_coeffs() -> &'static [Self];
    fn stirling_coeffs() -> &'static [Self];
    fn zeta_ints() -> &'static [Self];

    fn besselj0_num_coeffs() -> &'static [Self];
    fn besselj0_den_coeffs() -> &'static [Self];
    fn besselj0_pcos_coeffs() -> &'static [Self];
    fn besselj0_psin_coeffs() -> &'static [Self];

    fn besselj1_num_coeffs() -> &'static [Self];
    fn besselj1_den_coeffs() -> &'static [Self];
    fn besselj1_pcos_coeffs() -> &'static [Self];
    fn besselj1_psin_coeffs() -> &'static [Self];

    fn bessely0_num_coeffs() -> &'static [Self];
    fn bessely0_den_coeffs() -> &'static [Self];
    fn bessely0_pcos_coeffs() -> &'static [Self];
    fn bessely0_psin_coeffs() -> &'static [Self];

    fn bessely1_num_coeffs() -> &'static [Self];
    fn bessely1_den_coeffs() -> &'static [Self];
    fn bessely1_pcos_coeffs() -> &'static [Self];
    fn bessely1_psin_coeffs() -> &'static [Self];

    fn besseli0_small_coeffs() -> &'static [Self];
    fn besseli0_large_coeffs() -> &'static [Self];

    fn besseli1_small_coeffs() -> &'static [Self];
    fn besseli1_large_coeffs() -> &'static [Self];

    fn besselk0_small_coeffs() -> &'static [Self];
    fn besselk0_large_coeffs() -> &'static [Self];

    fn besselk1_small_coeffs() -> &'static [Self];
    fn besselk1_large_coeffs() -> &'static [Self];

    fn besselj_split() -> Self;
    fn bessel_miller_seed() -> Self;

    fn besselj0_root1() -> Self;
    fn besselj0_root2() -> Self;
    fn besselj1_root1() -> Self;
    fn besselj1_root2() -> Self;
}
