//! Extension trait exposing all math operations directly on the active
//! compile-time `FloatType`, so callers that don't need the full `Scalar`
//! wrapper can access every function with zero overhead.

use crate::float_ops::{
    FloatType, abs, acos, acosh, add, asin, asinh, assoc_legendre, atan, atan2, atanh, besseli,
    besselj, besselk, bessely, beta, cbrt, ceil, clone, cos, cosh, digamma, div, elliptic_e,
    elliptic_k, erf, erfc, exp, expm1, floor, fract, gamma, hermite, lambertw, lgamma, ln, log1p,
    mul, neg, polygamma, pow, round, signum, sin, sinh, spherical_harmonic, sqrt, sub, tan, tanh,
    tetragamma, trigamma, zeta, zeta_deriv,
};
use crate::int_math::IntType;

/// Extension trait providing every mathematical operation directly on the
/// compile-time active float backend (`f32`, `f64`, or `rug::Float`).
///
/// Integer parameters (Bessel orders, Lambert W branch, Polygamma order, etc.)
/// use [`IntType`] — the compile-time active integer backend — so there is
/// never a precision-losing cast from a fixed-width type.
pub trait AnafisMathExt {
    // --- Arithmetic ---
    /// Returns a copy of the value.
    #[must_use]
    fn clone(&self) -> Self;
    /// Adds `rhs` to `self`.
    #[must_use]
    fn add(&self, rhs: &Self) -> Self;
    /// Subtracts `rhs` from `self`.
    #[must_use]
    fn sub(&self, rhs: &Self) -> Self;
    /// Multiplies `self` by `rhs`.
    #[must_use]
    fn mul(&self, rhs: &Self) -> Self;
    /// Divides `self` by `rhs`.
    #[must_use]
    fn div(&self, rhs: &Self) -> Self;
    /// Negates `self`.
    #[must_use]
    fn neg(&self) -> Self;

    // --- Basic math ---
    /// Absolute value.
    #[must_use]
    fn abs(&self) -> Self;
    /// Sign of `self`: −1, 0, or +1.
    #[must_use]
    fn signum(&self) -> Self;
    /// Largest integer not greater than `self`.
    #[must_use]
    fn floor(&self) -> Self;
    /// Smallest integer not less than `self`.
    #[must_use]
    fn ceil(&self) -> Self;
    /// Rounds to the nearest integer.
    #[must_use]
    fn round(&self) -> Self;
    /// Fractional part of `self`.
    #[must_use]
    fn fract(&self) -> Self;
    /// `self` raised to the power `rhs`.
    #[must_use]
    fn pow(&self, rhs: &Self) -> Self;
    /// Square root.
    #[must_use]
    fn sqrt(&self) -> Self;
    /// Cube root.
    #[must_use]
    fn cbrt(&self) -> Self;

    // --- Trigonometric ---
    /// Sine.
    #[must_use]
    fn sin(&self) -> Self;
    /// Cosine.
    #[must_use]
    fn cos(&self) -> Self;
    /// Tangent.
    #[must_use]
    fn tan(&self) -> Self;
    /// Arcsine.
    #[must_use]
    fn asin(&self) -> Self;
    /// Arccosine.
    #[must_use]
    fn acos(&self) -> Self;
    /// Arctangent.
    #[must_use]
    fn atan(&self) -> Self;
    /// Two-argument arctangent `atan2(self, x)`.
    #[must_use]
    fn atan2(&self, x: &Self) -> Self;

    // --- Hyperbolic ---
    /// Hyperbolic sine.
    #[must_use]
    fn sinh(&self) -> Self;
    /// Hyperbolic cosine.
    #[must_use]
    fn cosh(&self) -> Self;
    /// Hyperbolic tangent.
    #[must_use]
    fn tanh(&self) -> Self;
    /// Inverse hyperbolic sine.
    #[must_use]
    fn asinh(&self) -> Self;
    /// Inverse hyperbolic cosine.
    #[must_use]
    fn acosh(&self) -> Self;
    /// Inverse hyperbolic tangent.
    #[must_use]
    fn atanh(&self) -> Self;

    // --- Exponential & logarithmic ---
    /// `e^self`.
    #[must_use]
    fn exp(&self) -> Self;
    /// `e^self − 1` (accurate near zero).
    #[must_use]
    fn expm1(&self) -> Self;
    /// Natural logarithm.
    #[must_use]
    fn ln(&self) -> Self;
    /// `ln(1 + self)` (accurate near zero).
    #[must_use]
    fn log1p(&self) -> Self;

    // --- Special functions (scalar) ---
    /// Error function `erf(self)`.
    #[must_use]
    fn erf(&self) -> Self;
    /// Complementary error function `erfc(self) = 1 − erf(self)`.
    #[must_use]
    fn erfc(&self) -> Self;
    /// Euler Gamma function `Γ(self)`.
    #[must_use]
    fn gamma(&self) -> Self;
    /// Natural logarithm of the Gamma function `ln|Γ(self)|`.
    #[must_use]
    fn lgamma(&self) -> Self;
    /// Digamma function `ψ(self) = Γ′(self)/Γ(self)`.
    #[must_use]
    fn digamma(&self) -> Self;
    /// Trigamma function `ψ¹(self)`.
    #[must_use]
    fn trigamma(&self) -> Self;
    /// Tetragamma function `ψ²(self)`.
    #[must_use]
    fn tetragamma(&self) -> Self;
    /// Complete elliptic integral of the first kind `K(self)`.
    #[must_use]
    fn elliptic_k(&self) -> Self;
    /// Complete elliptic integral of the second kind `E(self)`.
    #[must_use]
    fn elliptic_e(&self) -> Self;
    /// Riemann Zeta function `ζ(self)`.
    #[must_use]
    fn zeta(&self) -> Self;
    /// Beta function `B(self, b)`.
    #[must_use]
    fn beta(&self, b: &Self) -> Self;

    // --- Special functions with integer parameters ---
    /// Lambert W function. `order`: 0 = principal W₀, −1 = W₋₁ branch.
    #[must_use]
    fn lambertw(&self, order: IntType) -> Self;
    /// Bessel function of the first kind `J_n(self)`.
    #[must_use]
    fn besselj(&self, n: IntType) -> Self;
    /// Bessel function of the second kind `Y_n(self)`.
    #[must_use]
    fn bessely(&self, n: IntType) -> Self;
    /// Modified Bessel function of the first kind `I_n(self)`.
    #[must_use]
    fn besseli(&self, n: IntType) -> Self;
    /// Modified Bessel function of the second kind `K_n(self)`.
    #[must_use]
    fn besselk(&self, n: IntType) -> Self;
    /// Polygamma function `ψ^(n)(self)`.
    #[must_use]
    fn polygamma(&self, n: IntType) -> Self;
    /// n-th derivative of the Riemann Zeta function at `self`.
    #[must_use]
    fn zeta_deriv(&self, n: IntType) -> Self;
    /// Physicist's Hermite polynomial `H_n(self)`.
    #[must_use]
    fn hermite(&self, n: IntType) -> Self;
    /// Associated Legendre polynomial `P_l^m(self)`.
    #[must_use]
    fn assoc_legendre(&self, l: IntType, m: IntType) -> Self;
    /// Real spherical harmonic `Y_l^m(theta=self, phi)`.
    #[must_use]
    fn spherical_harmonic(&self, l: IntType, m: IntType, phi: &Self) -> Self;
}

impl AnafisMathExt for FloatType {
    // --- Arithmetic ---
    #[inline]
    fn clone(&self) -> Self {
        clone(self)
    }
    #[inline]
    fn add(&self, rhs: &Self) -> Self {
        add(self, rhs)
    }
    #[inline]
    fn sub(&self, rhs: &Self) -> Self {
        sub(self, rhs)
    }
    #[inline]
    fn mul(&self, rhs: &Self) -> Self {
        mul(self, rhs)
    }
    #[inline]
    fn div(&self, rhs: &Self) -> Self {
        div(self, rhs)
    }
    #[inline]
    fn neg(&self) -> Self {
        neg(self)
    }

    // --- Basic math ---
    #[inline]
    fn abs(&self) -> Self {
        abs(self)
    }
    #[inline]
    fn signum(&self) -> Self {
        signum(self)
    }
    #[inline]
    fn floor(&self) -> Self {
        floor(self)
    }
    #[inline]
    fn ceil(&self) -> Self {
        ceil(self)
    }
    #[inline]
    fn round(&self) -> Self {
        round(self)
    }
    #[inline]
    fn fract(&self) -> Self {
        fract(self)
    }
    #[inline]
    fn pow(&self, rhs: &Self) -> Self {
        pow(self, rhs)
    }
    #[inline]
    fn sqrt(&self) -> Self {
        sqrt(self)
    }
    #[inline]
    fn cbrt(&self) -> Self {
        cbrt(self)
    }

    // --- Trigonometric ---
    #[inline]
    fn sin(&self) -> Self {
        sin(self)
    }
    #[inline]
    fn cos(&self) -> Self {
        cos(self)
    }
    #[inline]
    fn tan(&self) -> Self {
        tan(self)
    }
    #[inline]
    fn asin(&self) -> Self {
        asin(self)
    }
    #[inline]
    fn acos(&self) -> Self {
        acos(self)
    }
    #[inline]
    fn atan(&self) -> Self {
        atan(self)
    }
    #[inline]
    fn atan2(&self, x: &Self) -> Self {
        atan2(self, x)
    }

    // --- Hyperbolic ---
    #[inline]
    fn sinh(&self) -> Self {
        sinh(self)
    }
    #[inline]
    fn cosh(&self) -> Self {
        cosh(self)
    }
    #[inline]
    fn tanh(&self) -> Self {
        tanh(self)
    }
    #[inline]
    fn asinh(&self) -> Self {
        asinh(self)
    }
    #[inline]
    fn acosh(&self) -> Self {
        acosh(self)
    }
    #[inline]
    fn atanh(&self) -> Self {
        atanh(self)
    }

    // --- Exponential & logarithmic ---
    #[inline]
    fn exp(&self) -> Self {
        exp(self)
    }
    #[inline]
    fn expm1(&self) -> Self {
        expm1(self)
    }
    #[inline]
    fn ln(&self) -> Self {
        ln(self)
    }
    #[inline]
    fn log1p(&self) -> Self {
        log1p(self)
    }

    // --- Special functions (scalar) ---
    #[inline]
    fn erf(&self) -> Self {
        erf(self)
    }
    #[inline]
    fn erfc(&self) -> Self {
        erfc(self)
    }
    #[inline]
    fn gamma(&self) -> Self {
        gamma(self)
    }
    #[inline]
    fn lgamma(&self) -> Self {
        lgamma(self)
    }
    #[inline]
    fn digamma(&self) -> Self {
        digamma(self)
    }
    #[inline]
    fn trigamma(&self) -> Self {
        trigamma(self)
    }
    #[inline]
    fn tetragamma(&self) -> Self {
        tetragamma(self)
    }
    #[inline]
    fn elliptic_k(&self) -> Self {
        elliptic_k(self)
    }
    #[inline]
    fn elliptic_e(&self) -> Self {
        elliptic_e(self)
    }
    #[inline]
    fn zeta(&self) -> Self {
        zeta(self)
    }
    #[inline]
    fn beta(&self, b: &Self) -> Self {
        beta(self, b)
    }

    // --- Special functions with integer parameters ---
    #[inline]
    fn lambertw(&self, order: IntType) -> Self {
        lambertw(&order, self)
    }
    #[inline]
    fn besselj(&self, n: IntType) -> Self {
        besselj(&n, self)
    }
    #[inline]
    fn bessely(&self, n: IntType) -> Self {
        bessely(&n, self)
    }
    #[inline]
    fn besseli(&self, n: IntType) -> Self {
        besseli(&n, self)
    }
    #[inline]
    fn besselk(&self, n: IntType) -> Self {
        besselk(&n, self)
    }
    #[inline]
    fn polygamma(&self, n: IntType) -> Self {
        polygamma(&n, self)
    }
    #[inline]
    fn zeta_deriv(&self, n: IntType) -> Self {
        zeta_deriv(&n, self)
    }
    #[inline]
    fn hermite(&self, n: IntType) -> Self {
        hermite(&n, self)
    }
    #[inline]
    fn assoc_legendre(&self, l: IntType, m: IntType) -> Self {
        assoc_legendre(&l, &m, self)
    }
    #[inline]
    fn spherical_harmonic(&self, l: IntType, m: IntType, phi: &Self) -> Self {
        spherical_harmonic(&l, &m, self, phi)
    }
}
