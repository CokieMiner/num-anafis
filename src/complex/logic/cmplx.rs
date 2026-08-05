use crate::scalar::Scalar;
use crate::traits::Number;
use core::cmp::Ordering;
use core::fmt::{Debug, Display, Formatter, Result};
use core::hash::{Hash, Hasher};
use core::ops::{Add, Div, Mul, Neg, Sub};
use crate::number::{s, r};

/// A complex number containing a real and an imaginary scalar.
#[derive(Clone, Debug)]
pub struct Complex {
    /// The real part.
    pub re: Scalar,
    /// The imaginary part.
    pub im: Scalar,
}

impl Complex {
    /// Creates a new complex number from real and imaginary parts.
    #[must_use]
    pub fn new(re: Scalar, im: Scalar) -> Self {
        Self { re, im }
    }

    /// Creates a complex zero.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            re: s(0),
            im: s(0),
        }
    }

    /// Creates a complex one.
    #[must_use]
    pub fn one() -> Self {
        Self {
            re: s(1),
            im: s(0),
        }
    }

    /// Creates the imaginary unit $i$.
    #[must_use]
    pub fn i() -> Self {
        Self {
            re: s(0),
            im: s(1),
        }
    }

    /// Returns the complex conjugate.
    #[must_use]
    pub fn conj(&self) -> Self {
        Self {
            re: self.re.clone(),
            im: -&self.im,
        }
    }

    /// Returns the squared absolute value $|z|^2$.
    #[must_use]
    pub fn abs_sq(&self) -> Scalar {
        &self.re * &self.re + &self.im * &self.im
    }
}

// Display
impl Display for Complex {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        if self.im.is_zero() {
            write!(f, "{}", self.re)
        } else if self.re.is_zero() {
            write!(f, "{}i", self.im)
        } else if self.im.is_negative() {
            // Unary negation to handle the sign properly
            write!(f, "{} - {}i", self.re, -&self.im)
        } else {
            write!(f, "{} + {}i", self.re, self.im)
        }
    }
}

impl PartialEq for Complex {
    fn eq(&self, other: &Self) -> bool {
        self.re == other.re && self.im == other.im
    }
}

impl PartialOrd for Complex {
    fn partial_cmp(&self, _other: &Self) -> Option<Ordering> {
        // Complex numbers are not totally ordered
        None
    }
}

impl Hash for Complex {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.re.hash(state);
        self.im.hash(state);
    }
}

// ============================================================================
// Arithmetic (Refs and Values)
// ============================================================================

impl Add<&Complex> for &Complex {
    type Output = Complex;
    fn add(self, rhs: &Complex) -> Complex {
        Complex::new(&self.re + &rhs.re, &self.im + &rhs.im)
    }
}

impl Sub<&Complex> for &Complex {
    type Output = Complex;
    fn sub(self, rhs: &Complex) -> Complex {
        Complex::new(&self.re - &rhs.re, &self.im - &rhs.im)
    }
}

impl Mul<&Complex> for &Complex {
    type Output = Complex;
    fn mul(self, rhs: &Complex) -> Complex {
        Complex::new(
            &self.re * &rhs.re - &self.im * &rhs.im,
            &self.re * &rhs.im + &self.im * &rhs.re,
        )
    }
}

impl Div<&Complex> for &Complex {
    type Output = Complex;
    fn div(self, rhs: &Complex) -> Complex {
        let denom = rhs.abs_sq();
        let conj = rhs.conj();
        let num = self * &conj;
        Complex::new(&num.re / &denom, &num.im / &denom)
    }
}

impl Neg for &Complex {
    type Output = Complex;
    fn neg(self) -> Complex {
        Complex::new(-&self.re, -&self.im)
    }
}

impl Add for Complex { type Output = Self; fn add(self, rhs: Self) -> Self { &self + &rhs } }
impl Sub for Complex { type Output = Self; fn sub(self, rhs: Self) -> Self { &self - &rhs } }
impl Mul for Complex { type Output = Self; fn mul(self, rhs: Self) -> Self { &self * &rhs } }
impl Div for Complex { type Output = Self; fn div(self, rhs: Self) -> Self { &self / &rhs } }
impl Neg for Complex { type Output = Self; fn neg(self) -> Self { -&self } }

impl Add<&Self> for Complex { type Output = Self; fn add(self, rhs: &Self) -> Self { &self + rhs } }
impl Sub<&Self> for Complex { type Output = Self; fn sub(self, rhs: &Self) -> Self { &self - rhs } }
impl Mul<&Self> for Complex { type Output = Self; fn mul(self, rhs: &Self) -> Self { &self * rhs } }
impl Div<&Self> for Complex { type Output = Self; fn div(self, rhs: &Self) -> Self { &self / rhs } }

// ============================================================================
// Number Trait
// ============================================================================

fn nan_complex() -> Complex {
    let scalar_nan = s(-1).sqrt(); // sqrt of -1 returns NaN on scalar backend
    Complex::new(scalar_nan.clone(), scalar_nan)
}

impl Number for Complex {
    // Unary functions
    fn sin(&self) -> Self {
        // sin(x+iy) = sin(x)cosh(y) + i cos(x)sinh(y)
        Self::new(
            self.re.sin() * self.im.cosh(),
            self.re.cos() * self.im.sinh()
        )
    }
    fn cos(&self) -> Self {
        // cos(x+iy) = cos(x)cosh(y) - i sin(x)sinh(y)
        Self::new(
            self.re.cos() * self.im.cosh(),
            -(self.re.sin() * self.im.sinh())
        )
    }
    fn tan(&self) -> Self {
        self.sin() / self.cos()
    }
    fn cot(&self) -> Self {
        self.cos() / self.sin()
    }
    fn sec(&self) -> Self {
        Self::one() / self.cos()
    }
    fn csc(&self) -> Self {
        Self::one() / self.sin()
    }

    // Inverse trig (stubs)
    fn asin(&self) -> Self { nan_complex() }
    fn acos(&self) -> Self { nan_complex() }
    fn atan(&self) -> Self { nan_complex() }
    fn acot(&self) -> Self { nan_complex() }
    fn asec(&self) -> Self { nan_complex() }
    fn acsc(&self) -> Self { nan_complex() }

    // Hyperbolic
    fn sinh(&self) -> Self {
        // sinh(x+iy) = sinh(x)cos(y) + i cosh(x)sin(y)
        Self::new(
            self.re.sinh() * self.im.cos(),
            self.re.cosh() * self.im.sin()
        )
    }
    fn cosh(&self) -> Self {
        // cosh(x+iy) = cosh(x)cos(y) + i sinh(x)sin(y)
        Self::new(
            self.re.cosh() * self.im.cos(),
            self.re.sinh() * self.im.sin()
        )
    }
    fn tanh(&self) -> Self {
        self.sinh() / self.cosh()
    }
    fn coth(&self) -> Self {
        self.cosh() / self.sinh()
    }
    fn sech(&self) -> Self {
        Self::one() / self.cosh()
    }
    fn csch(&self) -> Self {
        Self::one() / self.sinh()
    }

    // Inverse Hyperbolic (stubs)
    fn asinh(&self) -> Self { nan_complex() }
    fn acosh(&self) -> Self { nan_complex() }
    fn atanh(&self) -> Self { nan_complex() }
    fn acoth(&self) -> Self { nan_complex() }
    fn acsch(&self) -> Self { nan_complex() }
    fn asech(&self) -> Self { nan_complex() }

    // Exponential & Log
    fn exp(&self) -> Self {
        let r_exp = self.re.exp();
        Self::new(
            &r_exp * &self.im.cos(),
            &r_exp * &self.im.sin()
        )
    }
    fn expm1(&self) -> Self {
        self.exp() - Self::one()
    }
    fn exp_neg(&self) -> Self {
        (-self).exp()
    }
    fn ln(&self) -> Self {
        let r_scalar = self.abs_sq().sqrt();
        let theta = self.im.atan2(&self.re);
        Self::new(r_scalar.ln(), theta)
    }
    fn log1p(&self) -> Self {
        (self + &Self::one()).ln()
    }

    // Roots
    fn sqrt(&self) -> Self {
        if self.im.is_zero() {
            if self.re.is_negative() {
                Self::new(s(0), (-&self.re).sqrt())
            } else {
                Self::new(self.re.sqrt(), s(0))
            }
        } else {
            let r_scalar = self.abs_sq().sqrt();
            let theta = self.im.atan2(&self.re);
            let half = r(1, 2);
            let root_r = r_scalar.sqrt();
            let half_theta = &theta * &half;
            Self::new(
                &root_r * &half_theta.cos(),
                &root_r * &half_theta.sin()
            )
        }
    }
    fn cbrt(&self) -> Self {
        // Principal branch
        let r_scalar = self.abs_sq().sqrt();
        let theta = self.im.atan2(&self.re);
        let third = r(1, 3);
        let root_r = r_scalar.cbrt();
        let third_theta = &theta * &third;
        Self::new(
            &root_r * &third_theta.cos(),
            &root_r * &third_theta.sin()
        )
    }

    // Basic Math
    fn abs(&self) -> Self {
        Self::new(self.abs_sq().sqrt(), s(0))
    }
    fn signum(&self) -> Self {
        if self.is_zero() {
            Self::zero()
        } else {
            self / &Self::new(self.abs_sq().sqrt(), s(0))
        }
    }
    fn floor(&self) -> Self {
        Self::new(self.re.floor(), self.im.floor())
    }
    fn ceil(&self) -> Self {
        Self::new(self.re.ceil(), self.im.ceil())
    }
    fn round(&self) -> Self {
        Self::new(self.re.round(), self.im.round())
    }
    fn fract(&self) -> Self {
        self - &self.floor()
    }
    fn negate(&self) -> Self {
        -self
    }

    // Special Functions (Stubs)
    fn erf(&self) -> Self { nan_complex() }
    fn erfc(&self) -> Self { nan_complex() }
    fn gamma(&self) -> Self { nan_complex() }
    fn lgamma(&self) -> Self { nan_complex() }
    fn digamma(&self) -> Self { nan_complex() }
    fn trigamma(&self) -> Self { nan_complex() }
    fn tetragamma(&self) -> Self { nan_complex() }
    fn sinc(&self) -> Self {
        if self.is_zero() {
            Self::one()
        } else {
            self.sin() / self
        }
    }
    fn elliptic_k(&self) -> Self { nan_complex() }
    fn elliptic_e(&self) -> Self { nan_complex() }
    fn zeta(&self) -> Self { nan_complex() }
    fn exp_polar(&self) -> Self { nan_complex() }

    // Binary Functions
    fn atan2(&self, _x: &Self) -> Self { nan_complex() }
    fn log_base(&self, base: &Self) -> Self {
        self.ln() / base.ln()
    }
    fn pow(&self, exp: &Self) -> Self {
        (self.ln() * exp).exp()
    }
    fn besselj(&self, _order: &Self) -> Self { nan_complex() }
    fn bessely(&self, _order: &Self) -> Self { nan_complex() }
    fn besseli(&self, _order: &Self) -> Self { nan_complex() }
    fn besselk(&self, _order: &Self) -> Self { nan_complex() }
    fn polygamma(&self, _order: &Self) -> Self { nan_complex() }
    fn beta(&self, other: &Self) -> Self {
        (self.lgamma() + other.lgamma() - (self + other).lgamma()).exp()
    }
    fn zeta_deriv(&self, _order: &Self) -> Self { nan_complex() }
    fn lambertw(&self, _n: &Self) -> Self { nan_complex() }
    fn hermite(&self, _n: &Self) -> Self { nan_complex() }
    fn assoc_legendre(&self, _l: &Self, _m: &Self) -> Self { nan_complex() }
    fn spherical_harmonic(&self, _l: &Self, _m: &Self, _phi: &Self) -> Self { nan_complex() }

    // Properties
    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
    fn is_one(&self) -> bool {
        self.re.is_one() && self.im.is_zero()
    }
    fn is_neg_one(&self) -> bool {
        self.re.is_neg_one() && self.im.is_zero()
    }
    fn is_integer(&self) -> bool {
        self.im.is_zero() && self.re.is_integer()
    }
    fn is_negative(&self) -> bool {
        // Technically not a well ordered field, but maybe fallback to real part
        self.re.is_negative() && self.im.is_zero()
    }
    fn is_positive(&self) -> bool {
        self.re.is_positive() && self.im.is_zero()
    }
    fn is_finite(&self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
    fn to_float(&self) -> Self {
        Self::new(self.re.to_float(), self.im.to_float())
    }
    fn approx_eq_number(&self, other: &Self, tolerance: &Self) -> bool {
        let diff = (self - other).abs_sq();
        // compare to tol^2
        diff.total_cmp(&(tolerance.abs_sq())) == Ordering::Less
    }
    fn total_cmp(&self, _other: &Self) -> Ordering {
        // Complex numbers aren't totally ordered; just provide a fallback (e.g. by mag, then angle)
        // For trait bounds only
        Ordering::Equal
    }
    fn num_max(&self, _other: &Self) -> Self {
        nan_complex()
    }
    fn num_min(&self, _other: &Self) -> Self {
        nan_complex()
    }
}
