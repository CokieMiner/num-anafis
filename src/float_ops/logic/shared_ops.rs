//! Shared function implementations used by both `f64_ops` and `f32_ops`.
//!
//! This module is `#[macro_use]`d by each backend so the macro below is
//! expanded in the calling module's scope, picking up the local
//! `BackingFloat`, `IntType`, `math`, `rational`, and `int_math` bindings.

#[cfg(any(backend = "64", backend = "32"))]
macro_rules! impl_shared_float_ops {
    () => {
        #[inline]
        pub(in crate::float_ops) const fn clone(value: &BackingFloat) -> BackingFloat {
            *value
        }

        #[inline]
        pub(in crate::float_ops) fn set_precision(_bits: u32) -> bool {
            false
        }

        // --- Construction & conversion ---

        pub(in crate::float_ops) fn to_rational(
            value: &BackingFloat,
        ) -> Option<crate::rational_math::RationalType> {
            use num_traits::Float;
            if !value.is_finite() {
                return None;
            }
            let (mantissa, exponent, sign) = Float::integer_decode(*value);
            if mantissa == 0 {
                return Some(crate::rational_math::from_integer(crate::int_math::zero()));
            }

            #[allow(clippy::cast_possible_wrap, reason = "Float mantissa fits in i64")]
            let mut num = crate::int_math::from_i64(mantissa as i64)?;
            if sign < 0 {
                num = crate::int_math::neg(&num)?;
            }

            let result = if exponent >= 0 {
                let mut power = crate::int_math::from_i64(1)?;
                let base_2 = crate::int_math::from_i64(2)?;
                for _ in 0..exponent {
                    power = crate::int_math::mul(&power, &base_2)?;
                }
                crate::rational_math::from_integer(crate::int_math::mul(&num, &power)?)
            } else {
                let mut den = crate::int_math::from_i64(1)?;
                let base_2 = crate::int_math::from_i64(2)?;
                for _ in 0..(-exponent) {
                    den = crate::int_math::mul(&den, &base_2)?;
                }
                crate::rational_math::new(num, den)
            };
            Some(result)
        }

        pub(in crate::float_ops) fn to_string(value: &BackingFloat) -> String {
            ::alloc::format!("{value}")
        }

        #[cfg(feature = "serde")]
        pub(in crate::float_ops) fn from_str(value: &str) -> Option<BackingFloat> {
            value.parse().ok()
        }

        // --- Arithmetic ---

        #[inline]
        pub(in crate::float_ops) fn add(lhs: &BackingFloat, rhs: &BackingFloat) -> BackingFloat {
            *lhs + *rhs
        }

        #[inline]
        pub(in crate::float_ops) fn sub(lhs: &BackingFloat, rhs: &BackingFloat) -> BackingFloat {
            *lhs - *rhs
        }

        #[inline]
        pub(in crate::float_ops) fn mul(lhs: &BackingFloat, rhs: &BackingFloat) -> BackingFloat {
            *lhs * *rhs
        }

        #[inline]
        pub(in crate::float_ops) fn div(lhs: &BackingFloat, rhs: &BackingFloat) -> BackingFloat {
            *lhs / *rhs
        }

        #[inline]
        pub(in crate::float_ops) fn neg(value: &BackingFloat) -> BackingFloat {
            -*value
        }

        // --- Comparison ---

        #[inline]
        pub(in crate::float_ops) fn cmp(
            lhs: &BackingFloat,
            rhs: &BackingFloat,
        ) -> Option<Ordering> {
            lhs.partial_cmp(rhs)
        }

        // --- Properties ---

        #[inline]
        pub(in crate::float_ops) fn is_zero(value: &BackingFloat) -> bool {
            *value == 0.0
        }

        #[inline]
        pub(in crate::float_ops) fn is_one(value: &BackingFloat) -> bool {
            *value == 1.0
        }

        #[inline]
        pub(in crate::float_ops) fn is_neg_one(value: &BackingFloat) -> bool {
            *value == -1.0
        }

        #[inline]
        pub(in crate::float_ops) fn is_integer(value: &BackingFloat) -> bool {
            value.fract() == 0.0
        }

        #[inline]
        pub(in crate::float_ops) const fn is_finite(value: &BackingFloat) -> bool {
            value.is_finite()
        }

        #[inline]
        pub(in crate::float_ops) fn is_negative(value: &BackingFloat) -> bool {
            *value < 0.0
        }

        #[inline]
        pub(in crate::float_ops) fn is_positive(value: &BackingFloat) -> bool {
            *value > 0.0
        }

        #[inline]
        pub(in crate::float_ops) fn is_nan(value: &BackingFloat) -> bool {
            value.is_nan()
        }

        // --- Basic math ---

        #[inline]
        pub(in crate::float_ops) fn abs(value: &BackingFloat) -> BackingFloat {
            math::fabs(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn signum(value: &BackingFloat) -> BackingFloat {
            value.signum()
        }

        #[inline]
        pub(in crate::float_ops) fn floor(value: &BackingFloat) -> BackingFloat {
            math::floor(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn ceil(value: &BackingFloat) -> BackingFloat {
            math::ceil(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn round(value: &BackingFloat) -> BackingFloat {
            math::round(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn fract(value: &BackingFloat) -> BackingFloat {
            value.fract()
        }

        #[inline]
        pub(in crate::float_ops) fn pow(lhs: &BackingFloat, rhs: &BackingFloat) -> BackingFloat {
            math::pow(*lhs, *rhs)
        }

        #[inline]
        pub(in crate::float_ops) fn sqrt(value: &BackingFloat) -> BackingFloat {
            math::sqrt(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn cbrt(value: &BackingFloat) -> BackingFloat {
            math::cbrt(*value)
        }

        // --- Trigonometric ---

        #[inline]
        pub(in crate::float_ops) fn sin(value: &BackingFloat) -> BackingFloat {
            math::sin(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn cos(value: &BackingFloat) -> BackingFloat {
            math::cos(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn tan(value: &BackingFloat) -> BackingFloat {
            math::tan(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn asin(value: &BackingFloat) -> BackingFloat {
            math::asin(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn acos(value: &BackingFloat) -> BackingFloat {
            math::acos(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn atan(value: &BackingFloat) -> BackingFloat {
            math::atan(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn atan2(y: &BackingFloat, x: &BackingFloat) -> BackingFloat {
            math::atan2(*y, *x)
        }

        // --- Hyperbolic ---

        #[inline]
        pub(in crate::float_ops) fn sinh(value: &BackingFloat) -> BackingFloat {
            math::sinh(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn cosh(value: &BackingFloat) -> BackingFloat {
            math::cosh(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn tanh(value: &BackingFloat) -> BackingFloat {
            math::tanh(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn asinh(value: &BackingFloat) -> BackingFloat {
            math::asinh(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn acosh(value: &BackingFloat) -> BackingFloat {
            math::acosh(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn atanh(value: &BackingFloat) -> BackingFloat {
            math::atanh(*value)
        }

        // --- Exponential & logarithmic ---

        #[inline]
        pub(in crate::float_ops) fn exp(value: &BackingFloat) -> BackingFloat {
            math::exp(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn expm1(value: &BackingFloat) -> BackingFloat {
            math::expm1(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn ln(value: &BackingFloat) -> BackingFloat {
            math::ln(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn log1p(value: &BackingFloat) -> BackingFloat {
            math::log1p(*value)
        }

        // --- Bessel functions ---

        #[inline]
        pub(in crate::float_ops) fn besselj(n: &IntType, value: &BackingFloat) -> BackingFloat {
            super::special::besselj(*n, *value)
        }

        #[inline]
        pub(in crate::float_ops) fn bessely(n: &IntType, value: &BackingFloat) -> BackingFloat {
            super::special::bessely(*n, *value)
        }

        #[inline]
        pub(in crate::float_ops) fn besseli(n: &IntType, value: &BackingFloat) -> BackingFloat {
            super::special::besseli(*n, *value)
        }

        #[inline]
        pub(in crate::float_ops) fn besselk(n: &IntType, value: &BackingFloat) -> BackingFloat {
            super::special::besselk(*n, *value)
        }

        // --- Special functions ---

        #[inline]
        pub(in crate::float_ops) fn erf(value: &BackingFloat) -> BackingFloat {
            super::special::erf(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn erfc(value: &BackingFloat) -> BackingFloat {
            super::special::erfc(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn gamma(value: &BackingFloat) -> BackingFloat {
            super::special::gamma(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn lgamma(value: &BackingFloat) -> BackingFloat {
            super::special::lgamma(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn digamma(value: &BackingFloat) -> BackingFloat {
            super::special::digamma(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn trigamma(value: &BackingFloat) -> BackingFloat {
            super::special::trigamma(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn tetragamma(value: &BackingFloat) -> BackingFloat {
            super::special::tetragamma(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn polygamma(n: &IntType, value: &BackingFloat) -> BackingFloat {
            super::special::polygamma_n(*n, *value)
        }

        #[inline]
        pub(in crate::float_ops) fn zeta(value: &BackingFloat) -> BackingFloat {
            super::special::zeta(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn zeta_deriv(n: &IntType, value: &BackingFloat) -> BackingFloat {
            super::special::zeta_deriv(*n, *value)
        }

        #[inline]
        pub(in crate::float_ops) fn lambertw(
            order: &IntType,
            value: &BackingFloat,
        ) -> BackingFloat {
            if crate::int_math::is_zero(order) {
                super::special::lambertw0(*value)
            } else if crate::int_math::is_neg_one(order) {
                super::special::lambertwm1(*value)
            } else {
                super::special::SpecFloat::nan()
            }
        }

        #[inline]
        pub(in crate::float_ops) fn beta(a: &BackingFloat, b: &BackingFloat) -> BackingFloat {
            super::special::beta(*a, *b)
        }

        #[inline]
        pub(in crate::float_ops) fn elliptic_k(value: &BackingFloat) -> BackingFloat {
            super::special::elliptic_k(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn elliptic_e(value: &BackingFloat) -> BackingFloat {
            super::special::elliptic_e(*value)
        }

        #[inline]
        pub(in crate::float_ops) fn hermite(n: &IntType, value: &BackingFloat) -> BackingFloat {
            super::special::hermite(*n, *value)
        }

        #[inline]
        pub(in crate::float_ops) fn assoc_legendre(
            l: &IntType,
            m: &IntType,
            value: &BackingFloat,
        ) -> BackingFloat {
            super::special::assoc_legendre(*l, *m, *value)
        }

        #[inline]
        pub(in crate::float_ops) fn spherical_harmonic(
            l: &IntType,
            m: &IntType,
            theta: &BackingFloat,
            phi: &BackingFloat,
        ) -> BackingFloat {
            super::special::spherical_harmonic(*l, *m, *theta, *phi)
        }
    };
}
