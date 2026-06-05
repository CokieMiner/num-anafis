//! Shared function implementations used by both `i64_math` and `i32_math`.
//!
//! This module is `#[macro_use]`d by each backend so the macro below is
//! expanded in the calling module's scope, picking up the local
//! `BackingInt` binding.

/// Generates all int operations that are identical for i32 and i64.
///
/// The caller must define `BackingInt` and its own `from_i64` before
/// invoking this macro (the two types have different i64→Self conversions).
#[cfg(any(backend = "64", backend = "32"))]
macro_rules! impl_shared_int_ops {
    () => {
        #[inline]
        pub(in crate::int_math) const fn zero() -> BackingInt {
            0
        }

        #[inline]
        pub(in crate::int_math) const fn clone(value: &BackingInt) -> BackingInt {
            *value
        }

        #[inline]
        pub(in crate::int_math) const fn abs(value: &BackingInt) -> Option<BackingInt> {
            value.checked_abs()
        }

        #[inline]
        pub(in crate::int_math) const fn neg(value: &BackingInt) -> Option<BackingInt> {
            value.checked_neg()
        }

        #[inline]
        pub(in crate::int_math) const fn add(
            lhs: &BackingInt,
            rhs: &BackingInt,
        ) -> Option<BackingInt> {
            lhs.checked_add(*rhs)
        }

        #[inline]
        pub(in crate::int_math) const fn sub(
            lhs: &BackingInt,
            rhs: &BackingInt,
        ) -> Option<BackingInt> {
            lhs.checked_sub(*rhs)
        }

        #[inline]
        pub(in crate::int_math) const fn mul(
            lhs: &BackingInt,
            rhs: &BackingInt,
        ) -> Option<BackingInt> {
            lhs.checked_mul(*rhs)
        }

        #[inline]
        pub(in crate::int_math) fn cmp(lhs: &BackingInt, rhs: &BackingInt) -> Ordering {
            lhs.cmp(rhs)
        }

        #[inline]
        pub(in crate::int_math) const fn is_zero(value: &BackingInt) -> bool {
            *value == 0
        }

        #[inline]
        pub(in crate::int_math) const fn is_one(value: &BackingInt) -> bool {
            *value == 1
        }

        #[inline]
        pub(in crate::int_math) const fn is_neg_one(value: &BackingInt) -> bool {
            *value == -1
        }

        #[inline]
        pub(in crate::int_math) const fn is_negative(value: &BackingInt) -> bool {
            *value < 0
        }

        #[inline]
        pub(in crate::int_math) const fn is_positive(value: &BackingInt) -> bool {
            *value > 0
        }

        #[inline]
        pub(in crate::int_math) const fn is_even(value: &BackingInt) -> bool {
            *value % 2 == 0
        }

        #[inline]
        pub(in crate::int_math) const fn modulo(lhs: &BackingInt, rhs: &BackingInt) -> BackingInt {
            *lhs % *rhs
        }

        #[inline]
        #[allow(clippy::integer_division, reason = "Exact division is intended")]
        pub(in crate::int_math) const fn div_exact(
            lhs: &BackingInt,
            rhs: &BackingInt,
        ) -> BackingInt {
            *lhs / *rhs
        }

        #[inline]
        pub(in crate::int_math) fn gcd(lhs: &BackingInt, rhs: &BackingInt) -> BackingInt {
            use num_integer::Integer;
            lhs.gcd(rhs)
        }

        #[inline]
        pub(in crate::int_math) fn perfect_square(value: &BackingInt) -> Option<BackingInt> {
            use num_integer::Roots;
            if *value < 0 {
                return None;
            }
            let root = value.sqrt();
            (root * root == *value).then_some(root)
        }

        #[inline]
        pub(in crate::int_math) fn perfect_cube(value: &BackingInt) -> Option<BackingInt> {
            use num_integer::Roots;
            let root = value.cbrt();
            (root * root * root == *value).then_some(root)
        }

        pub(in crate::int_math) fn to_string(value: &BackingInt) -> String {
            value.to_string()
        }

        #[cfg(feature = "serde")]
        pub(in crate::int_math) fn from_str(value: &str) -> Option<BackingInt> {
            value.parse().ok()
        }
    };
}
