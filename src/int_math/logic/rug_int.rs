#![allow(
    clippy::unnecessary_wraps,
    reason = "Option<IntType> required for uniform API with overflow-checked backends"
)]
use alloc::string::{String, ToString};
use core::cmp::Ordering;
use rug::Integer;

pub(in crate::int_math) type BackingInt = Integer;

pub(in crate::int_math) const fn zero() -> BackingInt {
    Integer::new()
}

pub(in crate::int_math) fn clone(value: &BackingInt) -> BackingInt {
    value.clone()
}

pub(in crate::int_math) fn from_i64(value: i64) -> Option<BackingInt> {
    Some(Integer::from(value))
}

pub(in crate::int_math) fn abs(value: &BackingInt) -> Option<BackingInt> {
    Some(value.clone().abs())
}

pub(in crate::int_math) fn neg(value: &BackingInt) -> Option<BackingInt> {
    Some(Integer::from(-value))
}

pub(in crate::int_math) fn add(lhs: &BackingInt, rhs: &BackingInt) -> Option<BackingInt> {
    Some(Integer::from(lhs + rhs))
}

pub(in crate::int_math) fn sub(lhs: &BackingInt, rhs: &BackingInt) -> Option<BackingInt> {
    Some(Integer::from(lhs - rhs))
}

pub(in crate::int_math) fn mul(lhs: &BackingInt, rhs: &BackingInt) -> Option<BackingInt> {
    Some(Integer::from(lhs * rhs))
}

pub(in crate::int_math) fn cmp(lhs: &BackingInt, rhs: &BackingInt) -> Ordering {
    lhs.cmp(rhs)
}

pub(in crate::int_math) fn is_zero(value: &BackingInt) -> bool {
    *value == 0
}

pub(in crate::int_math) fn is_one(value: &BackingInt) -> bool {
    *value == 1
}

pub(in crate::int_math) fn is_neg_one(value: &BackingInt) -> bool {
    *value == -1
}

pub(in crate::int_math) fn is_negative(value: &BackingInt) -> bool {
    value.cmp0() == Ordering::Less
}

pub(in crate::int_math) fn is_positive(value: &BackingInt) -> bool {
    value.cmp0() == Ordering::Greater
}

pub(in crate::int_math) const fn is_even(value: &BackingInt) -> bool {
    value.is_even()
}

pub(in crate::int_math) fn modulo(lhs: &BackingInt, rhs: &BackingInt) -> BackingInt {
    Integer::from(lhs % rhs)
}

#[allow(clippy::integer_division, reason = "Exact division is intended")]
pub(in crate::int_math) fn div_exact(lhs: &BackingInt, rhs: &BackingInt) -> BackingInt {
    Integer::from(lhs / rhs)
}

/// GCD using GMP's highly optimized algorithm.
pub(in crate::int_math) fn gcd(lhs: &BackingInt, rhs: &BackingInt) -> BackingInt {
    let result = lhs.clone().gcd(rhs);
    if result == 0 {
        Integer::from(1)
    } else {
        result
    }
}

/// Perfect square check using GMP's `is_perfect_square` + `sqrt`.
pub(in crate::int_math) fn perfect_square(value: &BackingInt) -> Option<BackingInt> {
    if is_negative(value) {
        return None;
    }
    value.is_perfect_square().then(|| value.clone().sqrt())
}

/// Perfect cube check using GMP's `root` + verification.
pub(in crate::int_math) fn perfect_cube(value: &BackingInt) -> Option<BackingInt> {
    let abs_val = value.clone().abs();
    let (root, exact) = abs_val.root_rem(Integer::new(), 3);

    if exact == 0 {
        if is_negative(value) {
            Some(-root)
        } else {
            Some(root)
        }
    } else {
        None
    }
}

pub(in crate::int_math) fn to_string(value: &BackingInt) -> String {
    ToString::to_string(value)
}

#[cfg(feature = "serde")]
pub(in crate::int_math) fn from_str(value: &str) -> Option<BackingInt> {
    value.parse().ok()
}
