#![allow(
    clippy::trivially_copy_pass_by_ref,
    reason = "&T API required for non-Copy backend uniformity"
)]
use crate::int_math::IntType;
use core::cmp::Ordering;
use num_rational::Ratio;
extern crate alloc;
use alloc::format;
use alloc::string::{String, ToString};

pub(in crate::rational_math) type BackingRational = Ratio<IntType>;

#[inline]
pub(in crate::rational_math) fn from_integer(value: IntType) -> BackingRational {
    Ratio::from_integer(value)
}

#[inline]
pub(in crate::rational_math) fn new(num: IntType, den: IntType) -> BackingRational {
    Ratio::new(num, den)
}

#[inline]
pub(in crate::rational_math) const fn numer(value: &BackingRational) -> IntType {
    *value.numer()
}

#[inline]
pub(in crate::rational_math) const fn denom(value: &BackingRational) -> IntType {
    *value.denom()
}

#[inline]
pub(in crate::rational_math) fn add(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    lhs + rhs
}

#[inline]
pub(in crate::rational_math) fn sub(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    lhs - rhs
}

#[inline]
pub(in crate::rational_math) fn mul(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    lhs * rhs
}

#[inline]
pub(in crate::rational_math) fn div(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    lhs / rhs
}

#[inline]
pub(in crate::rational_math) fn neg(value: &BackingRational) -> BackingRational {
    -*value
}

#[inline]
pub(in crate::rational_math) fn cmp(lhs: &BackingRational, rhs: &BackingRational) -> Ordering {
    lhs.cmp(rhs)
}

#[inline]
pub(in crate::rational_math) fn is_integer(value: &BackingRational) -> bool {
    value.is_integer()
}

#[inline]
pub(in crate::rational_math) fn to_integer(value: &BackingRational) -> IntType {
    value.to_integer()
}

#[inline]
pub(in crate::rational_math) fn to_string(value: &BackingRational) -> String {
    if value.is_integer() {
        value.numer().to_string()
    } else {
        format!("{}/{}", value.numer(), value.denom())
    }
}
