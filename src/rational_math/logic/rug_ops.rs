use crate::int_math::IntType;
use core::cmp::Ordering;
use rug::Rational;
extern crate alloc;
use alloc::format;
use alloc::string::{String, ToString};

pub(in crate::rational_math) type BackingRational = Rational;

pub(in crate::rational_math) fn from_integer(value: IntType) -> BackingRational {
    Rational::from(value)
}

pub(in crate::rational_math) fn new(num: IntType, den: IntType) -> BackingRational {
    Rational::from((num, den))
}

pub(in crate::rational_math) fn numer(value: &BackingRational) -> IntType {
    value.numer().clone()
}

pub(in crate::rational_math) fn denom(value: &BackingRational) -> IntType {
    value.denom().clone()
}

pub(in crate::rational_math) fn add(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    Rational::from(lhs + rhs)
}

pub(in crate::rational_math) fn sub(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    Rational::from(lhs - rhs)
}

pub(in crate::rational_math) fn mul(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    Rational::from(lhs * rhs)
}

pub(in crate::rational_math) fn div(
    lhs: &BackingRational,
    rhs: &BackingRational,
) -> BackingRational {
    Rational::from(lhs / rhs)
}

pub(in crate::rational_math) fn neg(value: &BackingRational) -> BackingRational {
    Rational::from(-value)
}

pub(in crate::rational_math) fn cmp(lhs: &BackingRational, rhs: &BackingRational) -> Ordering {
    lhs.cmp(rhs)
}

pub(in crate::rational_math) const fn is_integer(value: &BackingRational) -> bool {
    value.is_integer()
}

pub(in crate::rational_math) fn to_integer(value: &BackingRational) -> IntType {
    IntType::from(value.numer() / value.denom())
}

pub(in crate::rational_math) fn to_string(value: &BackingRational) -> String {
    if is_integer(value) {
        value.numer().to_string()
    } else {
        format!("{}/{}", value.numer(), value.denom())
    }
}
