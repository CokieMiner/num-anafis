#![allow(
    clippy::missing_const_for_fn,
    clippy::trivially_copy_pass_by_ref,
    reason = "Delegation wrappers can't be const for all backends; &T API for non-Copy uniformity"
)]

use alloc::string::String;
use core::cmp::Ordering;

// Backend selection priority:
// 1. i64        (default)
// 2. backendrug   (rug::Integer — GMP-based arbitrary precision)
// 3. backend32     (i32 — memory-optimized)

use super::logic::backend;
use super::logic::backend::BackingInt;
#[cfg(feature = "serde")]
use super::logic::backend::from_str as backend_from_str;

/// The integer representation type selected by the active backend feature.
pub type IntType = BackingInt;

/// Macro that generates delegation functions for the integer backend.
/// Every function listed here MUST be implemented by every backend module,
/// or compilation will fail with a clear missing-function error.
macro_rules! delegate_int_ops {
    (
        $(
            fn $name:ident( $($arg:ident : $ty:ty),* ) -> $ret:ty;
        )*
    ) => {
        $(
            #[inline]
            pub fn $name( $($arg : $ty),* ) -> $ret {
                backend::$name( $($arg),* )
            }
        )*
    };
}

delegate_int_ops! {
    fn zero() -> IntType;
    fn clone(value: &IntType) -> IntType;
    fn from_i64(value: i64) -> Option<IntType>;
    fn abs(value: &IntType) -> Option<IntType>;
    fn neg(value: &IntType) -> Option<IntType>;
    fn add(lhs: &IntType, rhs: &IntType) -> Option<IntType>;
    fn sub(lhs: &IntType, rhs: &IntType) -> Option<IntType>;
    fn mul(lhs: &IntType, rhs: &IntType) -> Option<IntType>;
    fn cmp(lhs: &IntType, rhs: &IntType) -> Ordering;
    fn is_zero(value: &IntType) -> bool;
    fn is_one(value: &IntType) -> bool;
    fn is_neg_one(value: &IntType) -> bool;
    fn is_negative(value: &IntType) -> bool;
    fn is_positive(value: &IntType) -> bool;
    fn is_even(value: &IntType) -> bool;
    fn modulo(lhs: &IntType, rhs: &IntType) -> IntType;
    fn div_exact(lhs: &IntType, rhs: &IntType) -> IntType;
    fn gcd(lhs: &IntType, rhs: &IntType) -> IntType;
    fn perfect_square(value: &IntType) -> Option<IntType>;
    fn perfect_cube(value: &IntType) -> Option<IntType>;
    fn to_string(value: &IntType) -> String;
}

#[cfg(feature = "serde")]
#[inline]
pub fn from_str(value: &str) -> Option<IntType> {
    backend_from_str(value)
}
