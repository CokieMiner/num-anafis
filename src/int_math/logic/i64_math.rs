#![allow(
    clippy::trivially_copy_pass_by_ref,
    clippy::unnecessary_wraps,
    reason = "&T API for non-Copy backend uniformity; Option<IntType> for uniform API"
)]

use alloc::string::{String, ToString};
use core::cmp::Ordering;
pub(in crate::int_math) type BackingInt = i64;

#[inline]
pub(in crate::int_math) const fn from_i64(value: i64) -> Option<BackingInt> {
    Some(value)
}

impl_shared_int_ops!();
