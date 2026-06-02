#![allow(
    clippy::trivially_copy_pass_by_ref,
    reason = "&T API required for non-Copy backend uniformity"
)]

use alloc::string::{String, ToString};
use core::cmp::Ordering;

pub(in crate::int_math) type BackingInt = i32;

#[inline]
pub(in crate::int_math) fn from_i64(value: i64) -> Option<BackingInt> {
    i32::try_from(value).ok()
}

impl_shared_int_ops!();
