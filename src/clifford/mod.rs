//! Dense multivectors for Clifford (geometric) algebras with labeled generators.
//!
//! This module implements **Conformal Geometric Algebra** CGA Cl(4,1) with 5
//! base generators: `e1, e2, e3` (Euclidean 3-space), `e+` (extra positive),
//! and `e−` (extra negative).  All derived objects (quaternions, complex unit,
//! conformal origin/infinity) are constructed as products of these generators.

#![allow(
    clippy::indexing_slicing,
    clippy::module_name_repetitions,
    reason = "Indexing is bounds-checked by construction; type name repeats module name for discoverability."
)]

mod api;
mod logic;

pub use api::*;
