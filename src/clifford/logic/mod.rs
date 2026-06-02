mod arithmetic;
mod constructors;
mod core;
mod display;
mod fast;
mod large_mat;
mod matrix;
mod spectral;
mod types;

pub use constructors::{
    cga_gens, ci, e_minus, e_plus, e1, e2, e3, eps, inf, orig, pseudo3d, pseudo5d, qi, qj, qk, sj,
};
pub use fast::FastClifford;
pub use types::CliffordNumber;
pub use types::GeneratorSet;
pub use types::gen_id;
