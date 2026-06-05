//! Backend resolution: determines which numeric backend to compile.
//!
//! Priority:
//!   1. If exactly one backend feature is enabled → use it.
//!   2. If none is enabled → fall back to backend64 (i64/f64).
//!   3. If more than one is enabled → set all of them so `compile_error!` in
//!      `lib.rs` fires with a helpful message.
//!
//! This emits a custom `cfg(backend = "32"|"64"|"rug")` that the rest of the
//! crate uses for module selection, instead of raw `cfg(feature = "...")`.

fn main() {
    let b32 = std::env::var("CARGO_FEATURE_BACKEND32").is_ok();
    let b64 = std::env::var("CARGO_FEATURE_BACKEND64").is_ok();
    let brug = std::env::var("CARGO_FEATURE_BACKENDRUG").is_ok();

    let count = u8::from(b32) + u8::from(b64) + u8::from(brug);

    if count > 1 {
        // Multiple backends selected — set all so compile_error! triggers.
        if b32 {
            println!("cargo::rustc-cfg=backend=\"32\"");
        }
        if b64 {
            println!("cargo::rustc-cfg=backend=\"64\"");
        }
        if brug {
            println!("cargo::rustc-cfg=backend=\"rug\"");
        }
    } else if brug {
        println!("cargo::rustc-cfg=backend=\"rug\"");
    } else if b32 {
        println!("cargo::rustc-cfg=backend=\"32\"");
    } else {
        // No backend selected (or backend64 explicitly) → default to 64-bit.
        println!("cargo::rustc-cfg=backend=\"64\"");
    }

    // Declare valid cfg values so check-cfg (Rust 1.80+) doesn't warn.
    println!("cargo::rustc-check-cfg=cfg(backend, values(\"32\", \"64\", \"rug\"))");
}
