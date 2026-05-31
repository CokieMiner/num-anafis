# `num-anafis`

Core numerical and algebraic computation engine of the **SymbAnaFis** library. Provides a highly generic, multi-backend scalar math implementation coupled with a dense **Conformal Geometric Algebra (CGA)** engine powered by the Spectral Mapping Theorem.

---

## Feature Flags

| Feature     | Description                                                              | Default |
|-------------|--------------------------------------------------------------------------|---------|
| `backend64` | Native `i64`/`f64` arithmetic (fastest)                                  | **yes** |
| `backend32` | `i32`/`f32` (memory-optimized)                                           | no      |
| `backendrug`| GMP/MPFR arbitrary-precision via `rug`                                   | no      |
| `clifford`  | Clifford / CGA multivector engine                                        | no      |
| `serde`     | `Serialize`/`Deserialize` for all public types                           | no      |
| `python`    | PyO3 bindings for Python                                                 | no      |
| `std`       | `std` support (enabled automatically by `python`)                        | no      |

Backend features (`backend32`, `backend64`, `backendrug`) are **mutually exclusive**. Exactly one must be active.

---

## Architecture

```
┌──────────────────────────────────────────────────────┐
│                  Clifford Algebra                    │
│  CliffordNumber / FastClifford / GeneratorSet        │
│  Spectral: sin, cos, exp, ln → matrix → eigendecomp  │
├──────────────────────────────────────────────────────┤
│                     Scalar                           │
│  Int → Rational → Float (automatic promotion)        │
├──────────────────────────────────────────────────────┤
│   Backend: float_ops / int_math / rational_math      │
│   i64-f64 | i32-f32 | rug (MPFR/GMP)                 │
└──────────────────────────────────────────────────────┘
```

### 1. Backends

The `IntType`, `FloatType`, and `RationalType` type aliases abstract away precision. Switching from `f64` to 1000-bit MPFR requires **zero changes** to logic — just toggle the feature flag. `Scalar::epsilon()` is dynamically computed from the active backend's mantissa width.

### 2. `Scalar`

A dynamic exact/approximate hybrid number type that automatically promotes through the chain `Int → Rational → Float`:

- `s(42)` — exact integer
- `r(22, 7)` — exact rational
- `s(42).sqrt()` — returns exact `Int` if perfect square, `Float` otherwise
- `s(3).sin()` — falls back to `Float` automatically

Implements the full `Number` trait (trig, exp, log, gamma, zeta, Bessel, etc.).

### 3. Clifford Algebra (`clifford` feature)

Dense multivectors over an arbitrary generator set. Instead of hardcoding separate types for Complex, Quaternions, Dual numbers, or 5D CGA points, **everything** is a single `CliffordNumber` type parameterized by a `GeneratorSet`.

| Type | Description |
|------|-------------|
| `GeneratorSet` | Ordered set of `(id, metric)` pairs defining `Cl(p,q,r)` |
| `CliffordNumber` | Heap-allocated dense multivector (slow path for n>5) |
| `FastClifford<P,Q,R>` | Stack-allocated (`[Scalar; 1<<(P+Q+R)]`), compile-time generic |

#### Convenience constructors (Cl(4,1) CGA)

Pre-built generators for Conformal Geometric Algebra:

| Constructor | Square | Algebra |
|-------------|--------|---------|
| `e1()`, `e2()`, `e3()` | +1 | Euclidean 3D |
| `e_plus()` | +1 | Extra positive |
| `e_minus()` | -1 | Extra negative |
| `ci()` | -1 | Complex unit |
| `sj()` | +1 | Split-complex unit |
| `eps()` | 0 | Dual unit |
| `qi()`, `qj()`, `qk()` | -1 | Quaternion units |
| `orig()` | 0 | CGA origin `eₒ = ½(e₋ − e₊)` |
| `inf()` | 0 | CGA infinity `e∞ = e₋ + e₊` |

#### Spectral Mapping Theorem

The core innovation: any real function `f` (sin, cos, exp, ln, sqrt, etc.) is extended to multivectors via:

1. Map the multivector to a complex matrix (Pauli/Dirac embedding)
2. Compute Schur decomposition → eigenvalues on diagonal
3. Apply `f` to each eigenvalue
4. Reconstruct: `f(A) = V · f(D) · V⁻¹`
5. Map the result matrix back to a multivector

This means `exp(quaternion)`, `sin(CGA point)`, and `sqrt(dual number)` all go through the **same code path**.

#### `FastClifford<P,Q,R>`

Stack-only storage, heap-full compute. The `coeffs` array lives on the stack (no allocation), but all operations internally convert to `CliffordNumber`, compute, and convert back. Useful for small fixed algebras like Cl(3,0), Cl(4,1), or Cl(0,1).

---

## Quick Start

```rust
// backend64 is active by default
use num_anafis::{Scalar, s, r, Number};

let exact: Scalar = s(42);                 // Int
let rational: Scalar = r(22, 7);           // Rational
let approx: Scalar = s(2).sqrt();          // Float (√2 is irrational)
let pi: Scalar = s(3.141592653589793);     // Float
println!("sin(π/4) = {}", (pi / s(4)).sin());
```

```rust
// Add clifford to activate geometric algebra
// Cargo.toml: features = ["backend64", "clifford"]
use num_anafis::{CliffordNumber, Number, e1, e2, e3, ci, cga_gens, s};

let mv = &e1() * &s(3) + &(&e1() * &e2()) * &s(4);
println!("exp(mv)     = {}", mv.exp());
println!("sin(mv)     = {}", mv.sin());
println!("cosh(mv)    = {}", mv.cosh());

// Compare with known formula: exp(θ e₁e₂) = cos(θ) + sin(θ) e₁e₂
let rotor = (&(&e1() * &e2()) * &s(1.0472)).exp(); // exp(θ B), θ = π/3
println!("rotor       = {}", rotor.chop());
```

---

## Future Work

- **Interval arithmetic**: `[lo, hi]` with outward rounding for rigorous error bounds
- **Native Rust Arbitrary precision backend**: Maybe fork dashu and finish it, as it seems abandoned. 
---

## License

Apache License 2.0 — see [LICENSE](LICENSE).
