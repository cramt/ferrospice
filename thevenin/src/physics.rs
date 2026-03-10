//! Shared physical constants and utility functions for device models.
//!
//! Consolidates constants (EPSOX, EPSSI, CHARGE_Q, KBOQ, etc.) and safe exponential
//! functions that were previously duplicated across 10+ device model files.

// ─── Physical constants ────────────────────────────────────────────────────────

/// Permittivity of SiO₂ (F/m): ε₀ × 3.9
pub const EPSOX: f64 = 3.453133e-11;

/// Permittivity of Si (F/m): ε₀ × 11.7
pub const EPSSI: f64 = 1.03594e-10;

/// Elementary charge (C)
pub const CHARGE_Q: f64 = 1.60219e-19;

/// Boltzmann constant / elementary charge (V/K): k_B / q
pub const KBOQ: f64 = 8.617087e-5;

/// Silicon bandgap at 300 K (eV)
pub const EG300: f64 = 1.115;

// ─── Safe exponential (simple) ─────────────────────────────────────────────────

/// Upper clamp for the simple safe_exp (used by diode, BJT, JFET, MOSFET, MOS6, VBIC).
pub const EXP_LIMIT: f64 = 500.0;

/// Clamped exponential: prevents overflow for x > 500.
/// Used by simple device models (diode, BJT, JFET, MOSFET, MOS6, VBIC).
#[inline]
pub fn safe_exp(x: f64) -> f64 {
    if x > EXP_LIMIT {
        EXP_LIMIT.exp()
    } else {
        x.exp()
    }
}

// ─── Safe exponential (BSIM) ──────────────────────────────────────────────────

/// Upper clamp threshold for BSIM safe_exp.
pub const EXP_THRESHOLD: f64 = 34.0;

/// exp(34) — precomputed upper clamp value.
pub const MAX_EXP: f64 = 5.834617425e14;

/// exp(-34) — precomputed lower clamp value.
pub const MIN_EXP: f64 = 1.713908431e-15;

/// Clamped exponential for BSIM models: clamps both high (>34) and low (<-34).
/// Returns precomputed MAX_EXP/MIN_EXP to avoid overflow/underflow.
#[inline]
pub fn bsim_safe_exp(x: f64) -> f64 {
    if x > EXP_THRESHOLD {
        MAX_EXP
    } else if x < -EXP_THRESHOLD {
        MIN_EXP
    } else {
        x.exp()
    }
}

// ─── Safe exponential (SOI — extended range) ───────────────────────────────────

/// Upper clamp threshold for SOI extended exponential.
pub const EXPL_THRESHOLD: f64 = 100.0;

/// exp(100) — precomputed upper clamp value.
pub const MAX_EXPL: f64 = 2.688117142e43;

/// exp(-100) — precomputed lower clamp value.
pub const MIN_EXPL: f64 = 3.720075976e-44;

/// Safe exponential for SOI models (larger threshold than BSIM3).
/// Returns `(exp_val, derivative)` where derivative equals exp_val when in range,
/// MAX_EXPL at saturation, 0 at underflow.
#[inline]
pub fn soi_dexp(x: f64) -> (f64, f64) {
    if x > EXPL_THRESHOLD {
        (MAX_EXPL * (1.0 + x - EXPL_THRESHOLD), MAX_EXPL)
    } else if x < -EXPL_THRESHOLD {
        (MIN_EXPL, 0.0)
    } else {
        let e = x.exp();
        (e, e)
    }
}
