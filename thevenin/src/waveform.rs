//! Transient waveform evaluation for independent voltage/current sources.
//!
//! Implements PULSE, SIN, EXP, PWL, SFFM, and AM waveforms matching ngspice
//! semantics (see `ngspice-upstream/src/spicelib/devices/vsrc/vsrcload.c`).

use std::f64::consts::PI;

use thevenin_types::{Expr, Waveform};

/// Parameters needed for waveform default value computation.
#[derive(Debug, Clone, Copy)]
pub struct TranParams {
    /// The timestep from .tran (used as default for rise/fall times).
    pub tstep: f64,
    /// The final simulation time (used as default for SIN freq, SFFM freqs).
    pub tstop: f64,
}

/// Evaluate a waveform at a given time.
///
/// Returns the instantaneous value of the source at time `t`.
/// Uses `tran` parameters for default values when waveform parameters are omitted.
pub fn evaluate(wf: &Waveform, t: f64, tran: &TranParams) -> f64 {
    match wf {
        Waveform::Pulse {
            v1,
            v2,
            td,
            tr,
            tf,
            pw,
            per,
        } => eval_pulse(
            val(v1),
            val(v2),
            opt(td).unwrap_or(0.0),
            opt(tr).unwrap_or(tran.tstep).max(tran.tstep),
            opt(tf).unwrap_or(tran.tstep).max(tran.tstep),
            opt(pw).unwrap_or(0.0).max(0.0),
            opt(per),
            t,
            tran.tstep,
        ),
        Waveform::Sin {
            v0,
            va,
            freq,
            td,
            theta,
            phi,
        } => eval_sin(
            val(v0),
            val(va),
            opt(freq).unwrap_or(if tran.tstop > 0.0 {
                1.0 / tran.tstop
            } else {
                1.0
            }),
            opt(td).unwrap_or(0.0),
            opt(theta).unwrap_or(0.0),
            opt(phi).unwrap_or(0.0),
            t,
        ),
        Waveform::Exp {
            v1,
            v2,
            td1,
            tau1,
            td2,
            tau2,
        } => eval_exp(
            val(v1),
            val(v2),
            opt(td1).unwrap_or(tran.tstep),
            opt(tau1).unwrap_or(tran.tstep).max(tran.tstep),
            opt(td2).unwrap_or(opt(td1).unwrap_or(tran.tstep) + tran.tstep),
            opt(tau2).unwrap_or(tran.tstep).max(tran.tstep),
            t,
        ),
        Waveform::Pwl(points) => {
            let pairs: Vec<(f64, f64)> = points
                .iter()
                .map(|p| (val(&p.time), val(&p.value)))
                .collect();
            eval_pwl(&pairs, t)
        }
        Waveform::Sffm { v0, va, fc, fs, md } => {
            let fc_val = opt(fc).unwrap_or(if tran.tstop > 0.0 {
                5.0 / tran.tstop
            } else {
                5.0
            });
            let fs_val = opt(fs).unwrap_or(if tran.tstop > 0.0 {
                500.0 / tran.tstop
            } else {
                500.0
            });
            // Clamp modulation depth to [0, fc/fs].
            let max_md = if fs_val > 0.0 { fc_val / fs_val } else { 0.0 };
            let md_val = opt(md).unwrap_or(0.0).clamp(0.0, max_md);
            eval_sffm(val(v0), val(va), fc_val, fs_val, md_val, t)
        }
        Waveform::Am { va, vo, fc, fs, td } => eval_am(
            val(va),
            val(vo),
            val(fc),
            val(fs),
            opt(td).unwrap_or(0.0),
            t,
        ),
    }
}

/// Get the DC value of a waveform (value at t=0 for DC operating point).
///
/// For most waveforms, this is the initial/offset value. Returns `None` if
/// the waveform's DC contribution should come from the Source's `dc` field instead.
pub fn dc_value(wf: &Waveform) -> Option<f64> {
    match wf {
        Waveform::Pulse { v1, .. } => Some(val(v1)),
        Waveform::Sin { v0, va, phi, .. } => {
            let phi_rad = opt(phi).unwrap_or(0.0) * PI / 180.0;
            Some(val(v0) + val(va) * phi_rad.sin())
        }
        Waveform::Exp { v1, .. } => Some(val(v1)),
        Waveform::Pwl(points) => {
            if points.is_empty() {
                None
            } else {
                Some(val(&points[0].value))
            }
        }
        Waveform::Sffm { v0, .. } => Some(val(v0)),
        Waveform::Am { .. } => Some(0.0),
    }
}

// ---- Individual waveform evaluators ----

/// PULSE waveform.
///
/// Parameters: v1, v2, td, tr, tf, pw, per
/// Before td: v1
/// Rise: td to td+tr (linear from v1 to v2)
/// High: td+tr to td+tr+pw (v2)
/// Fall: td+tr+pw to td+tr+pw+tf (linear from v2 to v1)
/// After: v1 (repeats if periodic)
#[expect(clippy::too_many_arguments)]
fn eval_pulse(
    v1: f64,
    v2: f64,
    td: f64,
    tr: f64,
    tf: f64,
    pw: f64,
    per: Option<f64>,
    t: f64,
    tstep: f64,
) -> f64 {
    if t < td {
        return v1;
    }

    // Compute the period.
    let period = per.unwrap_or(tr + pw + tf).max(tr + pw + tf).max(tstep);

    // Fold time into the current period.
    let mut time = t - td;
    if period > 0.0 && time >= period {
        time -= period * (time / period).floor();
    }

    if time < tr {
        // Rising edge.
        v1 + (v2 - v1) * time / tr
    } else if time < tr + pw {
        // Pulse width (high).
        v2
    } else if time < tr + pw + tf {
        // Falling edge.
        v2 + (v1 - v2) * (time - tr - pw) / tf
    } else {
        // Rest of period (low).
        v1
    }
}

/// SIN waveform.
///
/// Before td: v0 + va * sin(phi)
/// After td: v0 + va * sin(2*pi*freq*(t-td) + phi) * exp(-theta*(t-td))
fn eval_sin(v0: f64, va: f64, freq: f64, td: f64, theta: f64, phi_deg: f64, t: f64) -> f64 {
    let phi_rad = phi_deg * PI / 180.0;

    if t <= td {
        v0 + va * phi_rad.sin()
    } else {
        let dt = t - td;
        let damping = if theta != 0.0 {
            (-theta * dt).exp()
        } else {
            1.0
        };
        v0 + va * (2.0 * PI * freq * dt + phi_rad).sin() * damping
    }
}

/// EXP waveform.
///
/// t <= td1: v1
/// td1 < t <= td2: v1 + (v2-v1)*(1 - exp(-(t-td1)/tau1))
/// t > td2: above + (v1-v2)*(1 - exp(-(t-td2)/tau2))
fn eval_exp(v1: f64, v2: f64, td1: f64, tau1: f64, td2: f64, tau2: f64, t: f64) -> f64 {
    if t <= td1 {
        v1
    } else if t <= td2 {
        v1 + (v2 - v1) * (1.0 - (-(t - td1) / tau1).exp())
    } else {
        v1 + (v2 - v1) * (1.0 - (-(t - td1) / tau1).exp())
            + (v1 - v2) * (1.0 - (-(t - td2) / tau2).exp())
    }
}

/// PWL (piecewise linear) waveform.
///
/// Linear interpolation between time-value pairs.
/// Before first point: first value. After last point: last value.
fn eval_pwl(points: &[(f64, f64)], t: f64) -> f64 {
    if points.is_empty() {
        return 0.0;
    }

    // Before first point.
    if t <= points[0].0 {
        return points[0].1;
    }

    // After last point — hold last value.
    if t >= points[points.len() - 1].0 {
        return points[points.len() - 1].1;
    }

    // Find the segment containing t.
    for i in 1..points.len() {
        if t <= points[i].0 {
            let (t0, v0) = points[i - 1];
            let (t1, v1) = points[i];
            let dt = t1 - t0;
            if dt <= 0.0 {
                return v1;
            }
            return v0 + (v1 - v0) * (t - t0) / dt;
        }
    }

    points[points.len() - 1].1
}

/// SFFM (single-frequency FM) waveform.
///
/// v0 + va * sin(2*pi*fc*t + md*sin(2*pi*fs*t))
fn eval_sffm(v0: f64, va: f64, fc: f64, fs: f64, md: f64, t: f64) -> f64 {
    v0 + va * (2.0 * PI * fc * t + md * (2.0 * PI * fs * t).sin()).sin()
}

/// AM (amplitude modulation) waveform.
///
/// Before td: 0
/// After td: va * (vo + sin(2*pi*fs*t)) * sin(2*pi*fc*t)
///
/// Note: ngspice AM parameters are (VO, VMO, VMA, FM, FC, TD, ...) but
/// ferrospice-netlist parses as (va, vo, fc, fs, td) — mapping:
/// va = amplitude, vo = offset, fc = carrier freq, fs = signal freq
fn eval_am(va: f64, vo: f64, fc: f64, fs: f64, td: f64, t: f64) -> f64 {
    if t < td {
        0.0
    } else {
        va * (vo + (2.0 * PI * fs * t).sin()) * (2.0 * PI * fc * t).sin()
    }
}

/// Collect breakpoint times from a waveform within [0, tstop].
///
/// Breakpoints are times where the waveform has a discontinuity in value or
/// first derivative (e.g., PULSE edges, PWL corners, EXP transitions).
/// The transient engine uses these to force timestep boundaries at these points.
pub fn breakpoints(wf: &Waveform, tran: &TranParams) -> Vec<f64> {
    let tstop = tran.tstop;
    let mut bps = Vec::new();

    match wf {
        Waveform::Pulse {
            td,
            tr,
            tf,
            pw,
            per,
            ..
        } => {
            let td_val = opt(td).unwrap_or(0.0);
            let tr_val = opt(tr).unwrap_or(tran.tstep).max(tran.tstep);
            let tf_val = opt(tf).unwrap_or(tran.tstep).max(tran.tstep);
            let pw_val = opt(pw).unwrap_or(0.0).max(0.0);
            let period = per
                .as_ref()
                .map(val)
                .unwrap_or(tr_val + pw_val + tf_val)
                .max(tr_val + pw_val + tf_val)
                .max(tran.tstep);

            // Breakpoints within one period relative to td:
            // td, td+tr, td+tr+pw, td+tr+pw+tf
            let edges = [0.0, tr_val, tr_val + pw_val, tr_val + pw_val + tf_val];

            let mut k = 0u64;
            loop {
                let base = td_val + k as f64 * period;
                if base > tstop {
                    break;
                }
                for &edge in &edges {
                    let t = base + edge;
                    if t >= 0.0 && t <= tstop {
                        bps.push(t);
                    }
                }
                k += 1;
                // Safety limit to avoid infinite loop with tiny periods.
                if k > 1_000_000 {
                    break;
                }
            }
        }
        Waveform::Pwl(points) => {
            for p in points {
                let t = val(&p.time);
                if t >= 0.0 && t <= tstop {
                    bps.push(t);
                }
            }
        }
        Waveform::Exp { td1, td2, .. } => {
            let td1_val = opt(td1).unwrap_or(tran.tstep);
            let td2_val = opt(td2).unwrap_or(td1_val + tran.tstep);
            if td1_val >= 0.0 && td1_val <= tstop {
                bps.push(td1_val);
            }
            if td2_val >= 0.0 && td2_val <= tstop {
                bps.push(td2_val);
            }
        }
        Waveform::Sin { td, .. } => {
            let td_val = opt(td).unwrap_or(0.0);
            if td_val > 0.0 && td_val <= tstop {
                bps.push(td_val);
            }
        }
        Waveform::Am { td, .. } => {
            let td_val = opt(td).unwrap_or(0.0);
            if td_val > 0.0 && td_val <= tstop {
                bps.push(td_val);
            }
        }
        Waveform::Sffm { .. } => {
            // SFFM is smooth everywhere, no breakpoints.
        }
    }

    bps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    bps.dedup_by(|a, b| (*a - *b).abs() < 1e-18);
    bps
}

// ---- Helpers ----

/// Extract a numeric value from an Expr. Panics on non-numeric (param exprs not supported yet).
fn val(expr: &Expr) -> f64 {
    match expr {
        Expr::Num(v) => *v,
        _ => 0.0, // Fallback for unresolved params
    }
}

/// Extract an optional numeric value.
fn opt(expr: &Option<Expr>) -> Option<f64> {
    expr.as_ref().map(val)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    fn tran() -> TranParams {
        TranParams {
            tstep: 1e-6,
            tstop: 1e-3,
        }
    }

    #[test]
    fn test_pulse_basic() {
        let wf = Waveform::Pulse {
            v1: Expr::Num(0.0),
            v2: Expr::Num(5.0),
            td: Some(Expr::Num(0.0)),
            tr: Some(Expr::Num(1e-6)),
            tf: Some(Expr::Num(1e-6)),
            pw: Some(Expr::Num(10e-6)),
            per: Some(Expr::Num(20e-6)),
        };
        let tp = tran();

        // Before delay: v1
        assert_eq!(evaluate(&wf, -1e-6, &tp), 0.0);

        // Middle of rise
        let v = evaluate(&wf, 0.5e-6, &tp);
        assert!((v - 2.5).abs() < 0.01, "mid-rise: {v}");

        // At pulse high
        assert_eq!(evaluate(&wf, 5e-6, &tp), 5.0);

        // At pulse end (start of fall)
        let v = evaluate(&wf, 11.5e-6, &tp);
        assert!((v - 2.5).abs() < 0.01, "mid-fall: {v}");

        // After fall (low)
        assert_eq!(evaluate(&wf, 15e-6, &tp), 0.0);

        // Next period: high again
        assert_eq!(evaluate(&wf, 25e-6, &tp), 5.0);
    }

    #[test]
    fn test_sin_basic() {
        let wf = Waveform::Sin {
            v0: Expr::Num(0.0),
            va: Expr::Num(1.0),
            freq: Some(Expr::Num(1000.0)),
            td: Some(Expr::Num(0.0)),
            theta: Some(Expr::Num(0.0)),
            phi: Some(Expr::Num(0.0)),
        };
        let tp = tran();

        // At t=0: sin(0) = 0
        assert!((evaluate(&wf, 0.0, &tp)).abs() < 1e-10);

        // At t=0.25ms: sin(pi/2) = 1
        let v = evaluate(&wf, 0.25e-3, &tp);
        assert!((v - 1.0).abs() < 1e-10, "quarter period: {v}");

        // At t=0.5ms: sin(pi) = 0
        let v = evaluate(&wf, 0.5e-3, &tp);
        assert!(v.abs() < 1e-10, "half period: {v}");

        // At t=0.75ms: sin(3*pi/2) = -1
        let v = evaluate(&wf, 0.75e-3, &tp);
        assert!((v + 1.0).abs() < 1e-10, "3/4 period: {v}");
    }

    #[test]
    fn test_sin_with_offset_and_phase() {
        let wf = Waveform::Sin {
            v0: Expr::Num(2.5),
            va: Expr::Num(1.0),
            freq: Some(Expr::Num(1000.0)),
            td: None,
            theta: None,
            phi: Some(Expr::Num(90.0)), // 90 degrees
        };
        let tp = tran();

        // At t=0 (before td=0): v0 + va*sin(90°) = 2.5 + 1.0 = 3.5
        let v = evaluate(&wf, 0.0, &tp);
        assert!((v - 3.5).abs() < 1e-10, "t=0 with phase: {v}");
    }

    #[test]
    fn test_sin_with_damping() {
        let wf = Waveform::Sin {
            v0: Expr::Num(0.0),
            va: Expr::Num(1.0),
            freq: Some(Expr::Num(1000.0)),
            td: Some(Expr::Num(0.0)),
            theta: Some(Expr::Num(1000.0)), // Heavy damping
            phi: None,
        };
        let tp = tran();

        // At t=0.25ms: sin(pi/2)*exp(-1000*0.00025) = 1 * exp(-0.25)
        let v = evaluate(&wf, 0.25e-3, &tp);
        let expected = (-0.25_f64).exp();
        assert!((v - expected).abs() < 1e-6, "damped: {v} vs {expected}");
    }

    #[test]
    fn test_exp_basic() {
        let wf = Waveform::Exp {
            v1: Expr::Num(0.0),
            v2: Expr::Num(5.0),
            td1: Some(Expr::Num(0.0)),
            tau1: Some(Expr::Num(1e-3)),
            td2: Some(Expr::Num(5e-3)),
            tau2: Some(Expr::Num(1e-3)),
        };
        let tp = tran();

        // At t=0: v1
        assert_eq!(evaluate(&wf, 0.0, &tp), 0.0);

        // At t=1ms: v1 + (v2-v1)*(1-exp(-1)) = 5*(1-0.368) ≈ 3.16
        let v = evaluate(&wf, 1e-3, &tp);
        let expected = 5.0 * (1.0 - (-1.0_f64).exp());
        assert!((v - expected).abs() < 1e-6, "rise: {v} vs {expected}");

        // At t=6ms (1ms after td2): should be decaying back
        let v = evaluate(&wf, 6e-3, &tp);
        assert!(v < 5.0, "should be decaying: {v}");
    }

    #[test]
    fn test_pwl_basic() {
        let wf = Waveform::Pwl(vec![
            thevenin_types::PwlPoint {
                time: Expr::Num(0.0),
                value: Expr::Num(0.0),
            },
            thevenin_types::PwlPoint {
                time: Expr::Num(1e-3),
                value: Expr::Num(5.0),
            },
            thevenin_types::PwlPoint {
                time: Expr::Num(2e-3),
                value: Expr::Num(5.0),
            },
            thevenin_types::PwlPoint {
                time: Expr::Num(3e-3),
                value: Expr::Num(0.0),
            },
        ]);
        let tp = tran();

        assert_eq!(evaluate(&wf, 0.0, &tp), 0.0);
        assert!((evaluate(&wf, 0.5e-3, &tp) - 2.5).abs() < 1e-10);
        assert_eq!(evaluate(&wf, 1e-3, &tp), 5.0);
        assert_eq!(evaluate(&wf, 1.5e-3, &tp), 5.0);
        assert!((evaluate(&wf, 2.5e-3, &tp) - 2.5).abs() < 1e-10);
        assert_eq!(evaluate(&wf, 3e-3, &tp), 0.0);
        // After last point: hold
        assert_eq!(evaluate(&wf, 5e-3, &tp), 0.0);
    }

    #[test]
    fn test_sffm_basic() {
        let wf = Waveform::Sffm {
            v0: Expr::Num(1.0),
            va: Expr::Num(2.0),
            fc: Some(Expr::Num(1000.0)),
            fs: Some(Expr::Num(100.0)),
            md: Some(Expr::Num(5.0)),
        };
        let tp = tran();

        // At t=0: v0 + va*sin(0 + md*sin(0)) = 1 + 2*sin(0) = 1.0
        let v = evaluate(&wf, 0.0, &tp);
        assert!((v - 1.0).abs() < 1e-10, "sffm at t=0: {v}");
    }

    #[test]
    fn test_am_basic() {
        let wf = Waveform::Am {
            va: Expr::Num(1.0),
            vo: Expr::Num(1.0),
            fc: Expr::Num(10000.0),
            fs: Expr::Num(1000.0),
            td: Some(Expr::Num(0.0)),
        };
        let tp = tran();

        // At t=0: va*(vo + sin(0))*sin(0) = 0
        let v = evaluate(&wf, 0.0, &tp);
        assert!(v.abs() < 1e-10, "am at t=0: {v}");
    }

    #[test]
    fn test_sin_matches_analytical() {
        // Acceptance criteria: SIN source output samples match analytical sin() within 1e-6
        let freq = 1000.0;
        let va = 3.3;
        let v0 = 1.65;
        let wf = Waveform::Sin {
            v0: Expr::Num(v0),
            va: Expr::Num(va),
            freq: Some(Expr::Num(freq)),
            td: Some(Expr::Num(0.0)),
            theta: Some(Expr::Num(0.0)),
            phi: Some(Expr::Num(0.0)),
        };
        let tp = tran();

        // Check at many time points
        for i in 0..100 {
            let t = i as f64 * 1e-5; // 0 to 1ms
            let expected = v0 + va * (2.0 * PI * freq * t).sin();
            let actual = evaluate(&wf, t, &tp);
            assert!(
                (actual - expected).abs() < 1e-6,
                "at t={t:.6e}: expected={expected:.10}, got={actual:.10}"
            );
        }
    }
}
