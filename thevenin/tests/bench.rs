//! Internal benchmarks for cross-target timing comparison (native x86 vs wasm32).
//!
//! Each benchmark times itself using `web_time::Instant` (maps to `performance.now()`
//! on wasm32-unknown-unknown) so that runner/browser startup overhead is excluded.
//!
//! Run with: cargo test --test bench -- --nocapture
//! Or wasm:  cargo test --test bench --target wasm32-unknown-unknown -- --nocapture
//!
//! Output format (parseable by scripts/bench-targets.sh):
//!   BENCH <name> <iterations> <total_ns> <per_iter_ns>

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

use web_time::Instant;

use thevenin::{simulate_ac, simulate_op, simulate_sens};
use thevenin_types::Netlist;

/// Run a closure `iters` times, returning (total_ns, per_iter_ns).
fn bench_fn(iters: u32, mut f: impl FnMut()) -> (u128, u128) {
    // Warmup
    f();

    let start = Instant::now();
    for _ in 0..iters {
        f();
    }
    let total = start.elapsed().as_nanos();
    (total, total / iters as u128)
}

/// Output a BENCH line visible on both native and wasm targets.
fn report(name: &str, iters: u32, total_ns: u128, per_iter_ns: u128) {
    let msg = format!("BENCH {name} {iters} {total_ns} {per_iter_ns}");
    // On wasm32, eprintln goes to console.error which the test runner only shows on failure.
    // Use web_sys::console::log_1 which the runner forwards to stdout.
    #[cfg(target_arch = "wasm32")]
    web_sys::console::log_1(&msg.into());
    #[cfg(not(target_arch = "wasm32"))]
    eprintln!("{msg}");
}

// ---------------------------------------------------------------------------
// Netlist parsing
// ---------------------------------------------------------------------------

const BSIM3_NMOS_PARAMS: &str = "\
+ binunit=1 paramchk=1 mobmod=1 capmod=3 acnqsmod=0 noimod=1 tnom=27
+ nch=1.7e+17 tox=1.5e-08 toxm=1.5e-08
+ wint=0.0 lint=0.0 ll=0.0 wl=0.0 lln=1.0 wln=1.0
+ lw=0.0 ww=0.0 lwn=1.0 wwn=1.0 lwl=0.0 wwl=0.0
+ xpart=0.0 xl=-30e-09
+ vth0=0.7 k1=0.5 k2=0.0 k3=80 k3b=0.0 w0=2.5e-06
+ dvt0=2.2 dvt1=0.53 dvt2=-0.032 nlx=1.74e-07
+ dvt0w=0.0 dvt1w=5.3e6 dvt2w=-0.032 dsub=0.56
+ xj=1.5e-07 ngate=0.0
+ cdsc=2.4e-04 cdscb=0.0 cdscd=0.0 cit=0.0
+ voff=-0.08 nfactor=1.0 eta0=0.08 etab=-0.07
+ vfb=-0.55 u0=670 ua=2.25e-09 ub=5.87e-19 uc=-4.65e-11
+ vsat=8e+04 a0=1.0 ags=0.0 a1=0.0 a2=1.0
+ b0=0.0 b1=0.0 keta=-0.047 dwg=0.0 dwb=0.0
+ pclm=1.3 pdiblc1=0.39 pdiblc2=0.0086 pdiblcb=0.0
+ drout=0.56 pvag=0.0 delta=0.01
+ pscbe1=4.24e+8 pscbe2=1e-05
+ rsh=10.0 rdsw=100.0 prwg=0.0 prwb=0.0 wr=1.0
+ alpha0=0.0 alpha1=0.0 beta0=30
+ cgbo=0.0 cgdl=2e-10 cgsl=2e-10 ckappa=0.6
+ acde=1.0 moin=15 noff=0.9 voffcv=0.02
+ kt1=-0.11 kt1l=0.0 kt2=-0.022 ute=-1.48
+ ua1=4.31e-09 ub1=-7.61e-18 uc1=-5.6e-11 prt=0.0 at=3.3e+04
+ ijth=0.1 js=0.0001 jsw=0.0
+ pb=1.0 cj=0.0005 mj=0.5
+ pbsw=1.0 cjsw=5e-10 mjsw=0.33
+ pbswg=1.0 cjswg=5e-10 mjswg=0.33
+ tpb=0.005 tcj=0.001 tpbsw=0.005 tcjsw=0.001
+ tpbswg=0.005 tcjswg=0.001 xti=3";

fn bsim3_dc_cir(vg: f64) -> String {
    format!(
        "\
BSIM3 DC Bench (Vg={vg})
.model nmod nmos level=8 version=3.3
{BSIM3_NMOS_PARAMS}
M1 d g 0 0 nmod W=10e-6 L=1e-6
Vgs g 0 {vg}
Vds d 0 0.0
Vbs b 0 0.0
.dc Vds 0.1 1.8 0.1
.end
"
    )
}

const LOWPASS_CIR: &str = "RC lowpass bench
V1 1 0 DC 0 AC 1
R1 1 2 1k
C1 2 0 1u
.ac DEC 100 10 100k
.end
";

const RESISTOR_DIVIDER_CIR: &str = "Resistor divider bench
V1 1 0 DC 10
R1 1 2 10k
R2 2 0 10k
.op
.end
";

const SENSITIVITY_CIR: &str = "Sensitivity bench
V1 1 0 DC 10
R1 1 2 1k
R2 2 0 2k
R3 2 3 3k
R4 3 0 4k
.sens V(2)
.end
";

// ---------------------------------------------------------------------------
// Benchmarks
// ---------------------------------------------------------------------------

#[test]
fn bench_parse_bsim3_netlist() {
    let cir = bsim3_dc_cir(1.0);
    let iters = 100;
    let (total, per) = bench_fn(iters, || {
        let _ = Netlist::parse(&cir).unwrap();
    });
    report("parse_bsim3_netlist", iters, total, per);
}

#[test]
fn bench_parse_lowpass_netlist() {
    let iters = 500;
    let (total, per) = bench_fn(iters, || {
        let _ = Netlist::parse(LOWPASS_CIR).unwrap();
    });
    report("parse_lowpass_netlist", iters, total, per);
}

#[test]
fn bench_op_resistor_divider() {
    let netlist = Netlist::parse(RESISTOR_DIVIDER_CIR).unwrap();
    let iters = 100;
    let (total, per) = bench_fn(iters, || {
        let _ = simulate_op(&netlist).unwrap();
    });
    report("op_resistor_divider", iters, total, per);
}

#[test]
fn bench_dc_sweep_bsim3() {
    let cir = bsim3_dc_cir(1.0);
    let netlist = Netlist::parse(&cir).unwrap();
    let iters = 10;
    let (total, per) = bench_fn(iters, || {
        let _ = simulate_op(&netlist).unwrap();
    });
    report("dc_sweep_bsim3", iters, total, per);
}

#[test]
fn bench_ac_lowpass() {
    let netlist = Netlist::parse(LOWPASS_CIR).unwrap();
    let iters = 20;
    let (total, per) = bench_fn(iters, || {
        let _ = simulate_ac(&netlist).unwrap();
    });
    report("ac_lowpass", iters, total, per);
}

#[test]
fn bench_sensitivity() {
    let netlist = Netlist::parse(SENSITIVITY_CIR).unwrap();
    let iters = 50;
    let (total, per) = bench_fn(iters, || {
        let _ = simulate_sens(&netlist).unwrap();
    });
    report("sensitivity", iters, total, per);
}
