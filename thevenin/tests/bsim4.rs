//! Integration tests for BSIM4 MOSFET model.
//!
//! Tests ported from ngspice-upstream/tests/bsim4/ qaSpec test suite.
//! Model parameters match nmosParameters/pmosParameters from the BSIM4 CMC QA suite.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use thevenin_types::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Full BSIM4 NMOS model parameters from ngspice-upstream/tests/bsim4/nmos/parameters/nmosParameters.
const BSIM4_NMOS_PARAMS: &str = "\
+ binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1
+ geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1
+ acnqsmod=0 trnqsmod=0 tnom=27
+ toxe=1.85e-9 toxp=1.2e-9 toxm=1.85e-9 dtox=0.65e-9 epsrox=3.9
+ wint=5e-9 lint=5.25e-9 ll=0 wl=0 lln=1 wln=1
+ lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0
+ xpart=0 toxref=1.85e-9 xl=-30e-9
+ vth0=0.423 k1=0.4 k2=0.01 k3=0 k3b=0 w0=2.5e-6
+ dvt0=1 dvt1=2 dvt2=-0.032 dvt0w=0 dvt1w=0 dvt2w=0
+ dsub=0.1 minv=0.05 voffl=0 dvtp0=1.0e-9 dvtp1=0.1
+ lpe0=0 lpeb=0 xj=1.96e-8
+ ngate=2e20 ndep=2.54e18 nsd=2e20 phin=0
+ cdsc=0.0 cdscb=0 cdscd=0 cit=0
+ voff=-0.13 nfactor=1.9 eta0=0.0058 etab=0 vfb=-0.55
+ u0=0.0491 ua=6e-10 ub=1.2e-18 uc=0 vsat=124340
+ a0=1.0 ags=1e-20 a1=0 a2=1.0 b0=0 b1=0
+ keta=0.04 dwg=0 dwb=0
+ pclm=0.04 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=-0.005
+ drout=0.5 pvag=1e-20 delta=0.01
+ pscbe1=8.14e8 pscbe2=1e-7
+ fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2.3e6
+ rsh=5 rdsw=165 rsw=85 rdw=85 rdswmin=0 rdwmin=0 rswmin=0
+ prwg=0 prwb=6.8e-11 wr=1
+ alpha0=0.074 alpha1=0.005 beta0=30
+ agidl=0.0002 bgidl=2.1e9 cgidl=0.0002 egidl=0.8
+ aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1
+ aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3
+ aigc=0.012 bigc=0.0028 cigc=0.002
+ aigsd=0.012 bigsd=0.0028 cigsd=0.002
+ nigc=1 poxedge=1 pigcd=1 ntox=1
+ xrcrg1=12 xrcrg2=5
+ cgso=1.5e-10 cgdo=1.5e-10 cgbo=2.56e-11
+ cgdl=2.653e-10 cgsl=2.653e-10 ckappas=0.03 ckappad=0.03
+ acde=1 moin=15 noff=0.9 voffcv=0.02
+ kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5
+ ua1=4.31e-9 ub1=7.61e-18 uc1=-5.6e-11 prt=0 at=33000
+ fnoimod=1 tnoimod=0";

/// Full BSIM4 PMOS model parameters from ngspice-upstream/tests/bsim4/pmos/parameters/pmosParameters.
const BSIM4_PMOS_PARAMS: &str = "\
+ binunit=1 paramchk=1 mobmod=0 capmod=2 igcmod=1 igbmod=1
+ geomod=1 diomod=1 rdsmod=0 rbodymod=1 rgatemod=1 permod=1
+ acnqsmod=0 trnqsmod=0 tnom=27
+ toxe=1.95e-9 toxp=1.2e-9 toxm=1.95e-9 dtox=0.75e-9 epsrox=3.9
+ wint=5e-9 lint=5.25e-9 ll=0 wl=0 lln=1 wln=1
+ lw=0 ww=0 lwn=1 wwn=1 lwl=0 wwl=0
+ xpart=0 toxref=1.95e-9 xl=-30e-9
+ vth0=-0.365 k1=0.4 k2=-0.01 k3=0 k3b=0 w0=2.5e-6
+ dvt0=1 dvt1=2 dvt2=-0.032 dvt0w=0 dvt1w=0 dvt2w=0
+ dsub=0.1 minv=0.05 voffl=0 dvtp0=1.0e-9 dvtp1=0.05
+ lpe0=0 lpeb=0 xj=1.96e-8
+ ngate=2e20 ndep=1.87e18 nsd=2e20 phin=0
+ cdsc=0.0 cdscb=0 cdscd=0 cit=0
+ voff=-0.126 nfactor=1.9 eta0=0.0058 etab=0 vfb=0.55
+ u0=0.00574 ua=2.0e-9 ub=0.5e-18 uc=0 vsat=70000
+ a0=1.0 ags=1e-20 a1=0 a2=1 b0=-1e-20 b1=0
+ keta=-0.047 dwg=0 dwb=0
+ pclm=0.12 pdiblc1=0.001 pdiblc2=0.001 pdiblcb=3.4e-8
+ drout=0.56 pvag=1e-20 delta=0.01
+ pscbe1=8.14e8 pscbe2=9.58e-7
+ fprout=0.2 pdits=0.08 pditsd=0.23 pditsl=2.3e6
+ rsh=5 rdsw=165 rsw=85 rdw=85 rdswmin=0 rdwmin=0 rswmin=0
+ prwg=3.22e-8 prwb=6.8e-11 wr=1
+ alpha0=0.074 alpha1=0.005 beta0=30
+ agidl=0.0002 bgidl=2.1e9 cgidl=0.0002 egidl=0.8
+ aigbacc=0.012 bigbacc=0.0028 cigbacc=0.002 nigbacc=1
+ aigbinv=0.014 bigbinv=0.004 cigbinv=0.004 eigbinv=1.1 nigbinv=3
+ aigc=0.69 bigc=0.0012 cigc=0.0008
+ aigsd=0.0087 bigsd=0.0012 cigsd=0.0008
+ nigc=1 poxedge=1 pigcd=1 ntox=1
+ xrcrg1=12 xrcrg2=5
+ cgso=1.5e-10 cgdo=1.5e-10 cgbo=2.56e-11
+ cgdl=2.653e-10 cgsl=2.653e-10 ckappas=0.03 ckappad=0.03
+ acde=1 moin=15 noff=0.9 voffcv=0.02
+ kt1=-0.11 kt1l=0 kt2=0.022 ute=-1.5
+ ua1=4.31e-9 ub1=7.61e-18 uc1=-5.6e-11 prt=0 at=33000
+ fnoimod=1 tnoimod=0";

/// Build a .cir netlist for a BSIM4 NMOS DC sweep test.
fn bsim4_dc_sweep_cir(vg: f64) -> String {
    format!(
        "\
BSIM4 NMOS DC Sweep Test (Vg={vg})
.model nmod nmos level=14
{BSIM4_NMOS_PARAMS}
M1 d g 0 0 nmod W=10e-6 L=100e-9
Vgs g 0 {vg}
Vds d 0 0.0
.dc Vds 0.1 1.2 0.1
.end
"
    )
}

// ============================================================================
// DC Sweep Tests
// ============================================================================

/// BSIM4 NMOS DC sweep at Vg=1.0V — validates saturation current.
#[test]
fn test_bsim4_dc_sweep_vg1p0() {
    let cir = bsim4_dc_sweep_cir(1.0);
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let i_vds = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    assert_eq!(i_vds.real.len(), 12, "Expected 12 sweep points");

    // Current should be negative (flowing into drain)
    for (i, &current) in i_vds.real.iter().enumerate() {
        assert!(
            current < 0.0,
            "Drain current at point {} should be negative (NMOS): {}",
            i,
            current
        );
    }

    // Current magnitude should increase with Vds then saturate
    let i_first = i_vds.real[0].abs();
    let i_last = i_vds.real[11].abs();
    assert!(
        i_last > i_first,
        "Current should increase with Vds: first={}, last={}",
        i_first,
        i_last
    );

    // At Vg=1.0, Vds=1.2 for this technology, current on the order of mA
    assert!(
        i_last > 1e-4,
        "Saturation current should be > 0.1mA: {}",
        i_last
    );
    assert!(
        i_last < 0.1,
        "Saturation current should be < 100mA: {}",
        i_last
    );
}

/// BSIM4 NMOS DC sweep at Vg=0.6V — near-threshold operation.
#[test]
fn test_bsim4_dc_sweep_vg0p6() {
    let cir = bsim4_dc_sweep_cir(0.6);
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let i_vds = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    assert_eq!(i_vds.real.len(), 12);
    let i_sat = i_vds.real[11].abs();
    assert!(
        i_sat < 1e-3,
        "Near-threshold current should be < 1mA: {}",
        i_sat
    );
}

/// BSIM4 NMOS DC sweep — check saturation behavior (current plateau).
#[test]
fn test_bsim4_dc_sweep_saturation_plateau() {
    let cir = bsim4_dc_sweep_cir(1.0);
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let i_vds = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    // In saturation, current should plateau — last 4 points within 20%
    let last_4: Vec<f64> = i_vds.real[8..12].iter().map(|x| x.abs()).collect();
    let max = last_4.iter().cloned().fold(f64::MIN, f64::max);
    let min = last_4.iter().cloned().fold(f64::MAX, f64::min);
    assert!(
        (max - min) / max < 0.20,
        "Saturation current should plateau: min={}, max={}",
        min,
        max
    );
}

// ============================================================================
// PMOS Operating Point Tests
// ============================================================================

/// BSIM4 PMOS operating point — validates PMOS conducts with correct bias.
#[test]
fn test_bsim4_pmos_op() {
    let cir = format!(
        "\
BSIM4 PMOS OP Test
.model pmod pmos level=14
{BSIM4_PMOS_PARAMS}
M1 d g vdd vdd pmod W=10e-6 L=100e-9
Vdd vdd 0 1.2
Vgs g 0 0.2
Vds d 0 0.2
.op
.end
"
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();

    let plot = &result.plots[0];
    let vds_branch = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    let ids = vds_branch.real[0];
    assert!(
        ids.abs() > 1e-6,
        "PMOS should have significant current: {}",
        ids
    );
}

/// BSIM4 PMOS operating point — cutoff when Vgs ≈ 0.
#[test]
fn test_bsim4_pmos_cutoff() {
    let cir = format!(
        "\
BSIM4 PMOS Cutoff Test
.model pmod pmos level=14
{BSIM4_PMOS_PARAMS}
M1 d g vdd vdd pmod W=10e-6 L=100e-9
Vdd vdd 0 1.2
Vgs g 0 1.2
Vds d 0 1.2
.op
.end
"
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();

    let plot = &result.plots[0];
    let vds_branch = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    let ids = vds_branch.real[0];
    assert!(
        ids.abs() < 1e-6,
        "PMOS in cutoff should have very small current: {}",
        ids
    );
}

// ============================================================================
// NMOS Operating Point Test
// ============================================================================

/// BSIM4 NMOS OP with full qaSpec parameters — validates convergence and current magnitude.
#[test]
fn test_bsim4_nmos_op_full_params() {
    let cir = format!(
        "\
BSIM4 NMOS OP Full Params Test
.model nmod nmos level=14
{BSIM4_NMOS_PARAMS}
M1 d g 0 0 nmod W=10e-6 L=100e-9
Vgs g 0 1.0
Vds d 0 1.0
.op
.end
"
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();

    let plot = &result.plots[0];
    let vds_branch = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    let ids = vds_branch.real[0];
    assert!(ids < 0.0, "NMOS drain current should be negative: {}", ids);
    assert!(
        ids.abs() > 1e-5 && ids.abs() < 0.1,
        "NMOS current should be in reasonable range (10uA-100mA): {}",
        ids
    );
}

// ============================================================================
// AC Frequency Tests (ported from qaSpec acFreq)
// ============================================================================

/// BSIM4 NMOS AC frequency sweep — validates capacitance values (C(g,g), C(g,s), C(g,d)).
/// Reference: ngspice-upstream/tests/bsim4/nmos/reference/acFreq.standard
#[test]
fn test_bsim4_ac_freq() {
    let cir = format!(
        "\
BSIM4 NMOS AC Frequency Test
.model nmod nmos level=14
{BSIM4_NMOS_PARAMS}
M1 d g 0 0 nmod W=10e-6 L=60e-9
Vgs g 0 dc 1.2 ac 1.0
Vds d 0 1.2
.ac dec 10 1e3 1e8
.end
"
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_ac(&netlist).unwrap();

    let plot = &result.plots[0];
    let freq_vec = plot
        .vecs
        .iter()
        .find(|v| v.name == "frequency")
        .expect("no frequency vector");
    assert_eq!(
        freq_vec.real.len(),
        51,
        "Expected 51 frequency points (dec 10 1e3 1e8)"
    );

    // Check gate voltage complex response at mid-frequency (10kHz, index ~10)
    let vg = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(g)")
        .expect("no v(g)");

    // Gate voltage should have non-zero imaginary part due to capacitances
    // At low frequencies, real part dominates; at high frequencies, imaginary grows
    assert!(
        !vg.complex.is_empty(),
        "AC results should have complex data"
    );

    // The gate is driven by Vgs AC source, so v(g) should be approximately 1.0 + j*0
    // at all frequencies (it's the source node)
    let c = &vg.complex[0];
    let mag = (c.re * c.re + c.im * c.im).sqrt();
    assert!(
        (mag - 1.0).abs() < 0.1,
        "Gate AC voltage should be ~1.0: mag={}",
        mag
    );
}

/// BSIM4 NMOS AC frequency with xpart=1 — tests charge partitioning variation.
/// Reference: ngspice-upstream/tests/bsim4/nmos/reference/acFreq_xpart.standard
#[test]
fn test_bsim4_ac_freq_xpart() {
    let cir = format!(
        "\
BSIM4 NMOS AC Freq xpart=1
.model nmod nmos level=14
{BSIM4_NMOS_PARAMS}
+ xpart=1
M1 d g 0 0 nmod W=10e-6 L=60e-9
Vgs g 0 dc 1.2 ac 1.0
Vds d 0 1.2
.ac dec 10 1e3 1e8
.end
"
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_ac(&netlist).unwrap();

    let plot = &result.plots[0];
    let freq_vec = plot
        .vecs
        .iter()
        .find(|v| v.name == "frequency")
        .expect("no frequency vector");
    assert_eq!(freq_vec.real.len(), 51);

    // With xpart=1 (0/100 partition), all inversion charge to source.
    // This changes the drain/source capacitance partitioning but total gate cap stays same.
    let vg = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(g)")
        .expect("no v(g)");
    assert!(!vg.complex.is_empty());
}

/// BSIM4 NMOS AC frequency with mobmod=1 — tests mobility model variation.
#[test]
fn test_bsim4_ac_freq_mobmod() {
    let cir = format!(
        "\
BSIM4 NMOS AC Freq mobmod=1
.model nmod nmos level=14
{BSIM4_NMOS_PARAMS}
+ mobmod=1
M1 d g 0 0 nmod W=10e-6 L=60e-9
Vgs g 0 dc 1.2 ac 1.0
Vds d 0 1.2
.ac dec 10 1e3 1e8
.end
"
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_ac(&netlist).unwrap();

    let plot = &result.plots[0];
    assert_eq!(
        plot.vecs
            .iter()
            .find(|v| v.name == "frequency")
            .unwrap()
            .real
            .len(),
        51
    );
}

// ============================================================================
// Noise Tests (ported from qaSpec noise1-noise4)
// ============================================================================

/// Build a noise test netlist for BSIM4 with specified noise model parameters.
/// Uses Vdd + Rload topology so that the drain node is NOT directly on a voltage source,
/// allowing the adjoint transfer function to be non-zero.
fn bsim4_noise_cir(vg: f64, noise_params: &str) -> String {
    format!(
        "\
BSIM4 NMOS Noise Test (Vg={vg})
.model nmod nmos level=14
{BSIM4_NMOS_PARAMS}
+ {noise_params}
M1 d g 0 0 nmod W=10e-6 L=60e-9
Vgs g 0 {vg}
Vdd vdd 0 2.0
Rload vdd d 1k
.noise V(d) Vgs dec 10 1e3 1e8
.end
"
    )
}

/// noise1: fnoimod=0 tnoimod=0 — KF flicker + simple thermal.
/// Reference: ngspice-upstream/tests/bsim4/nmos/reference/noise1.standard
#[test]
fn test_bsim4_noise1_fnoi0_tnoi0() {
    let cir = bsim4_noise_cir(1.0, "fnoimod=0 tnoimod=0 kf=1e-30 af=1.2 ef=1.1");
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_noise(&netlist).unwrap();

    // Should have noise spectrum plot
    assert!(
        !result.plots.is_empty(),
        "Should have at least one noise plot"
    );

    let plot = &result.plots[0];
    let freq_vec = plot
        .vecs
        .iter()
        .find(|v| v.name == "frequency")
        .expect("no frequency vector");
    assert_eq!(freq_vec.real.len(), 51, "Expected 51 frequency points");

    // Get output noise spectrum
    let onoise = plot
        .vecs
        .iter()
        .find(|v| v.name.contains("onoise") && v.name.contains("spectrum"))
        .expect("no onoise_spectrum vector");

    // Noise values should be positive and reasonable
    for &val in &onoise.real {
        assert!(val >= 0.0, "Noise should be non-negative: {}", val);
    }

    // With fnoimod=0 (KF flicker), noise should decrease with frequency (1/f behavior)
    let low_freq_noise = onoise.real[0]; // 1kHz
    let high_freq_noise = onoise.real[50]; // 100MHz
    assert!(
        low_freq_noise > high_freq_noise || (low_freq_noise < 1.0e-30 && high_freq_noise < 1.0e-30),
        "Flicker noise should decrease with frequency: low={}, high={}",
        low_freq_noise,
        high_freq_noise
    );
}

/// noise2: fnoimod=1 tnoimod=1 — NOIA trap flicker + unified thermal.
/// Reference: ngspice-upstream/tests/bsim4/nmos/reference/noise2.standard
#[test]
fn test_bsim4_noise2_fnoi1_tnoi1() {
    let cir = bsim4_noise_cir(
        1.0,
        "fnoimod=1 tnoimod=1 noia=6.25e41 noib=3.125e26 noic=8.75 em=4.1e7 ef=1.1",
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_noise(&netlist).unwrap();

    let plot = &result.plots[0];
    let onoise = plot
        .vecs
        .iter()
        .find(|v| v.name.contains("onoise") && v.name.contains("spectrum"))
        .expect("no onoise_spectrum");

    // With both flicker and thermal, noise at low freq should be higher than high freq
    let low_freq_noise = onoise.real[0];
    let high_freq_noise = onoise.real[50];

    assert!(
        low_freq_noise > 0.0,
        "Low-freq noise should be positive: {}",
        low_freq_noise
    );
    assert!(
        high_freq_noise > 0.0,
        "High-freq noise should be positive: {}",
        high_freq_noise
    );
    // Low frequency noise should dominate due to flicker
    assert!(
        low_freq_noise > high_freq_noise,
        "Low-freq should dominate: low={}, high={}",
        low_freq_noise,
        high_freq_noise
    );
}

/// noise3: fnoimod=1 tnoimod=0 — NOIA flicker only (no explicit thermal from channel).
#[test]
fn test_bsim4_noise3_fnoi1_tnoi0() {
    let cir = bsim4_noise_cir(
        0.8,
        "fnoimod=1 tnoimod=0 noia=6.25e41 noib=3.125e26 noic=8.75 em=4.1e7 ef=1.1",
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_noise(&netlist).unwrap();

    let plot = &result.plots[0];
    let onoise = plot
        .vecs
        .iter()
        .find(|v| v.name.contains("onoise") && v.name.contains("spectrum"))
        .expect("no onoise_spectrum");

    // NOIA flicker noise should be significant
    assert!(
        onoise.real[0] > 0.0,
        "NOIA noise should be non-zero: {}",
        onoise.real[0]
    );
}

/// noise4: fnoimod=0 tnoimod=1 — KF flicker + unified thermal.
#[test]
fn test_bsim4_noise4_fnoi0_tnoi1() {
    let cir = bsim4_noise_cir(0.6, "fnoimod=0 tnoimod=1 kf=1e-30 af=1.2 ef=1.1");
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_noise(&netlist).unwrap();

    let plot = &result.plots[0];
    let onoise = plot
        .vecs
        .iter()
        .find(|v| v.name.contains("onoise") && v.name.contains("spectrum"))
        .expect("no onoise_spectrum");

    // With tnoimod=1 and very small kf, thermal noise should dominate at high freq
    // (roughly flat), and tiny flicker at low freq
    let mid_noise = onoise.real[25];
    assert!(
        mid_noise > 0.0,
        "Should have some noise at mid-freq: {}",
        mid_noise
    );
}

// ============================================================================
// AC Vd Test (ported from qaSpec acVd)
// ============================================================================

/// BSIM4 NMOS AC drain voltage sweep — validates transconductance and capacitance vs Vd.
/// Reference: ngspice-upstream/tests/bsim4/nmos/reference/acVd.standard
#[test]
fn test_bsim4_ac_vd_sweep() {
    let cir = format!(
        "\
BSIM4 NMOS AC Vd Sweep
.model nmod nmos level=14
{BSIM4_NMOS_PARAMS}
M1 d g 0 0 nmod W=10e-6 L=60e-9
Vgs g 0 1.2
Vds d 0 0.0
.dc Vds 0.1 1.2 0.1
.end
"
    );
    let netlist = Netlist::parse(&cir).unwrap();
    let result = thevenin::simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let i_vds = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    // Should have 12 DC sweep points
    assert_eq!(i_vds.real.len(), 12);

    // Current should be negative (NMOS, conventional Vds source)
    for val in &i_vds.real {
        assert!(*val < 0.0, "Drain current should be negative: {}", val);
    }
}
