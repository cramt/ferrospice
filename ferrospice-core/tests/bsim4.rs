//! Integration tests for BSIM4 MOSFET model.
//!
//! Tests ported from ngspice-upstream/tests/bsim4/ qaSpec test suite.
//! Model parameters match nmosParameters/pmosParameters from the BSIM4 CMC QA suite.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use ferrospice_netlist::Netlist;
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
    let result = ferrospice_core::simulate_dc(&netlist).unwrap();

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
    let result = ferrospice_core::simulate_dc(&netlist).unwrap();

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
    let result = ferrospice_core::simulate_dc(&netlist).unwrap();

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
    let result = ferrospice_core::simulate_op(&netlist).unwrap();

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
    let result = ferrospice_core::simulate_op(&netlist).unwrap();

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
    let result = ferrospice_core::simulate_op(&netlist).unwrap();

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
