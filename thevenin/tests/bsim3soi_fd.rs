//! Integration tests for BSIM3SOI-FD MOSFET model (level 55).
//!
//! Tests ported from ngspice-upstream/tests/bsim3soifd/ test suite.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use thevenin::{simulate_dc, simulate_op};
use thevenin_types::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// NMOS BSIM3SOI-FD model parameters from ngspice-upstream/tests/bsim3soifd/nmosfd.mod
const NMOS_SOI_FD_PARAMS: &str = "\
+ TNOM=27 TOX=4.5E-09 TSI=5e-8 TBOX=8E-08
+ MOBMOD=0 CAPMOD=3 SHMOD=0
+ PARAMCHK=0 WINT=0 LINT=-2E-08
+ VTH0=.52 K1=.39 K2=.1 K3=0
+ KB1=.95 K3B=2.2 NLX=7.2E-08
+ DVT0=.55 DVT1=.28 DVT2=-1.4
+ DVT0W=0 DVT1W=0 DVT2W=0
+ NCH=3.3E+17 NSUB=1E+15 NGATE=1E+20
+ DVBD0=60.0 DVBD1=1.1 VBSA=0.0
+ KB3=2.2 DELP=0.02
+ ABP=0.9 MXC=0.9 ADICE0=0.93
+ KBJT1=1.0E-08 EDL=.0000005
+ NDIODE=1.13 NTUN=14.0
+ ISBJT=2e-6 ISDIF=1e-6 ISTUN=0.0 ISREC=1e-5
+ XBJT=0.01 XDIF=0.01 XREC=0.01 XTUN=0.001
+ U0=352 UA=1.3E-11 UB=1.7E-18 UC=-4E-10
+ W0=1.16E-06 AGS=.25 A1=0 A2=1
+ B0=.01 B1=10
+ RDSW=700 PRWG=0 PRWB=-.2 WR=1
+ RBODY=0.0 RBSH=0.0
+ A0=1.4 KETA=-.67 VSAT=135000
+ DWG=0 DWB=0
+ ALPHA0=0.0 ALPHA1=1.5 BETA0=20.5
+ AII=1.2 BII=0.1e-7 CII=0.8 DII=0.6
+ VOFF=-.14 NFACTOR=.7 CDSC=.00002 CDSCB=0
+ CDSCD=0 CIT=0
+ PCLM=2.9 PVAG=12 PDIBLC1=.18 PDIBLC2=.004
+ PDIBLCB=-.234 DROUT=.2
+ DELTA=.01 ETA0=.01 ETAB=0
+ DSUB=.3 RTH0=.006
+ CLC=.0000001 CLE=.6 CF=1E-20 CKAPPA=.6
+ CGDL=1E-20 CGSL=1E-20 KT1=-.3 KT1L=0
+ KT2=.022 UTE=-1.5 UA1=4.31E-09 UB1=-7.61E-18
+ UC1=-5.6E-11 PRT=760 AT=22400
+ CGSO=1e-10 CGDO=1e-10 CJSWG=5e-10 TT=3e-10
+ ASD=0.3 CSDESW=1e-12";

/// PMOS BSIM3SOI-FD model parameters from ngspice-upstream/tests/bsim3soifd/pmosfd.mod
const PMOS_SOI_FD_PARAMS: &str = "\
+ TNOM=27 TOX=4.5E-09 TSI=5e-8 TBOX=8E-08
+ MOBMOD=2 CAPMOD=3 SHMOD=0
+ WINT=0 LINT=-2E-08
+ VTH0=-.52 K1=.39 K2=.1 K3=0
+ KB1=.95 K3B=2.2 NLX=7.2E-08
+ DVT0=.55 DVT1=.28 DVT2=-1.4
+ DVT0W=0 DVT1W=0 DVT2W=0
+ NCH=3.0E+17 NSUB=1E+15 NGATE=1E+20
+ DVBD0=60.0 DVBD1=1.1 VBSA=-0.2
+ KB3=2.2 DELP=0.02
+ ABP=0.9 MXC=0.9 ADICE0=0.93
+ KBJT1=1.0E-08 EDL=.0000005
+ NDIODE=1.13 NTUN=14.0
+ ISBJT=0.0 ISDIF=1e-6 ISTUN=0.0 ISREC=0.0
+ XBJT=0.01 XDIF=0.01 XREC=0.01 XTUN=0.001
+ U0=145 UA=1.3E-11 UB=1.7E-18 UC=-4E-10
+ W0=1.16E-06 AGS=.25 A1=0 A2=1
+ B0=.01 B1=10
+ RDSW=700 PRWG=0 PRWB=-.2 WR=1
+ RBODY=0.0 RBSH=0.0
+ A0=1.4 KETA=-.67 VSAT=75000
+ DWG=0 DWB=0
+ ALPHA0=0.0 ALPHA1=1e-4 BETA0=19
+ AII=1.25 BII=0.1e-7 CII=0.8 DII=0.6
+ VOFF=-.14 NFACTOR=.7 CDSC=.00002 CDSCB=0
+ CDSCD=0 CIT=0
+ PCLM=2.9 PVAG=12 PDIBLC1=.18 PDIBLC2=.004
+ PDIBLCB=-.234 DROUT=.2
+ DELTA=.01 ETA0=.01 ETAB=0
+ DSUB=.3 RTH0=.006
+ CLC=.0000001 CLE=.6 CF=1E-20 CKAPPA=.6
+ CGDL=1E-20 CGSL=1E-20 KT1=-.3 KT1L=0
+ KT2=.022 UTE=-1.5 UA1=4.31E-09 UB1=-7.61E-18
+ UC1=-5.6E-11 PRT=760 AT=22400
+ CGSO=1e-10 CGDO=1e-10 CJSWG=5e-10 TT=3e-10
+ ASD=0.3 CSDESW=1e-12";

/// Test that a single NMOS FD SOI transistor can reach a DC operating point.
/// Uses the floating body configuration (4 terminals: d g s e).
#[test]
fn bsim3soi_fd_nmos_op() {
    let cir = format!(
        "\
BSIM3SOI-FD NMOS OP Test
.model n1 nmos level=55
{NMOS_SOI_FD_PARAMS}
m1 d g 0 e n1 W=10u L=0.25u
Vgs g 0 1.0
Vds d 0 1.5
Ve e 0 0.0
.op
.end
"
    );
    let netlist = Netlist::parse(&cir).expect("parse");
    let result = simulate_op(&netlist);
    assert!(result.is_ok(), "OP failed: {:?}", result.err());
}

/// Test NMOS FD SOI with different bias points to verify operation regions.
#[test]
fn bsim3soi_fd_nmos_bias_points() {
    for (vgs, vds) in [(0.5, 0.1), (1.0, 1.0), (1.0, 1.5), (1.5, 1.5)] {
        let cir = format!(
            "\
BSIM3SOI-FD NMOS Bias Vgs={vgs} Vds={vds}
.model n1 nmos level=55
{NMOS_SOI_FD_PARAMS}
m1 d g 0 e n1 W=10u L=0.25u
Vgs g 0 {vgs}
Vds d 0 {vds}
Ve e 0 0.0
.op
.end
"
        );
        let netlist = Netlist::parse(&cir).expect("parse");
        let result = simulate_op(&netlist);
        assert!(
            result.is_ok(),
            "OP failed at Vgs={vgs}, Vds={vds}: {:?}",
            result.err()
        );
    }
}

/// Test PMOS FD SOI operating point.
#[test]
fn bsim3soi_fd_pmos_op() {
    let cir = format!(
        "\
BSIM3SOI-FD PMOS OP Test
.model p1 pmos level=55
{PMOS_SOI_FD_PARAMS}
m1 d g vdd e p1 W=20u L=0.25u
Vgs g vdd 0
Vds d vdd -1.5
Vdd vdd 0 2.5
Ve e 0 1.25
.op
.end
"
    );
    let netlist = Netlist::parse(&cir).expect("parse");
    let result = simulate_op(&netlist);
    assert!(result.is_ok(), "PMOS OP failed: {:?}", result.err());
}

/// SOI FD Inverter test from ngspice-upstream/tests/bsim3soifd/inv2.cir pattern.
/// Tests both NMOS and PMOS FD SOI together in a CMOS inverter.
#[test]
#[ignore = "SOI FD inverter convergence needs source-stepping improvements"]
fn bsim3soi_fd_inverter_op() {
    let cir = format!(
        "\
SOI FD Inverter - floating body
.model n1 nmos level=55
{NMOS_SOI_FD_PARAMS}
.model p1 pmos level=55
{PMOS_SOI_FD_PARAMS}
vin in 0 dc 0.0
vdd dd 0 dc 2.5
vss ss 0 dc 0
ve e 0 dc 1.25
m1 out in dd e p1 W=20u L=0.25u
m2 out in ss e n1 W=10u L=0.25u
.op
.end
"
    );
    let netlist = Netlist::parse(&cir).expect("parse");
    let result = simulate_op(&netlist);
    assert!(result.is_ok(), "Inverter OP failed: {:?}", result.err());

    let op = result.unwrap();
    let plot = &op.plots[0];
    let v_out = plot
        .vecs
        .iter()
        .find(|sv| sv.name == "v(out)")
        .expect("v(out)");
    assert!(
        v_out.real[0] > 2.0,
        "Inverter output should be high with vin=0, got {}",
        v_out.real[0]
    );
}

/// Test 5-terminal FD SOI MOSFET (with body contact).
#[test]
fn bsim3soi_fd_5_terminal() {
    let cir = format!(
        "\
BSIM3SOI-FD 5-terminal NMOS
.model n1 nmos level=55
{NMOS_SOI_FD_PARAMS}
m1 d g 0 e b n1 W=10u L=0.25u NBC=1
Vgs g 0 1.0
Vds d 0 1.5
Ve e 0 0.0
Vb b 0 0.0
.op
.end
"
    );
    let netlist = Netlist::parse(&cir).expect("parse");
    let result = simulate_op(&netlist);
    assert!(result.is_ok(), "5-terminal OP failed: {:?}", result.err());
}

/// Test DC sweep for FD SOI MOSFET.
#[test]
#[ignore = "SOI FD DC sweep convergence needs source-stepping improvements"]
fn bsim3soi_fd_dc_sweep() {
    let cir = format!(
        "\
SOI FD NMOS DC Sweep
.model n1 nmos level=55
{NMOS_SOI_FD_PARAMS}
m1 d g 0 e n1 W=10u L=0.25u
Vgs g 0 1.5
Vds d 0 0.0
Ve e 0 0.0
.dc Vds 0.1 2.5 0.1
.end
"
    );
    let netlist = Netlist::parse(&cir).expect("parse");
    let result = simulate_dc(&netlist);
    assert!(result.is_ok(), "DC sweep failed: {:?}", result.err());

    let plot = &result.unwrap().plots[0];
    let vsweep = plot
        .vecs
        .iter()
        .find(|v| v.name == "v-sweep")
        .expect("v-sweep");
    assert_eq!(vsweep.real.len(), 25, "Expected 25 sweep points");
}
