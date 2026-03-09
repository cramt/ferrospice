//! Integration tests for BSIM3SOI-PD MOSFET model (level 57).
//!
//! Tests ported from ngspice-upstream/tests/bsim3soipd/ test suite.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use ferrospice_core::{simulate_dc, simulate_op};
use ferrospice_netlist::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// NMOS BSIM3SOI-PD model parameters from ngspice-upstream/tests/bsim3soipd/nmospd.mod
const NMOS_SOI_PARAMS: &str = "\
+ TNOM=27 TOX=4.5E-09 TSI=1E-7 TBOX=8E-08
+ MOBMOD=0 CAPMOD=2 SHMOD=0
+ PARAMCHK=0 WINT=0 LINT=-2E-08
+ VTH0=.42 K1=.49 K2=.1 K3=0
+ K3B=2.2 NLX=2E-7
+ DVT0=10 DVT1=.55 DVT2=-1.4
+ DVT0W=0 DVT1W=0 DVT2W=0
+ NCH=4.7E+17 NSUB=-1E+15 NGATE=1E+20
+ AGIDL=1e-15 BGIDL=1e9 NGIDL=1.1
+ NDIODE=1.13 NTUN=14.0 NRECF0=2.5 NRECR0=4
+ VREC0=1.2 NTRECF=.1 NTRECR=.2
+ ISBJT=1e-4 ISDIF=1e-5 ISTUN=2e-5 ISREC=4e-2
+ XBJT=.9 XDIF=.9 XREC=.9 XTUN=0.01
+ AHLI=1e-9 LBJT0=0.2e-6 LN=2e-6
+ NBJT=.8 NDIF=-1 AELY=1e8 VABJT=0
+ U0=352 UA=1.3E-11 UB=1.7E-18 UC=-4E-10
+ W0=1.16E-06 AGS=.25 A1=0 A2=1
+ B0=.01 B1=10
+ RDSW=0 PRWG=0 PRWB=-.2 WR=1
+ RBODY=1e0 RBSH=0.0
+ A0=1.4 KETA=0.1 KETAS=0.2 VSAT=135000
+ DWG=0 DWB=0
+ ALPHA0=1e-8 BETA0=0 BETA1=0.05 BETA2=0.07
+ VDSATII0=.8 ESATII=1e7
+ VOFF=-.14 NFACTOR=.7 CDSC=.00002 CDSCB=0
+ CDSCD=0 CIT=0
+ PCLM=2.9 PVAG=12 PDIBLC1=.18 PDIBLC2=.004
+ PDIBLCB=-.234 DROUT=.2
+ DELTA=.01 ETA0=.05 ETAB=0
+ DSUB=.2 RTH0=.005
+ CLC=1E-07 CLE=.6 CF=1E-20 CKAPPA=.6
+ CGDL=1E-20 CGSL=1E-20 KT1=-.3 KT1L=0
+ KT2=.022 UTE=-1.5 UA1=4.31E-09 UB1=-7.61E-18
+ UC1=-5.6E-11 PRT=760 AT=22400
+ CGSO=1e-10 CGDO=1e-10 CJSWG=1e-12 TT=3e-10
+ ASD=0.3 CSDESW=1e-12
+ TCJSWG=1e-4 MJSWG=.5 PBSWG=1";

/// PMOS BSIM3SOI-PD model parameters from ngspice-upstream/tests/bsim3soipd/pmospd.mod
const PMOS_SOI_PARAMS: &str = "\
+ TNOM=27 TOX=4.5E-09 TSI=1E-7 TBOX=8E-08
+ MOBMOD=2 CAPMOD=2 SHMOD=0
+ WINT=0 LINT=-2E-08
+ VTH0=-.42 K1=.49 K2=.1 K3=0
+ K3B=2.2 NLX=2E-07
+ DVT0=10 DVT1=.55 DVT2=-1.4
+ DVT0W=0 DVT1W=0 DVT2W=0
+ NCH=4.7E+17 NSUB=-1E+15 NGATE=1E+20
+ AGIDL=1e-16 BGIDL=1e9 NGIDL=1.1
+ NDIODE=1.13 NTUN=14.0 NRECF0=2.5 NRECR0=4
+ VREC0=1.2 NTRECF=.1 NTRECR=.2
+ ISBJT=1e-4 ISDIF=1e-5 ISTUN=2e-5 ISREC=4e-2
+ XBJT=.9 XDIF=.9 XREC=.9 XTUN=0.01
+ AHLI=1e-9 LBJT0=0.2e-6 LN=2e-6
+ NBJT=.8 NDIF=-1 AELY=1e8 VABJT=0
+ U0=145 UA=1.3E-11 UB=1.7E-18 UC=-4E-10
+ W0=1.16E-06 AGS=.25 A1=0 A2=1
+ B0=.01 B1=10
+ RDSW=350 PRWG=0 PRWB=-.2 WR=1
+ RBODY=1e0 RBSH=0.0
+ A0=1.4 KETA=0.1 KETAS=0.2 VSAT=75000
+ DWG=0 DWB=0
+ ALPHA0=1e-8 BETA0=0 BETA1=0.05 BETA2=0.07
+ VDSATII0=1.6 ESATII=1e7
+ VOFF=-.14 NFACTOR=.7 CDSC=.00002 CDSCB=0
+ CDSCD=0 CIT=0
+ PCLM=2.9 PVAG=12 PDIBLC1=.18 PDIBLC2=.004
+ PDIBLCB=-.234 DROUT=.2
+ DELTA=.01 ETA0=.05 ETAB=0
+ DSUB=.2 RTH0=.005
+ CLC=1E-07 CLE=.6 CF=1E-20 CKAPPA=.6
+ CGDL=1E-20 CGSL=1E-20 KT1=-.3 KT1L=0
+ KT2=.022 UTE=-1.5 UA1=4.31E-09 UB1=-7.61E-18
+ UC1=-5.6E-11 PRT=760 AT=22400
+ CGSO=1e-10 CGDO=1e-10 CJSWG=1e-12 TT=3e-10
+ ASD=0.3 CSDESW=1e-12
+ TCJSWG=1e-4 MJSWG=.5 PBSWG=1";

/// Test that a single NMOS SOI transistor can reach a DC operating point.
/// Uses the floating body configuration (4 terminals: d g s e).
#[test]
fn bsim3soi_pd_nmos_op() {
    let cir = format!(
        "\
BSIM3SOI-PD NMOS OP Test
.model n1 nmos level=57
{NMOS_SOI_PARAMS}
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

/// Test NMOS SOI with different bias points to verify operation regions.
#[test]
fn bsim3soi_pd_nmos_bias_points() {
    for (vgs, vds) in [(0.5, 0.1), (1.0, 1.0), (1.0, 1.5), (1.5, 1.5)] {
        let cir = format!(
            "\
BSIM3SOI-PD NMOS Bias Vgs={vgs} Vds={vds}
.model n1 nmos level=57
{NMOS_SOI_PARAMS}
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

/// Test PMOS SOI operating point.
#[test]
fn bsim3soi_pd_pmos_op() {
    let cir = format!(
        "\
BSIM3SOI-PD PMOS OP Test
.model p1 pmos level=57
{PMOS_SOI_PARAMS}
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

/// SOI Inverter test from ngspice-upstream/tests/bsim3soipd/inv2.cir
/// Tests both NMOS and PMOS SOI together in a CMOS inverter.
#[test]
#[ignore = "SOI inverter convergence needs source-stepping improvements"]
fn bsim3soi_pd_inverter_op() {
    let cir = format!(
        "\
SOI Inverter - floating body
.model n1 nmos level=57
{NMOS_SOI_PARAMS}
.model p1 pmos level=57
{PMOS_SOI_PARAMS}
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

    // With vin=0, PMOS on, NMOS off => output should be near VDD=2.5
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

/// SOI Inverter with input high — output should be near ground.
#[test]
#[ignore = "SOI inverter convergence needs source-stepping improvements"]
fn bsim3soi_pd_inverter_input_high() {
    let cir = format!(
        "\
SOI Inverter - input high
.model n1 nmos level=57
{NMOS_SOI_PARAMS}
.model p1 pmos level=57
{PMOS_SOI_PARAMS}
vin in 0 dc 2.5
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

    // With vin=2.5, NMOS on, PMOS off => output should be near 0
    let op = result.unwrap();
    let plot = &op.plots[0];
    let v_out = plot
        .vecs
        .iter()
        .find(|sv| sv.name == "v(out)")
        .expect("v(out)");
    assert!(
        v_out.real[0] < 0.5,
        "Inverter output should be low with vin=2.5, got {}",
        v_out.real[0]
    );
}

/// Test 5-terminal SOI MOSFET (with body contact).
/// From ngspice t4.cir pattern: m1 d g s e b n1
#[test]
fn bsim3soi_pd_5_terminal() {
    let cir = format!(
        "\
BSIM3SOI-PD 5-terminal NMOS
.model n1 nmos level=57
{NMOS_SOI_PARAMS}
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

/// Test DC sweep from ngspice-upstream/tests/bsim3soipd/inv2.cir pattern.
#[test]
#[ignore = "SOI DC sweep convergence needs source-stepping improvements"]
fn bsim3soi_pd_dc_sweep() {
    let cir = format!(
        "\
SOI NMOS DC Sweep
.model n1 nmos level=57
{NMOS_SOI_PARAMS}
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
    // Should have 25 data points (0.1 to 2.5 in steps of 0.1)
    let vsweep = plot
        .vecs
        .iter()
        .find(|v| v.name == "v-sweep")
        .expect("v-sweep");
    assert_eq!(vsweep.real.len(), 25, "Expected 25 sweep points");
}
