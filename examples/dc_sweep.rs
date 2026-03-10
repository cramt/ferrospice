//! DC sweep: I-V curve of a diode.
//!
//! Sweeps voltage source V1 from -0.5V to 0.8V and prints the diode current.

use thevenin::simulate_dc;
use thevenin_types::Netlist;

fn main() {
    let netlist = Netlist::parse(
        "\
Diode IV Curve
V1 anode 0 0
D1 anode 0 DMOD
.model DMOD D IS=1e-14 N=1.0
.dc V1 -0.5 0.8 0.1
.end
",
    )
    .expect("failed to parse netlist");

    let result = simulate_dc(&netlist).expect("simulation failed");

    let plot = &result.plots[0];
    let sweep = &plot.vecs[0]; // sweep variable (v1)
    // Branch current of V1 (current through the voltage source = diode current)
    let current_vec = plot
        .vecs
        .iter()
        .find(|v| v.name.contains("branch"))
        .expect("no branch current found");

    println!("=== Diode I-V Curve ===");
    println!("{:>10}  {:>14}", "V(anode)", "I(diode)");
    println!("{:>10}  {:>14}", "--------", "--------");
    for (v, i) in sweep.real.iter().zip(current_vec.real.iter()) {
        println!("{v:>10.3}  {i:>14.6e}");
    }
}
