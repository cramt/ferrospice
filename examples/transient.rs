//! Transient analysis: pulse response of an RC circuit.
//!
//! A 1V pulse drives a 1kΩ + 10µF RC network.
//! Time constant τ = R·C = 10ms.

use thevenin::simulate_tran;
use thevenin_types::Netlist;

fn main() {
    let netlist = Netlist::parse(
        "\
RC Pulse Response
V1 in 0 PULSE(0 1 0 1n 1n 50m 100m)
R1 in out 1k
C1 out 0 10u
.tran 1m 100m
.end
",
    )
    .expect("failed to parse netlist");

    let result = simulate_tran(&netlist).expect("simulation failed");

    let plot = &result.plots[0];
    let time = &plot.vecs[0]; // time
    let vout = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(out)")
        .expect("v(out) not found");

    println!("=== RC Pulse Response ===");
    println!("{:>10}  {:>10}", "Time (s)", "V(out)");
    println!("{:>10}  {:>10}", "--------", "------");

    // Print every 10th point to keep output manageable
    for (i, (t, v)) in time.real.iter().zip(vout.real.iter()).enumerate() {
        if i % 10 == 0 {
            println!("{t:>10.4e}  {v:>10.6}");
        }
    }
}
