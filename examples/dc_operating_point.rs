//! DC operating point analysis of a simple resistive voltage divider.
//!
//! Circuit:
//!   V1 (1V) -> R1 (1kΩ) -> node "mid" -> R2 (2kΩ) -> GND
//!
//! Expected: V(mid) = 1V * 2k/(1k+2k) ≈ 0.667V

use ferrospice_core::simulate_op;
use ferrospice_netlist::Netlist;

fn main() {
    let netlist = Netlist::parse(
        "\
Voltage Divider
V1 in 0 1.0
R1 in mid 1k
R2 mid 0 2k
.op
.end
",
    )
    .expect("failed to parse netlist");

    let result = simulate_op(&netlist).expect("simulation failed");

    println!("=== DC Operating Point ===");
    for plot in &result.plots {
        for vec in &plot.vecs {
            if let Some(&val) = vec.real.first() {
                println!("  {:>20} = {:>12.6} ", vec.name, val);
            }
        }
    }
}
