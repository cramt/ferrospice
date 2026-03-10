//! AC frequency response of a low-pass RC filter.
//!
//! Circuit: V1 (AC=1) -> R1 (1kΩ) -> node "out" -> C1 (1µF) -> GND
//! Cutoff frequency: f_c = 1/(2π·R·C) ≈ 159 Hz

use thevenin::simulate_ac;
use thevenin_types::Netlist;

fn main() {
    let netlist = Netlist::parse(
        "\
RC Low-Pass Filter
V1 in 0 DC 0 AC 1
R1 in out 1k
C1 out 0 1u
.ac dec 10 1 100k
.end
",
    )
    .expect("failed to parse netlist");

    let result = simulate_ac(&netlist).expect("simulation failed");

    let plot = &result.plots[0];
    let freq = &plot.vecs[0]; // frequency
    let vout = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(out)")
        .expect("v(out) not found");

    println!("=== RC Low-Pass Filter Frequency Response ===");
    println!(
        "{:>12}  {:>12}  {:>10}",
        "Freq (Hz)", "|V(out)|", "Phase (°)"
    );
    println!(
        "{:>12}  {:>12}  {:>10}",
        "---------", "--------", "---------"
    );
    for (i, f) in freq.real.iter().enumerate() {
        let re = vout.complex[i].re;
        let im = vout.complex[i].im;
        let mag = (re * re + im * im).sqrt();
        let phase_deg = im.atan2(re).to_degrees();
        println!("{f:>12.1}  {mag:>12.6}  {phase_deg:>10.2}");
    }
}
