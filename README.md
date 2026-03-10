# thevenin

A SPICE circuit simulator written in Rust, targeting both native and WebAssembly.

Thevenin is a from-scratch rewrite of [ngspice](https://ngspice.sourceforge.io/) in idiomatic Rust. It parses standard SPICE netlists and runs simulations with the same algorithms — Modified Nodal Analysis, Newton-Raphson iteration, sparse direct solvers — but with Rust's type safety and `wasm32` support.

## Crates

| Crate | Description |
|-------|-------------|
| [`thevenin`](thevenin/) | Simulation engine — MNA assembly, Newton-Raphson solver, analysis drivers |
| [`thevenin-types`](thevenin-types/) | Netlist parser and IR types for the ngspice SPICE dialect |

## Supported analyses

- **DC operating point** (`.op`)
- **DC sweep** (`.dc`)
- **AC small-signal** (`.ac`)
- **Transient** (`.tran`)
- **Noise** (`.noise`)
- **Sensitivity** (`.sens`)
- **Transfer function** (`.tf`)
- **Pole-zero** (`.pz`)

## Device models

- Resistors, capacitors, inductors
- Independent voltage/current sources (DC, AC, pulse, sin, PWL, exp, AM, SFFM)
- Dependent sources (VCVS, VCCS, CCVS, CCCS)
- Diodes
- BJTs (Gummel-Poon)
- MOSFETs (level 1-3)
- JFETs
- BSIM3v3, BSIM4
- BSIM3SOI (FD, PD, DD)
- VBIC

## Quick start

```rust
use thevenin::simulate_op;
use thevenin_types::Netlist;

let netlist = Netlist::parse("
Voltage Divider
V1 in 0 1.0
R1 in mid 1k
R2 mid 0 2k
.op
.end
").unwrap();

let result = simulate_op(&netlist).unwrap();
for plot in &result.plots {
    for vec in &plot.vecs {
        if let Some(&val) = vec.real.first() {
            println!("{}: {:.6}", vec.name, val);
        }
    }
}
```

## Building

```bash
# native
cargo build
cargo test --workspace

# wasm32
cargo test --workspace --target wasm32-unknown-unknown
```

## License

BSD-3-Clause
