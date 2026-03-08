# PRD: Ferrospice Simulation Engine

## Introduction

Build a SPICE circuit simulation engine in Rust as a library crate (`ferrospice`). This is a ground-up rewrite of ngspice's simulation core, using the existing `ferrospice-netlist` crate for parsing. The engine must produce identical numerical results to ngspice for all supported analyses.

**The project is done when ferrospice passes the entire ngspice base test suite** (113 `.cir` tests in `ngspice-upstream/tests/`). The test suite includes resistance circuits, filters, transient analysis, device models (BJT, MOSFET, JFET, MESFET, transmission lines), pole-zero analysis, sensitivity analysis, and regression tests.

This is a **library-first** project. Ferrospice is shipped as a Rust crate for embedding in other tools. A CLI may be added for convenience and test-running, but the public API surface is the library.

## Goals

- Implement Modified Nodal Analysis (MNA) matrix assembly and solving
- Support all core SPICE analyses: `.op`, `.dc`, `.tran`, `.ac`, `.noise`, `.tf`, `.sens`, `.pz`
- Implement device models: R, C, L, V, I, D, Q (BJT), M (MOSFET), J (JFET), transmission lines
- Pass all 113 ngspice base test suite `.cir` files with output matching `.out` reference files
- Expose a clean, idiomatic Rust library API
- Test-driven: port ngspice test cases as Rust integration tests before implementing features

## User Stories

### Phase 1: Foundation (MNA + DC Operating Point)

### US-001: Create ferrospice-core crate with sparse matrix types
**Description:** As a developer, I need a sparse matrix representation suitable for MNA so that circuit equations can be assembled and solved.

**Acceptance Criteria:**
- [ ] New crate `ferrospice-core` in the workspace
- [ ] Sparse matrix struct supporting dynamic insertion and triplet-form assembly
- [ ] LU factorization and solve (use `faer` or similar — pick what fits, pure Rust preferred)
- [ ] Unit tests: assemble a 3x3 system, solve it, verify result within 1e-12 tolerance
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-002: MNA matrix assembly from netlist
**Description:** As a developer, I need to convert a parsed `Netlist` into an MNA equation system so that the circuit can be solved.

**Acceptance Criteria:**
- [ ] Function that takes a `ferrospice_netlist::Netlist` and produces an MNA system (matrix + RHS vector)
- [ ] Node name → matrix index mapping (ground node "0" is the reference, not in matrix)
- [ ] Stamps for: resistor, independent voltage source (adds branch equation), independent current source
- [ ] Unit test: voltage divider (V1=5V, R1=1k, R2=1k) → MNA system has correct dimensions and entries
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-003: DC operating point solver (.op)
**Description:** As a developer, I need a DC operating point solver so that node voltages and branch currents can be computed for a DC circuit.

**Acceptance Criteria:**
- [ ] `simulate_op(netlist: &Netlist) -> SimResult` function in the public API
- [ ] Returns node voltages and voltage source branch currents
- [ ] Test: voltage divider V1=5V, R1=1k, R2=1k → V(mid) = 2.5V within 1e-9
- [ ] Test: series resistors with current source → correct node voltages
- [ ] Port `ngspice-upstream/tests/resistance/res_simple.cir` as a Rust integration test
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-004: Capacitor and inductor device stamps
**Description:** As a developer, I need C and L device stamps so that reactive elements are represented in the MNA system (for DC: C=open, L=short).

**Acceptance Criteria:**
- [ ] Capacitor stamp: open circuit in DC (no contribution to DC matrix)
- [ ] Inductor stamp: short circuit in DC (zero-resistance branch equation)
- [ ] Test: RC circuit DC op → capacitor acts as open, voltage across C equals source voltage through divider
- [ ] Test: RL circuit DC op → inductor acts as short, full current flows
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-005: Port ngspice resistance test suite
**Description:** As a developer, I need the ngspice resistance tests running against ferrospice to validate basic DC solving.

**Acceptance Criteria:**
- [ ] Integration test module that loads `.cir` files from `ngspice-upstream/tests/resistance/`
- [ ] Parse each `.cir` with `ferrospice_netlist`, simulate with ferrospice, compare output to `.out`
- [ ] Output comparison uses the same filtering as ngspice's `check.sh` (ignore metadata, compare numerical values)
- [ ] All 3 resistance tests pass (res_simple, res_tc, res_temp — or skip temp-dependent ones with `#[ignore]` and a TODO)
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 2: DC Sweep + Diode (Nonlinear)

### US-006: DC sweep analysis (.dc)
**Description:** As a developer, I need DC sweep analysis so that a source can be swept across a range and the operating point computed at each step.

**Acceptance Criteria:**
- [ ] `.dc src start stop step` sweeps the named source and computes OP at each point
- [ ] Results returned as `SimResult` with one data point per sweep step
- [ ] Test: sweep V1 from 0 to 5V in 1V steps across a resistor → 6 data points with correct I(V1)
- [ ] Double sweep (`.dc V1 ... V2 ...`) supported
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-007: Newton-Raphson nonlinear solver
**Description:** As a developer, I need a Newton-Raphson iteration loop so that nonlinear devices (diodes, transistors) can be solved.

**Acceptance Criteria:**
- [ ] NR iteration with configurable convergence tolerance (ABSTOL, RELTOL, VNTOL matching ngspice defaults)
- [ ] Iteration limit (default 100, matching ngspice's ITL1)
- [ ] Convergence check: both absolute and relative criteria on node voltages and branch currents
- [ ] Source stepping or Gmin stepping fallback when NR fails to converge
- [ ] Test: solve a simple nonlinear equation (e.g., diode I-V) via NR, verify convergence in < 20 iterations
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-008: Diode device model
**Description:** As a developer, I need a diode device model so that nonlinear circuits with PN junctions can be simulated.

**Acceptance Criteria:**
- [ ] Shockley diode equation: I = IS * (exp(V/Vt) - 1)
- [ ] `.model` parameters: IS, N, RS, BV (breakdown), CJO (junction capacitance for later)
- [ ] NR stamp: companion model (conductance + current source) linearized at each iteration
- [ ] Test: diode with 1V source and 1k resistor → forward voltage ~0.6-0.7V
- [ ] Test: diode I-V curve via `.dc` sweep matches ngspice within 1e-6 relative tolerance
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 3: Transient Analysis

### US-009: Transient analysis engine (.tran)
**Description:** As a developer, I need transient (time-domain) analysis so that circuits with time-varying sources and reactive elements can be simulated.

**Acceptance Criteria:**
- [ ] Backward Euler integration method (simplest, matches ngspice's default starter)
- [ ] Trapezoidal integration method (ngspice's default after startup)
- [ ] Capacitor companion model: stamp updates each timestep based on integration method
- [ ] Inductor companion model: stamp updates each timestep based on integration method
- [ ] Fixed timestep initially (tstep from `.tran` command)
- [ ] Test: RC step response — charge curve matches analytical V(t) = V0*(1-exp(-t/RC)) within 1% at 5 time constants
- [ ] Test: LC oscillator — frequency matches 1/(2*pi*sqrt(LC)) within 1%
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-010: Time-varying source waveforms
**Description:** As a developer, I need transient waveform evaluation (PULSE, SIN, PWL, etc.) so that time-dependent sources drive the circuit during `.tran`.

**Acceptance Criteria:**
- [ ] PULSE waveform: v1, v2, td, tr, tf, pw, per — all parameters
- [ ] SIN waveform: v0, va, freq, td, theta, phi
- [ ] PWL waveform: piecewise linear interpolation between time-value pairs
- [ ] EXP, SFFM, AM waveforms
- [ ] Test: PULSE source into RC circuit → output matches ngspice `res_simple.cir` transient output
- [ ] Test: SIN source → output samples match analytical sin() within 1e-6
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-011: Adaptive timestep control
**Description:** As a developer, I need adaptive timestep control so that transient analysis is accurate and efficient.

**Acceptance Criteria:**
- [ ] Local truncation error (LTE) estimation using the difference between BE and trap results
- [ ] Timestep adjustment: reduce if LTE > RELTOL, increase if LTE << RELTOL
- [ ] Respect tmax from `.tran` if specified
- [ ] Breakpoint handling: force timestep at source discontinuities (PULSE edges, PWL corners)
- [ ] Test: PULSE into RC → adaptive stepping places steps near edges, coasts during flat regions
- [ ] Port `ngspice-upstream/tests/transient/fourbitadder.cir` as integration test (can be `#[ignore]` until BJT model exists)
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 4: AC Small-Signal Analysis

### US-012: AC analysis (.ac)
**Description:** As a developer, I need AC small-signal analysis so that frequency response can be computed.

**Acceptance Criteria:**
- [ ] Linearize circuit at DC operating point
- [ ] Build complex-valued MNA matrix (G + jwC)
- [ ] Sweep frequency: DEC, OCT, LIN modes
- [ ] AC source magnitude and phase applied correctly
- [ ] Results: complex node voltages at each frequency point
- [ ] Test: RC lowpass filter, -3dB point at f = 1/(2*pi*RC) within 1%
- [ ] Port `ngspice-upstream/tests/filters/` test as integration test
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 5: Semiconductor Device Models

### US-013: BJT (bipolar junction transistor) model
**Description:** As a developer, I need a BJT device model so that bipolar transistor circuits can be simulated.

**Acceptance Criteria:**
- [ ] Ebers-Moll model (level 1, matching ngspice's default BJT)
- [ ] NPN and PNP support
- [ ] `.model` parameters: BF, BR, IS, VAF, VAR, IKF, IKR, ISE, ISC, NE, NC, RB, RC, RE, CJE, CJC, CJS, TF, TR
- [ ] NR stamps for collector and base currents
- [ ] Test: common-emitter amplifier DC bias point matches ngspice within 1e-6
- [ ] Port `ngspice-upstream/tests/general/` tests that use BJTs
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-014: MOSFET level 1 model
**Description:** As a developer, I need a basic MOSFET model so that MOS circuits can be simulated.

**Acceptance Criteria:**
- [ ] Shichman-Hodges model (level 1)
- [ ] NMOS and PMOS support
- [ ] `.model` parameters: VTO, KP, LAMBDA, PHI, GAMMA, CBD, CBS, CGSO, CGDO, CGBO
- [ ] Operating regions: cutoff, linear, saturation with correct equations
- [ ] NR stamps for drain current
- [ ] Test: NMOS inverter transfer curve matches ngspice within 1e-6
- [ ] Port `ngspice-upstream/tests/mos6/` tests (or simplest MOSFET tests available)
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-015: JFET model
**Description:** As a developer, I need a JFET device model.

**Acceptance Criteria:**
- [ ] Standard JFET model matching ngspice's level 1
- [ ] N-channel and P-channel support
- [ ] `.model` parameters: VTO, BETA, LAMBDA, RD, RS, CGS, CGD, IS
- [ ] Test: JFET amplifier bias matches ngspice
- [ ] Port `ngspice-upstream/tests/jfet/` test
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 6: Subcircuits + Controlled Sources

### US-016: Subcircuit expansion (.subckt / X)
**Description:** As a developer, I need subcircuit instantiation so that hierarchical circuits can be flattened and simulated.

**Acceptance Criteria:**
- [ ] Flatten nested `.subckt` definitions into a flat netlist before MNA assembly
- [ ] Port mapping: subcircuit ports connected to caller's nets
- [ ] `PARAMS:` parameter substitution in subcircuit instances
- [ ] Nested subcircuits (subcircuit containing another subcircuit call)
- [ ] Test: voltage divider as subcircuit, instantiated and simulated, matches direct circuit
- [ ] Port `ngspice-upstream/tests/regression/subckt-processing/` tests
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-017: Controlled sources (E, F, G, H)
**Description:** As a developer, I need linear controlled sources so that dependent sources can be simulated.

**Acceptance Criteria:**
- [ ] VCVS (E): voltage-controlled voltage source
- [ ] CCCS (F): current-controlled current source
- [ ] VCCS (G): voltage-controlled current source
- [ ] CCVS (H): current-controlled voltage source
- [ ] MNA stamps for each type
- [ ] Test: op-amp modeled as high-gain VCVS, inverting amplifier gain = -R2/R1
- [ ] Port `ngspice-upstream/tests/general/` tests that use controlled sources
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 7: Advanced Analyses

### US-018: Noise analysis (.noise)
**Description:** As a developer, I need noise analysis to compute spectral noise density.

**Acceptance Criteria:**
- [ ] Thermal noise for resistors: 4kTR
- [ ] Shot noise for diodes/BJTs: 2qI
- [ ] Flicker noise (1/f) for supported devices
- [ ] Output: noise spectral density at each frequency point
- [ ] Test: resistor thermal noise matches 4kTR analytically
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-019: Transfer function (.tf) and sensitivity (.sens) analyses
**Description:** As a developer, I need .tf and .sens analyses.

**Acceptance Criteria:**
- [ ] `.tf output input` computes small-signal transfer function, input/output resistance
- [ ] `.sens output` computes sensitivity of output to each component
- [ ] Port `ngspice-upstream/tests/sensitivity/` test
- [ ] Port `ngspice-upstream/tests/regression/sens/` tests
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-020: Pole-zero analysis (.pz)
**Description:** As a developer, I need pole-zero analysis.

**Acceptance Criteria:**
- [ ] Compute poles and zeros of the transfer function
- [ ] Port `ngspice-upstream/tests/polezero/` tests (6 tests)
- [ ] Port `ngspice-upstream/tests/regression/pz/` test
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 8: Advanced Device Models

### US-021: BSIM3 MOSFET model
**Description:** As a developer, I need the BSIM3 model for accurate short-channel MOSFET simulation.

**Acceptance Criteria:**
- [ ] BSIM3v3 model implementation (or port from ngspice C code)
- [ ] NMOS and PMOS
- [ ] Port `ngspice-upstream/tests/bsim3/` tests
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-022: BSIM4 MOSFET model
**Description:** As a developer, I need the BSIM4 model.

**Acceptance Criteria:**
- [ ] BSIM4 model implementation
- [ ] Port `ngspice-upstream/tests/bsim4/` tests
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-023: BSIM-SOI models
**Description:** As a developer, I need BSIM3-SOI device models (PD, FD, DD variants).

**Acceptance Criteria:**
- [ ] BSIM3-SOI Partially Depleted, Fully Depleted, and Double Diffused models
- [ ] Port `ngspice-upstream/tests/bsim3soipd/`, `bsim3soifd/`, `bsim3soidd/` tests (18 tests total)
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-024: VBIC bipolar model
**Description:** As a developer, I need the VBIC bipolar model.

**Acceptance Criteria:**
- [ ] VBIC model implementation
- [ ] Port `ngspice-upstream/tests/vbic/` tests (6 tests)
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-025: MESA, HFET, MOS6 and remaining device models
**Description:** As a developer, I need the remaining device models to complete test suite coverage.

**Acceptance Criteria:**
- [ ] MESA FET model — port `ngspice-upstream/tests/mesa/` (11 tests)
- [ ] HFET model — port `ngspice-upstream/tests/hfet/` (2 tests)
- [ ] MOS6 model — port `ngspice-upstream/tests/mos6/` (2 tests)
- [ ] MESFET model — port `ngspice-upstream/tests/mes/` (1 test)
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-026: Transmission line models
**Description:** As a developer, I need transmission line models (LTRA, TXL, CPL).

**Acceptance Criteria:**
- [ ] Lossy transmission line (LTRA) model
- [ ] Coupled multiconductor line (TXL/CPL) models
- [ ] Port `ngspice-upstream/tests/transmission/` tests (6 tests)
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### Phase 9: Regression Suite + CLI

### US-027: Port full regression test suite
**Description:** As a developer, I need all regression tests passing to ensure parser, control flow, and edge cases work.

**Acceptance Criteria:**
- [ ] Port `ngspice-upstream/tests/regression/misc/` (13 tests)
- [ ] Port `ngspice-upstream/tests/regression/lib-processing/` (7 tests)
- [ ] Port `ngspice-upstream/tests/regression/parser/` (5 tests)
- [ ] Port `ngspice-upstream/tests/regression/temper/` (4 tests)
- [ ] Port `ngspice-upstream/tests/regression/model/` (3 tests)
- [ ] Port `ngspice-upstream/tests/regression/func/` (1 test)
- [ ] All 41 regression tests pass
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-028: Test harness matching ngspice's check.sh
**Description:** As a developer, I need an automated test harness that runs all 113 `.cir` files and compares output to `.out` references, matching ngspice's filtering rules.

**Acceptance Criteria:**
- [ ] Rust integration test or binary that iterates `ngspice-upstream/tests/**/*.cir`
- [ ] For each `.cir`: parse, simulate, format output, filter, diff against `.out`
- [ ] Filtering matches `check.sh` rules (strip metadata, normalize floats)
- [ ] Summary report: N/113 tests passing
- [ ] All 113 tests pass
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

### US-029: Minimal CLI for running netlists
**Description:** As a user, I want a CLI to run SPICE netlists so I can use ferrospice from the command line.

**Acceptance Criteria:**
- [ ] `ferrospice --batch input.cir` runs all analyses and prints results to stdout
- [ ] Output format matches ngspice's batch output (for test compatibility)
- [ ] CLI is a thin wrapper around the library — all logic lives in library crates
- [ ] `cargo clippy --workspace -- -D warnings` passes
- [ ] `cargo test --workspace` passes

## Functional Requirements

- FR-1: Library crate `ferrospice` exposes `simulate(netlist: &Netlist) -> Result<SimResult>` as the primary API
- FR-2: All analysis types (`.op`, `.dc`, `.tran`, `.ac`, `.noise`, `.tf`, `.sens`, `.pz`) are supported
- FR-3: All linear devices (R, C, L, V, I, E, F, G, H, K) have correct MNA stamps
- FR-4: All nonlinear devices (D, Q, M, J) use Newton-Raphson with companion models
- FR-5: Subcircuit expansion handles nesting, parameter passing, and port mapping
- FR-6: Transient analysis supports BE and trapezoidal integration with adaptive timestep
- FR-7: AC analysis builds complex MNA matrix and sweeps frequency
- FR-8: Output format for CLI matches ngspice batch mode (for test suite compatibility)
- FR-9: Numerical results match ngspice within tolerances used by `check.sh` diff comparison
- FR-10: `.model` and `.param` directives are fully supported

## Non-Goals (Out of Scope)

- XSpice mixed-signal / code models (digital gates, state machines, etc.)
- Interactive mode / control language (`.control` / `.endc` blocks)
- GUI or plotting
- OSDI device model loading
- Shared library / C FFI interface
- BSIM-CMG (FinFET) or other bleeding-edge models not in the base test suite
- Monte Carlo / statistical analysis
- Parallel / multi-threaded simulation

## Technical Considerations

- **Workspace structure:** Split into subcrates along natural boundaries:
  - `ferrospice-netlist` — parser (already exists)
  - `ferrospice-core` — sparse matrix, MNA assembly, solver infrastructure
  - `ferrospice-devices` — device model implementations (stamps, companion models)
  - `ferrospice` — top-level library tying it together + CLI binary
- **Linear algebra:** Use `faer` or similar pure-Rust sparse solver. Avoid C dependencies.
- **Device model trait:** Define a `Device` trait with methods for DC stamp, AC stamp, transient stamp, noise contribution. Each device type implements this trait.
- **Test infrastructure:** Build a test harness that can load `.cir`/`.out` pairs from ngspice-upstream and run them as Rust integration tests.
- **Reference code:** ngspice C source in `ngspice-upstream/` is the authoritative reference for equations, convergence algorithms, and expected behavior.

## Success Metrics

- **All 113 ngspice base test suite `.cir` files produce output matching their `.out` reference files** using the same filtering and comparison rules as ngspice's `check.sh`
- Library API is ergonomic: simulate a circuit in < 10 lines of Rust
- No `unsafe` code outside of clearly justified performance-critical sections
- Zero `serde` or `syn` dependencies — use `facet` and `unsynn`

## Open Questions

- Should we support `.control` blocks for the regression tests that use self-validating control scripts, or rewrite those as native Rust assertions?
- Which BSIM3/4 variant parameters are actually exercised by the test suite? (May be able to implement a subset)
- Should temperature-dependent parameters (`.temp`, TC1/TC2) be in Phase 1 or deferred?
- Mutual inductance (K element) — needed for any base tests? If not, defer.
