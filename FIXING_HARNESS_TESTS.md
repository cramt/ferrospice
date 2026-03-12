# Fixing Ignored Harness Tests

Systematic methodology for diagnosing and fixing `#[ignore]`d tests in
`thevenin/tests/harness.rs`.  Each test compares thevenin's batch output
against the reference `.out` file from `ngspice-upstream/tests/`.

---

## 1. Pick a test

Choose one test (or a group with the same ignore reason).  Prefer tests whose
ignore reason sounds like a numerical bug over tests that require whole missing
features.

**Triage by category (most → least tractable):**

| Category | Count | Approach |
|---|---|---|
| Numerical offset / accuracy | few | Compare equations term-by-term against C |
| tran: singular matrix | 2 | Check device stamp completeness, initial conditions |
| tran: initial values wrong | 1 | Compare DC OP → first tran step transition |
| tran: output mismatch | 1 | Diff output at each timepoint, find where it diverges |
| transient timestep (US-055) | 12 | Timestep control / output interpolation |
| device info output (US-061) | 10 | Missing `.print` columns (showmod, element params) |
| AC complex formatting (US-058) | 2 | Output formatter doesn't emit complex columns |
| .plot ASCII art (US-058) | 1 | Missing `.plot` renderer |
| BSIM3SOI accuracy (US-059) | 15 | Model equation bugs — big effort |
| sensitivity param mismatch | 1 | Missing BJT sensitivity params |
| parameter expressions | 1 | Resistor array `{expr}` not evaluated |
| PZ numerical accuracy | 1 | Eigenvalue solver bug for inductors |
| MESA non-default temp | 1 | Temperature scaling bug in MESA model |
| BSIM1/BSIM2 (US-052/053) | 2 | Entire model not implemented |
| XSPICE (US-056) | 3 | Entire subsystem not implemented |
| .control scripting | 35 | Entire interpreter not implemented |
| TEMPER keyword | 4 | Needs expression-in-parameter + .control |

## 2. Run the test and capture the diff

```bash
nix develop --command cargo test --package thevenin \
  --test harness TESTNAME -- --ignored --nocapture 2>&1 | head -80
```

The harness prints both the expected (filtered) and actual (filtered) output,
plus the first mismatch.  Key things to note:

- **Which column fails** (col 0 = first output variable, etc.)
- **Constant offset vs relative error vs divergence** — constant offset
  suggests a parasitic term (gmin, leakage); relative error suggests a formula
  bug; divergence suggests convergence failure or wrong region selection.
- **At what sweep/time point it starts failing** — if only the first point is
  wrong, it's likely an initial-condition or default-value issue.

## 3. Understand the circuit

Read the `.cir` file from `thevenin/tests/fixtures/<subdir>/` (or
`ngspice-upstream/tests/<subdir>/`).  Identify:

- What analysis: `.dc`, `.tran`, `.ac`, `.pz`, `.noise`, `.sens`?
- What devices: which model type and level?
- What is being measured: `.print` variables?
- What `.model` parameters are set (non-default)?

## 4. Locate the relevant code

The ngspice C source is authoritative.  Key directories:

```
ngspice-upstream/src/spicelib/devices/<model>/     # device models
ngspice-upstream/src/spicelib/analysis/             # analysis drivers
ngspice-upstream/src/maths/ni/                      # NR iteration (niiter.c)
ngspice-upstream/src/maths/sparse/                  # sparse matrix (spsmp.c)
```

Device model files follow a naming convention:
- `<model>load.c`  — NR load function (stamps Jacobian + RHS)
- `<model>defs.h`  — instance/model structs, parameter IDs
- `<model>temp.c`  — temperature-dependent preprocessing
- `<model>dset.c`  — parameter setup / defaults
- `<model>acld.c`  — AC small-signal load

Our Rust equivalents live in `thevenin/src/<model>.rs`.

### Rust codebase architecture: companion → stamp → device_stamp

Every nonlinear device follows the same pattern in our codebase:

1. **`<model>.rs`** — contains two key functions:
   - `<model>_companion(inst, vgs, vgd, gmin) → <Model>Companion` — pure
     computation: given terminal voltages, returns linearised conductances
     (`gm`, `gds`, `ggs`, `ggd`), currents (`cd`, `cg`), and capacitances.
     This is the Rust equivalent of the computation section in
     `<model>load.c`.
   - `stamp_<model>_with_voltages(comp, inst, vgs, vgd, matrix, rhs)` —
     applies the companion model to the MNA matrix and RHS vector (Norton
     equivalent current sources + Y-matrix conductance stamps).  This is the
     equivalent of the `load:` label section in `<model>load.c`.

2. **`device_stamp.rs`** — the glue layer.  `DeviceVoltageState::stamp_devices()`
   is called at every NR iteration.  For each device type it:
   - Extracts terminal voltages from the solution vector
   - Applies voltage limiting (`pnjlim`, `fetlim`)
   - Calls the `_companion()` function
   - Calls the `stamp_…()` function
   This is the shared NR load logic used by both `simulate.rs` (DC OP / DC
   sweep) and `transient.rs`.

3. **`simulate.rs`** — DC analysis driver.  Creates the `load` closure that
   copies base linear stamps, adds LTRA/TXL/CPL DC equations, and calls
   `stamp_devices()`.  Passes this closure to `newton_raphson_solve()`.

4. **`transient.rs`** — transient analysis driver.  Same pattern but also
   stamps reactive elements (capacitor/inductor companion models) and handles
   timestep control.

When debugging a device model bug, you almost always only need to look at
`<model>.rs`.  When debugging a stamping or solver bug, look at
`device_stamp.rs` and `newton.rs`.

### The output pipeline: format → filter → compare

The harness test pipeline (`harness.rs`) works as follows:

1. **`format_batch_output(netlist, result)`** (`output.rs:16`) — takes the
   `SimResult` and produces text output mimicking ngspice `--batch` mode.
   Emits title, temperature line, data tables with headers.  This is where
   "device info output (US-061)" failures originate — ngspice prints extra
   sections like `.OP` element parameters, `.SHOWMOD` model info, and small-
   signal parameters that our formatter doesn't emit.

2. **`filter_output(text)`** (`output.rs:547`) — applies the same filtering as
   ngspice's `check.sh` FILTER regex.  Removes lines containing keywords like
   `"Circuit"`, `"Index"`, `"Date"`, `"---"`, etc.  Both expected and actual
   output are filtered before comparison.

3. **`compare_filtered(expected, actual)`** (`output.rs:654`) — compares the
   two filtered outputs.  Non-numeric tokens are compared exactly.  Numeric
   tokens use relative tolerance 1e-4 and absolute tolerance 1e-15.  If line
   counts differ (e.g. different timestep counts), falls back to
   interpolation-aware comparison that linearly interpolates actual data at
   the expected data's x-coordinates.

**Common output-related failure modes:**
- **"device info output (US-061)"** — the `.out` file contains device
  parameter sections (element info, operating point details) that we don't
  emit.  After filtering, these sections still leave residual lines that don't
  appear in our output.  Fix: add the relevant output sections to
  `format_batch_output`, or extend `filter_output` if the lines should be
  stripped.
- **"AC complex output formatting (US-058)"** — AC analysis outputs complex
  numbers (magnitude + phase, or real + imaginary) that our formatter doesn't
  handle.  The `SimVector.complex` field exists but
  `format_batch_output` doesn't emit it in the right format.
- **"transient timestep (US-055)"** — our transient produces different
  timepoints than ngspice.  The interpolation fallback in `compare_filtered`
  tries to handle this, but large timestep differences cause interpolation
  errors.  Fix: improve timestep control in `transient.rs` to match ngspice's
  adaptive stepping algorithm.

## 5. Systematic diff against ngspice C

This is the core debugging technique.  For a numerical bug:

### 5a. Isolate the failing quantity

From the test output, identify the variable (e.g. `vids#branch`) and the
operating point where the mismatch occurs.  Build a minimal mental model of
what currents/voltages should flow.

### 5b. Compare device equations line-by-line

Open the ngspice `<model>load.c` and our `<model>.rs` side by side.  Check:

1. **Junction / diode currents** — reverse bias approximation, forward bias
   exponential, `csat` scaling by area.
2. **Drain/collector current** — region selection (cutoff / linear / saturation),
   formula for `cdrain`, `gm`, `gds`.
3. **Conductance signs and stamps** — the Y-matrix stamp pattern.  Compare
   every `*(here->...Ptr) += m * (...)` line in C against the corresponding
   `matrix.add(row, col, value)` in Rust.
4. **Current accounting** — `cd = cdrain - cgd`, `cg = cgs + cgd`, Norton
   equivalents `ceqgd`, `ceqgs`, `cdreq`.
5. **Default parameter values** — check `<model>dset.c` or the `init` function.
   A wrong default (especially for `b`, `lambda`, `alpha`, `is`) shifts results.

### 5c. Check the solver / infrastructure

Some bugs are not in the device model but in the framework:

- **Gmin handling** — device models include `CKTgmin` in junction conductances.
  The solver diagonal `CKTdiagGmin` is separate (see `spsmp.c:LoadGmin`).
  In ngspice, `CKTdiagGmin` starts at 0 for DC analysis and is only elevated
  during Gmin stepping.  Our `NrOptions.diag_gmin` mirrors this — DC sweep
  sets it to 0, other paths keep it at `gmin`.
- **Voltage limiting** — `DEVpnjlim`, `DEVfetlim` in ngspice vs our `pnjlim`,
  `fetlim`.  Wrong limiting causes NR oscillation.
- **Temperature** — many models scale parameters at `TEMP != TNOM`.  Check
  `<model>temp.c`.  Our code may skip temperature adjustment entirely.
- **Charge model / capacitances** — for transient, the charge integration
  (trapezoidal or Gear) must match.  Wrong capacitance → wrong transient.
- **Output formatting** — the harness compares filtered text output.  If
  column ordering, number formatting, or header text differs, the comparison
  fails even if the numbers are correct.

## 6. Write a focused unit test

Before fixing, write a small test that reproduces the bug at the device-model
level (not the full harness).  For example:

```rust
#[test]
fn mesfet_subthreshold_junction_current() {
    let inst = make_test_instance();
    let comp = mesfet_companion(&inst, -3.0, -3.1, 1e-12);
    // Expected: cd ≈ 3.114e-12 (from ngspice reference)
    assert_abs_diff_eq!(comp.cd, 3.114e-12, epsilon = 1e-14);
}
```

This makes the fix testable independently of the full simulation pipeline.

## 7. Fix and verify

Apply the fix.  Then:

```bash
# 1. Verify the target test passes
nix develop --command cargo test --package thevenin \
  --test harness TESTNAME -- --ignored --nocapture

# 2. Run the FULL test suite — check for regressions
nix develop --command cargo test --package thevenin 2>&1 | grep FAILED

# 3. Clippy
nix develop --command cargo clippy --workspace -- -D warnings
```

**Regression risk areas:**
- Changing `newton.rs` (NR solver) affects every nonlinear circuit.
- Changing a device model affects all tests using that model.
- Changing output formatting affects every harness test.

## 8. Un-ignore the test

In `harness.rs`, change:
```rust
harness_test!(name, "path/file.cir", ignore = "reason");
```
to:
```rust
harness_test!(name, "path/file.cir");
```

## 9. Check if sibling tests also pass

Many ignore reasons are shared across a group.  After fixing the root cause,
try the whole group:

```bash
nix develop --command cargo test --package thevenin \
  --test harness -- --ignored --nocapture 2>&1 | grep -E "ok|FAILED"
```

---

## Worked Example: MES subthreshold offset

**Symptom:** `harness_mes_subth` — constant +2e-13 offset in drain current
across entire subthreshold DC sweep.

**Diagnosis:**
1. Circuit: MESFET DC sweep, Vgs from -3V to 0V, measuring `vids#branch`.
2. In deep subthreshold, drain current ≈ 0; measured current is dominated by
   reverse-biased gate-drain junction leakage ≈ `gmin × V(drain)`.
3. Expected ≈ 3.114e-12, got ≈ 3.314e-12.  Diff = 2e-13 = `gmin × 0.2V`
   (roughly the voltage at drain_prime).
4. Root cause: our NR solver added `options.gmin` to every matrix diagonal
   (newton.rs line 108), but ngspice's `CKTdiagGmin` is 0 for DC analysis.
   Device models already include gmin internally in junction conductances, so
   the diagonal addition was double-counting.

**Fix:** Added `NrOptions.diag_gmin` field.  DC sweep sets `diag_gmin = 0`;
other paths keep it at `gmin` for backward compatibility.  The direct NR
attempt in `newton_raphson_solve` uses `diag_gmin` instead of `gmin` for the
diagonal addition.

**Verification:** MES test passes.  All 465+ existing tests still pass.
No LTRA/transient regressions.

---

## Important prior fix: `diag_gmin` (affects all numerical comparisons)

Before investigating any numerical offset in a DC sweep, be aware of this
already-applied fix.  The `NrOptions` struct has two separate gmin fields:

- **`gmin`** (default 1e-12) — used by device models in junction conductance
  equations.  Always active.  This is ngspice's `CKTgmin`.
- **`diag_gmin`** (default 1e-12) — added from every node to ground by the NR
  solver during matrix factorization.  This is ngspice's `CKTdiagGmin`.

In ngspice, `CKTdiagGmin` starts at **0** and is only elevated during Gmin
stepping.  For DC sweep (`simulate_dc` in `simulate.rs`), we set
`diag_gmin = 0` to match.  Other code paths (transient, standalone OP) keep
the default `diag_gmin = gmin` for backward compatibility.

If you see a constant ~1e-13 to ~1e-12 offset in a DC sweep test, check
whether the code path is using `diag_gmin = 0`.  If a new analysis driver is
added, it should set `diag_gmin` appropriately based on whether ngspice would
use `CKTdiagGmin = 0` for that context.
