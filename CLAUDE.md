# Thevenin

A Rust rewrite of [ngspice](https://ngspice.sourceforge.io/), the open-source mixed-level/mixed-signal electronic circuit simulator.

## Project Structure

- `ngspice-upstream/` - Cloned ngspice C source (reference implementation, do not modify)
- `src/` - Rust source code
- `scripts/ralph/` - Ralph autonomous agent loop for incremental implementation

## Development

**Always run commands through `nix develop --command ...`** so that flake.nix stays honest
and no dependency works by accident from the host environment.

```bash
nix develop --command cargo build
nix develop --command cargo test
nix develop --command cargo clippy --workspace -- -D warnings
nix develop --command cargo fmt --check
```

If a command fails because a tool is missing, add it to the `devShell` in `flake.nix`
rather than installing it on the host.

## Ralph Loop

Run the ralph loop for autonomous implementation:
```bash
./scripts/ralph/ralph.sh
```

### Test fixtures
`ngspice-upstream/` is gitignored, so any test fixture files sourced from it (e.g. `.cir`,
`.spice`, model files) **must be copied** into a tracked fixtures directory (e.g.
`tests/fixtures/` or `<crate>/tests/fixtures/`) before they are referenced by tests. This
ensures tests work on a clean clone without needing the upstream repo checked out.

## Project Goals

### Test-driven development
Write tests first, implement second. Port test cases from ngspice wherever possible — the
`ngspice-upstream/tests/` directory and example circuits in `ngspice-upstream/examples/` are
gold. If ngspice has a test for a behavior, steal it. The ngspice source is the authoritative
reference for correctness.

### Correctness over performance
Get it right first. Optimize later, guided by profiling, not intuition. A slow correct
simulator is infinitely more useful than a fast wrong one. Do not reach for `unsafe`,
hand-rolled SIMD, or clever bit tricks until profiling proves they're needed.

### Make impossible states irrepresentable
Encode invariants in the type system. Use enums over stringly-typed fields, newtypes over
bare primitives, and builder patterns or typestate where construction has constraints. If a
function can't fail for a given input type, the type should make that obvious. Prefer
compile-time guarantees over runtime checks.

### Prefer facet + unsynn over serde + syn
Use `facet` for reflection/serialization and `unsynn` for any proc-macro or parsing work.
Do not pull in `serde` or `syn` unless there is a hard dependency that requires them.

### Split into subcrates freely
Don't let the workspace stay monolithic. Extract crates along natural boundaries — parser,
IR, device models, solvers, output formatting, etc. Small focused crates compile faster,
test in isolation, and enforce API boundaries. Use a Cargo workspace from the start.

## Conventions

- Follow the original ngspice architecture where it makes sense, but use idiomatic Rust
- Use `thiserror` for error types
- Use `nalgebra` or `ndarray` for matrix operations (SPICE is heavily linear algebra)
- Prefer safe Rust; use `unsafe` only when absolutely necessary for performance-critical paths
- The ngspice source in `ngspice-upstream/` is the authoritative reference for behavior
- Never include Co-Authored-By lines in git commits
