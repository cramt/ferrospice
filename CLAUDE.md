# Ferrospice

A Rust rewrite of [ngspice](https://ngspice.sourceforge.io/), the open-source mixed-level/mixed-signal electronic circuit simulator.

## Project Structure

- `ngspice-upstream/` - Cloned ngspice C source (reference implementation, do not modify)
- `src/` - Rust source code
- `scripts/ralph/` - Ralph autonomous agent loop for incremental implementation

## Development

- Build: `cargo build`
- Test: `cargo test`
- Check: `cargo clippy -- -D warnings`
- Format: `cargo fmt --check`

## Ralph Loop

Run the ralph loop for autonomous implementation:
```bash
./scripts/ralph/ralph.sh --tool claude
```

## Conventions

- Follow the original ngspice architecture where it makes sense, but use idiomatic Rust
- Use `thiserror` for error types
- Use `nalgebra` or `ndarray` for matrix operations (SPICE is heavily linear algebra)
- Prefer safe Rust; use `unsafe` only when absolutely necessary for performance-critical paths
- Write tests alongside implementation, referencing ngspice's expected behavior
- The ngspice source in `ngspice-upstream/` is the authoritative reference for behavior
- Never include Co-Authored-By lines in git commits
