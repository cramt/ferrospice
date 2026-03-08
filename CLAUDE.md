# Ferrospice

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
