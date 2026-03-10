# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0](https://github.com/cramt/thevenin/releases/tag/ferrospice-v0.1.0) - 2026-03-10

### Added

- US-045 - Examples directory with runnable examples
- US-030 - VBIC bipolar model (level 4)
- US-029 - BSIM3SOI-DD MOSFET model (level 56)
- US-027 - BSIM3SOI-PD MOSFET model (level 57)
- add cross-target benchmarking suite (native x86 vs wasm32)
- US-026 - BSIM4 MOSFET model — AC, charge model, noise, and test suite
- US-024 - BSIM4 MOSFET model — parameter parsing and temperature preprocessing
- US-023 - BSIM3 AC, charge model, and test suite
- US-002 - MNA matrix assembly from netlist
- US-001 - Create ferrospice-core crate with sparse matrix types

### Fixed

- fix tests

### Other

- add release-plz and CI workflows, add crates.io metadata
- rename ferrospice-core to thevenin and ferrospice-netlist to thevenin-types
- a
- update progress log and PRD for US-045 examples implementation
- issue tags
- update progress log and PRD for US-028 BSIM3SOI-FD implementation
- update progress log and PRD for architecture pass (iteration 5)
- update PRD and progress log for US-025 completion
- tell ralph agent to always push after committing
- add multi-pass system to ralph loop
- Merge branch 'main' into ralph/simulation-engine
- cleanup
- better loop
- update PRD and progress for US-020
- update PRD and progress for US-019
- a
- update PRD and progress for US-018
- update PRD and progress for US-017
- update PRD and progress for US-016
- update PRD and progress for US-015
- update PRD and progress for US-014
- update PRD and progress for US-013
- update PRD and progress for US-012
- update PRD and progress for US-011
- update PRD and progress for US-010
- update PRD and progress for US-009
- update PRD and progress for US-008
- update PRD and progress for US-007
- update PRD and progress for US-006
- update PRD and progress for US-005
- update PRD and progress for US-004
- update PRD and progress for US-003
- Use portable shebang in ralph.sh
- Rename spice-netlist to ferrospice-netlist, add PRD and ralph config
- Add flake.nix with crane, rust-overlay, and dev shell
- Initial project setup for ferrospice - ngspice rewrite in Rust
