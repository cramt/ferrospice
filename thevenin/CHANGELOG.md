# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0](https://github.com/cramt/thevenin/compare/thevenin-v0.1.0...thevenin-v0.2.0) - 2026-03-11

### Added

- US-057 - PZ complex pair formatting and sensitivity table pagination
- US-055 - Interpolation-aware transient output comparison
- US-054 - Fix DC sweep variable column and add numeric tolerance comparison
- US-065 - XSPICE A-element parser with bracketed port groups
- US-039 - Test harness matching ngspice check.sh
- US-038 - Regression tests for models, subcircuit processing, and misc
- US-037 - Regression tests for parser, func, and lib-processing
- US-036 - TXL and CPL transmission line models
- US-035 - Lossy transmission line model (LTRA)
- US-034 HFET and MESFET (Statz/Curtice) models
- US-032 MESA FET model (Ytterdal/Lee/Shur/Fjeldly GaAs MESFET)
- US-033 MOS6 MOSFET model (Sakurai-Newton n-th power)
- VBIC AC/noise analysis, charge model, and 10 passing tests (US-031)

### Other

- ignore failing harness tests with issue references
- extract shared physics constants and safe_exp into physics.rs
- rename project from ferrospice to thevenin across all files
- extract shared NR device stamping into device_stamp module
- release v0.1.0

## [0.1.0](https://github.com/cramt/thevenin/releases/tag/thevenin-v0.1.0) - 2026-03-10

### Fixed

- add version to thevenin-types dependency for crates.io publishing

### Other

- format resistance and transient tests
- copy ngspice test fixtures into repo so CI doesn't need ngspice-upstream
- add release-plz and CI workflows, add crates.io metadata
- rename crates: core crate is now `thevenin`, parser crate is now `thevenin-types`
