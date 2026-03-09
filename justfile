# Ferrospice development commands

# Run benchmark comparison: native x86 vs wasm32
bench:
    nix develop --command bash scripts/bench-targets.sh

# Run benchmark comparison, output JSON
bench-json:
    nix develop --command bash scripts/bench-targets.sh --json
