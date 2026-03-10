#!/usr/bin/env bash
# bench-targets.sh — Run internal benchmarks on native x86 and wasm32, compare results.
#
# Usage:
#   ./scripts/bench-targets.sh           # Run both targets, show comparison
#   ./scripts/bench-targets.sh --json    # Output JSON instead of table
#
# Requires: cargo, wasm-bindgen-test-runner, chromium (all provided by nix develop)
set -euo pipefail

JSON_MODE=false
[[ "${1:-}" == "--json" ]] && JSON_MODE=true

BOLD='\033[1m'
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
RESET='\033[0m'

run_bench() {
    local label="$1"
    shift
    # Run and capture stderr (where BENCH lines go)
    cargo test -p thevenin --test bench "$@" -- --nocapture 2>&1 | grep '^BENCH ' || true
}

parse_results() {
    # Input: lines like "BENCH <name> <iters> <total_ns> <per_iter_ns>"
    # Output: associative array via declare -p
    local -n out="$1"
    while IFS=' ' read -r _ name iters total_ns per_iter_ns; do
        out["$name"]="$per_iter_ns"
    done
}

fmt_ns() {
    local ns="$1"
    if (( ns >= 1000000000 )); then
        awk "BEGIN { printf \"%.2fs\", $ns / 1000000000 }"
    elif (( ns >= 1000000 )); then
        awk "BEGIN { printf \"%.2fms\", $ns / 1000000 }"
    elif (( ns >= 1000 )); then
        awk "BEGIN { printf \"%.2fμs\", $ns / 1000 }"
    else
        printf "%dns" "$ns"
    fi
}

echo -e "${BOLD}=== Thevenin Benchmark: native x86 vs wasm32 ===${RESET}"
echo ""

# --- Pre-build both targets so compilation doesn't count ---
echo -e "${CYAN}Building native...${RESET}"
cargo test -p thevenin --test bench --no-run 2>&1 | tail -1
echo -e "${CYAN}Building wasm32...${RESET}"
cargo test -p thevenin --test bench --target wasm32-unknown-unknown --no-run 2>&1 | tail -1
echo ""

# --- Run benchmarks ---
echo -e "${CYAN}Running native benchmarks...${RESET}"
native_raw=$(run_bench "native")

echo -e "${CYAN}Running wasm32 benchmarks...${RESET}"
wasm_raw=$(run_bench "wasm32" --target wasm32-unknown-unknown)

echo ""

# --- Parse results ---
declare -A native_results
declare -A wasm_results

while IFS=' ' read -r _ name iters total_ns per_iter_ns; do
    native_results["$name"]="$per_iter_ns"
done <<< "$native_raw"

while IFS=' ' read -r _ name iters total_ns per_iter_ns; do
    wasm_results["$name"]="$per_iter_ns"
done <<< "$wasm_raw"

# --- Collect all benchmark names ---
declare -A all_names
for k in "${!native_results[@]}" "${!wasm_results[@]}"; do
    all_names["$k"]=1
done

sorted_names=($(printf '%s\n' "${!all_names[@]}" | sort))

if $JSON_MODE; then
    # --- JSON output ---
    echo "["
    first=true
    for name in "${sorted_names[@]}"; do
        native_ns="${native_results[$name]:-null}"
        wasm_ns="${wasm_results[$name]:-null}"
        ratio="null"
        if [[ "$native_ns" != "null" && "$wasm_ns" != "null" && "$native_ns" != "0" ]]; then
            ratio=$(awk "BEGIN { printf \"%.2f\", $wasm_ns / $native_ns }")
        fi
        $first || echo ","
        first=false
        printf '  {"name":"%s","native_ns":%s,"wasm_ns":%s,"ratio":%s}' \
            "$name" "$native_ns" "$wasm_ns" "$ratio"
    done
    echo ""
    echo "]"
    exit 0
fi

# --- Table output ---
printf "${BOLD}%-30s %14s %14s %12s${RESET}\n" "Benchmark" "Native" "Wasm32" "Ratio"
printf "%-30s %14s %14s %12s\n" \
    "$(printf '%.0s─' {1..30})" \
    "$(printf '%.0s─' {1..14})" \
    "$(printf '%.0s─' {1..14})" \
    "$(printf '%.0s─' {1..12})"

for name in "${sorted_names[@]}"; do
    native_ns="${native_results[$name]:-}"
    wasm_ns="${wasm_results[$name]:-}"

    native_fmt="—"
    wasm_fmt="—"
    ratio_fmt="—"

    [[ -n "$native_ns" ]] && native_fmt=$(fmt_ns "$native_ns")
    [[ -n "$wasm_ns" ]] && wasm_fmt=$(fmt_ns "$wasm_ns")

    if [[ -n "$native_ns" && -n "$wasm_ns" && "$native_ns" != "0" ]]; then
        ratio=$(awk "BEGIN { printf \"%.2f\", $wasm_ns / $native_ns }")
        if awk "BEGIN { exit !($wasm_ns > $native_ns) }"; then
            ratio_fmt=$(printf "${RED}%sx${RESET}" "$ratio")
        else
            ratio_fmt=$(printf "${GREEN}%sx${RESET}" "$ratio")
        fi
    fi

    printf "%-30s %14s %14s %12b\n" "$name" "$native_fmt" "$wasm_fmt" "$ratio_fmt"
done

echo ""
echo -e "${BOLD}Ratio = wasm/native (>1 means wasm is slower)${RESET}"
