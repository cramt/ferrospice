{
  description = "Thevenin - ngspice rewrite in Rust";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    crane = {
      url = "github:ipetkov/crane";
    };
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
    rust-overlay,
    crane,
  }:
    flake-utils.lib.eachDefaultSystem (
      system: let
        overlays = [(import rust-overlay)];
        pkgs = import nixpkgs {inherit system overlays;};

        rustToolchain = pkgs.rust-bin.stable.latest.default.override {
          targets = ["wasm32-unknown-unknown" "wasm32-unknown-emscripten" "wasm32-wasip1"];
          extensions = ["rust-src" "rust-analyzer" "clippy" "rustfmt"];
        };

        craneLib = (crane.mkLib pkgs).overrideToolchain rustToolchain;

        src = craneLib.cleanCargoSource ./.;

        commonCraneArgs = {
          inherit src;
          strictDeps = true;
        };

        cargoArtifacts = craneLib.buildDepsOnly commonCraneArgs;

        thevenin = craneLib.buildPackage (commonCraneArgs
          // {
            inherit cargoArtifacts;
            doCheck = true;
          });

        ci-build = pkgs.writeShellScriptBin "ci-build" ''
          set -euo pipefail
          echo "=== Building thevenin ==="
          ${pkgs.nix}/bin/nix build .#thevenin --print-build-logs
          echo ""
          echo "=== Build complete ==="
        '';

        test-wasm = pkgs.writeShellScriptBin "test-wasm" ''
          set -euo pipefail
          echo "=== Running tests in WebAssembly (Headless Chrome) ==="
          ${rustToolchain}/bin/cargo test --target wasm32-unknown-unknown
        '';

        update-deps = pkgs.writeShellScriptBin "update-deps" ''
          set -euo pipefail
          echo "=== Updating all dependencies ==="

          echo ""
          echo "--- Nix flake inputs ---"
          ${pkgs.nix}/bin/nix flake update

          echo ""
          echo "--- Cargo dependencies ---"
          ${rustToolchain}/bin/cargo update

          echo ""
          echo "=== All dependencies updated ==="
        '';
      in {
        packages = {
          default = thevenin;
          inherit thevenin ci-build update-deps test-wasm;
        };

        apps.default = flake-utils.lib.mkApp {
          drv = thevenin;
        };

        apps.ci-build = flake-utils.lib.mkApp {
          drv = ci-build;
        };

        apps.update-deps = flake-utils.lib.mkApp {
          drv = update-deps;
        };

        apps.test-wasm = flake-utils.lib.mkApp {
          drv = test-wasm;
        };

        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            rustToolchain
            pkg-config
            openssl
            git
            jq
            wasmtime
            wasm-bindgen-cli
            nodejs
            chromium
            chromedriver
            cargo-nextest
          ];

          shellHook = ''
            export RUST_BACKTRACE=1
            export RUST_LOG=info
            echo "=== Thevenin dev environment ==="
            echo "  rustc: $(rustc --version)"
            echo ""
            echo "  Build:   cargo build"
            echo "  Test:    cargo test"
            echo "  Wasm:    cargo test --target wasm32-unknown-unknown"
            echo "  Check:   cargo clippy --workspace -- -D warnings"
          '';
        };
      }
    );
}
