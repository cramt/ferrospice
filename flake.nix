{
  description = "Ferrospice - ngspice rewrite in Rust";

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
          extensions = ["rust-src" "rust-analyzer" "clippy" "rustfmt"];
        };

        craneLib = (crane.mkLib pkgs).overrideToolchain rustToolchain;

        src = craneLib.cleanCargoSource ./.;

        commonCraneArgs = {
          inherit src;
          strictDeps = true;
        };

        cargoArtifacts = craneLib.buildDepsOnly commonCraneArgs;

        ferrospice = craneLib.buildPackage (commonCraneArgs // {
          inherit cargoArtifacts;
          doCheck = true;
        });

        ci-build = pkgs.writeShellScriptBin "ci-build" ''
          set -euo pipefail
          echo "=== Building ferrospice ==="
          ${pkgs.nix}/bin/nix build .#ferrospice --print-build-logs
          echo ""
          echo "=== Build complete ==="
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
          default = ferrospice;
          inherit ferrospice ci-build update-deps;
        };

        apps.default = flake-utils.lib.mkApp {
          drv = ferrospice;
        };

        apps.ci-build = flake-utils.lib.mkApp {
          drv = ci-build;
        };

        apps.update-deps = flake-utils.lib.mkApp {
          drv = update-deps;
        };

        devShells.default = pkgs.mkShell {
          packages = with pkgs; [
            rustToolchain
            pkg-config
            openssl
            git
            jq
          ];

          shellHook = ''
            export RUST_BACKTRACE=1
            export RUST_LOG=info
            echo "=== Ferrospice dev environment ==="
            echo "  rustc: $(rustc --version)"
            echo ""
            echo "  Build:   cargo build"
            echo "  Test:    cargo test"
            echo "  Check:   cargo clippy --workspace -- -D warnings"
          '';
        };
      }
    );
}
