{
  description = "groebner";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay.url = "github:oxalica/rust-overlay";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      nixpkgs,
      rust-overlay,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ (import rust-overlay) ];
        };
        # Pinned by rust-toolchain.toml, which rustup users and CI read too.
        rustToolchain = pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
      in
      {
        devShells.default = pkgs.mkShell {
          packages = [
            (rustToolchain.override {
              extensions = [
                "rust-src"
                "rust-analyzer"
              ];
            })
            # The Python bindings (`python/`): setuptools-rust builds the
            # cdylib with the `ffi` feature, ruff and mypy check the wrapper.
            pkgs.uv
            pkgs.ruff
            (pkgs.python3.withPackages (ps: [
              ps.setuptools
              ps.setuptools-rust
              ps.mypy
            ]))
          ];
        };
      }
    );
}
