# MPSToolkit Example Redesign

**Goal:** Replace the current ad hoc example set with a curated notebook-first tutorial suite that introduces the entire public API in a structured way.

## Scope

The redesign replaces all existing files under `examples/` and `examples/notebooks/` with a smaller, coherent set:

- one overview notebook
- focused notebooks for TEBD, TDVP, observables/bases, and ScarFinder
- matching lightweight Julia scripts for direct terminal execution

The examples remain finite-`MPS` only and reflect the current supported package surface.

## Example Structure

New scripts:

- `examples/tebd_spin1_xy.jl`
- `examples/tdvp_spin1_xy.jl`
- `examples/observables_and_bases.jl`
- `examples/scarfinder_pxp_tebd.jl`
- `examples/scarfinder_pxp_tdvp.jl`

New notebooks:

- `examples/notebooks/00_overview.ipynb`
- `examples/notebooks/01_evolution_tebd.ipynb`
- `examples/notebooks/02_evolution_tdvp.ipynb`
- `examples/notebooks/03_observables_and_bases.ipynb`
- `examples/notebooks/04_scarfinder_tebd.ipynb`
- `examples/notebooks/05_scarfinder_tdvp.ipynb`

## Teaching Strategy

The notebooks should progress from low-level to high-level:

1. package map and namespace overview
2. TEBD evolution details
3. TDVP evolution details
4. observables and Pauli-basis helpers
5. TEBD-based ScarFinder workflow
6. TDVP-based ScarFinder workflow

Each notebook should explain arguments and keyword controls in plain language, not just run code.

## Model Choices

Use a small stable model set:

- spin-1 XY for evolution tutorials
- PXP for ScarFinder tutorials

This avoids model churn while still covering the important public API.

## Documentation Changes

`README.md` and `docs/examples.md` should become navigational indices for the new example set. Historical filenames and references to the removed examples should disappear so the docs match the repository exactly.
