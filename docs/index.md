# MPSToolkit Documentation

`MPSToolkit.jl` is a low-level finite-`MPS` package for evolution, observables, and ScarFinder-style search workflows built on top of `ITensors.jl` and `ITensorMPS.jl`.

The package is organized around one shared algorithm and finite OBC backend-specific evolution/truncation methods:

- dense-gate TEBD
- MPO-based TDVP
- explicit bond-dimension projection
- selector-based refinement on projected trajectories
- Pauli-basis helpers for spin-1/2 operator work

Start here:

- [Architecture](/Users/ren/Codex/MPSToolkit/docs/architecture.md)
- [API](/Users/ren/Codex/MPSToolkit/docs/api.md)
- [Examples](/Users/ren/Codex/MPSToolkit/docs/examples.md)

## Current Limits

- only finite OBC `MPS` is supported
- finite-ring PBC is not implemented
- TEBD applies dense local gates on explicit finite schedules
- `TDVPEvolution` currently accepts MPO generators only
- DAOE/operator-space evolution is not implemented yet
