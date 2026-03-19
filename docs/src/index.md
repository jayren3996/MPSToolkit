# MPSToolkit.jl

`MPSToolkit.jl` is a finite-`MPS` toolkit built on top of `ITensors.jl` and `ITensorMPS.jl`. It is designed for users who want explicit control over tensor-network workflows instead of one monolithic driver.

The package keeps the important moving parts visible:

- TEBD is represented as dense local gates plus an explicit schedule
- TDVP is represented as an MPO generator plus explicit solver settings
- ScarFinder is a short loop that calls public evolution, projection, energy matching, and selector routines
- operator-space workflows keep the Pauli-basis representation explicit
- Chebyshev tools expose both moment generation and spectral reconstruction

## What MPSToolkit Is Good For

Use `MPSToolkit.jl` when you want to:

- build finite-chain TEBD workflows from dense local Hamiltonians
- run MPO-based TDVP on finite open-boundary `MPS`
- study explicit projection-and-refinement loops such as ScarFinder
- work in Pauli operator space for open-system or projector-based calculations
- compute Chebyshev moments and reconstruct spectral functions

## Design Philosophy

The package is intentionally modular. Higher-level workflows are thin wrappers around reusable lower-level functions, so the same pieces can be recombined for custom projects, experiments, and debugging.

That means the manual is organized by subsystem rather than by one giant "run everything" interface. If you already know what kind of computation you want, it is usually best to jump straight to the corresponding manual page.

## Current Limits

- finite OBC `MPS` is the main supported state class
- periodic chains are not a general supported mode
- `TDVPEvolution` currently expects MPO generators
- DMT is currently implemented for operator-space workflows

## Navigation

- [Getting Started](getting-started.md)
  Minimal setup, installation, and first examples.
- [Architecture](manual/architecture.md)
  How the codebase is organized and how the main building blocks relate.
- [TEBD And TDVP](manual/tebd-tdvp.md)
  Dense-gate TEBD helpers and MPO-based TDVP.
- [ScarFinder](manual/scarfinder.md)
  Explicit projection workflows, selector refinement, and ScarFinder-specific step-count guidance.
- [Operator Space](manual/operator-space.md)
  Pauli-basis states, gates, and operator-space evolution tools.
- [DAOE](manual/daoe.md)
  DAOE and FDAOE projector construction.
- [DMT](manual/dmt.md)
  Operator-space DMT truncation and scheduling.
- [Chebyshev](manual/chebyshev.md)
  Moment generation, energy-window cutoff, and spectral reconstruction.
- [Examples](examples.md)
  Script and notebook entry points.
- [API Reference](api.md)
  Source-backed API docs grouped by subsystem.
