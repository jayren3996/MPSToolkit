# MPSToolkit.jl

`MPSToolkit.jl` is a finite-`MPS` toolkit built on top of `ITensors.jl` and `ITensorMPS.jl`.

The package is organized around explicit tensor-network building blocks rather than one opaque driver. You can use the higher-level workflows when they help, but the lower-level evolution, truncation, observable, operator-space, and Chebyshev pieces remain directly callable.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jayren3996/MPSToolkit.git")
```

## What It Covers

- dense-gate TEBD on finite OBC chains
- helper-driven TEBD setup from local Hamiltonians
- MPO-based TDVP
- ScarFinder projection and selector workflows
- operator-space DMT with scheduled gate evolution
- Pauli-basis helpers for operator-space calculations
- DAOE / FDAOE projector MPOs
- Chebyshev moments and spectral reconstruction
- entanglement, energy, and fidelity diagnostics

## Current Limits

- finite OBC `MPS` is the main supported state class
- periodic chains are not a general supported mode
- `TDVPEvolution` currently expects MPO generators
- DMT is currently implemented for operator-space workflows

## Navigation

- [Getting Started](getting-started.md)
- [Architecture](manual/architecture.md)
- [TEBD And TDVP](manual/tebd-tdvp.md)
- [ScarFinder](manual/scarfinder.md)
- [Operator Space](manual/operator-space.md)
- [DAOE](manual/daoe.md)
- [DMT](manual/dmt.md)
- [Chebyshev](manual/chebyshev.md)
- [Examples](examples.md)
- [API Reference](api.md)
