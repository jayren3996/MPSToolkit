# MPSToolkit.jl

`MPSToolkit.jl` is a finite-`MPS` toolkit built on top of [`ITensors.jl`](https://itensor.github.io/ITensors.jl/stable/) and `ITensorMPS.jl`. The package is organized around explicit tensor-network building blocks rather than one opaque driver, so evolution, projection, observables, operator-space tools, and spectral routines stay directly callable.

[![Documentation - Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jayren3996.github.io/MPSToolkit/dev/)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jayren3996/MPSToolkit.git")
```

## Quick Start

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 6
sites = siteinds("S=1/2", nsites)
psi = MPS(sites, n -> isodd(n) ? "Up" : "Dn")

evolution = tebd_strang_evolution(
  nsites,
  0.05;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.8),
  maxdim=32,
  cutoff=1e-12,
)

evolve!(psi, evolution)
println(expect(psi, "Sz"))
```

## For AI Agents

This README is intentionally written as an agent-oriented landing page. If you are making changes in the codebase, these are the main ideas to keep in working memory.

### Project Mental Model

- `MPSToolkit` prefers explicit workflows over hidden orchestration.
- Evolution is represented by small config objects such as `LocalGateEvolution`, `DMTGateEvolution`, and `TDVPEvolution`.
- High-level workflows like ScarFinder are deliberately thin wrappers around public building blocks like `evolve!`, `project!`, `match_energy!`, and selector scoring.
- Most user-facing routines mutate `MPS` objects in place and return the same object for convenience.

### Source Map

- `src/evolution/`
  TEBD and TDVP configuration types plus the concrete `evolve!` implementations.
- `src/scarfinder/`
  Projection settings, selector types, explicit ScarFinder loop, and post-step energy matching.
- `src/observables/`
  Energy density, entanglement entropy, entanglement spectrum, and fidelity-style diagnostics.
- `src/operator_space/`
  Pauli-basis helpers, DMT, and DAOE / FDAOE projectors.
- `src/chebyshev/`
  Chebyshev moments, energy-window projection, damping kernels, and spectral reconstruction.
- `src/models/` and `src/bases/`
  Small dense helper matrices used by examples and helper constructors.

### Important Conventions

- `evolve!`, `project!`, and `energy_density` are dispatch points. Check the concrete methods before changing behavior.
- Dense local operators are interpreted through the local site dimension of the input `MPS`.
- Operator-space code assumes the local Pauli ordering `(I, X, Y, Z)` unless a docstring says otherwise.
- Finite OBC `MPS` is the main supported state class. Periodic behavior is only supported in narrow helper cases.
- Most top-level workflows are intentionally composable. If a change would hide a low-level primitive behind a large driver, it is probably moving away from the intended design.

### ScarFinder-Specific Caveat

- ScarFinder now treats an effective evolution step count of `1` as a bad main-loop setting.
- If `scarfinder_step!`, `scarfinder!`, or `trajectory_refine!` receive `LocalGateEvolution`, `DMTGateEvolution`, or `TDVPEvolution` with effective steps `1`, they emit a warning and internally use `10` for that ScarFinder call.
- This rule is local to ScarFinder. Global TEBD, DMT, and TDVP constructors still keep their original defaults.
- Internal energy-correction substeps inside `match_energy!` still use single-step updates on purpose.

### Good Starting Files For Common Tasks

- Add or change TEBD behavior:
  `src/evolution/types.jl`, `src/evolution/tebd.jl`, `docs/src/manual/tebd-tdvp.md`
- Change ScarFinder behavior:
  `src/scarfinder/types.jl`, `src/scarfinder/dispatch.jl`, `src/scarfinder/selectors.jl`, `src/scarfinder/algorithm.jl`, `docs/src/manual/scarfinder.md`
- Change operator-space helpers:
  `src/operator_space/helpers.jl`, `src/operator_space/dmt.jl`, `src/operator_space/daoe.jl`
- Change spectral tooling:
  `src/chebyshev/types.jl`, `src/chebyshev/moments.jl`, `src/chebyshev/reconstruction.jl`

## Main Functionality

- dense-gate TEBD and helper-driven schedule construction from local Hamiltonians
- MPO-based TDVP for finite OBC `MPS`
- ScarFinder workflows with explicit projection, energy targeting, selector-driven refinement, and ScarFinder-specific step-count validation
- Pauli-basis operator-space tools for coherent and open-system evolution
- operator-space DMT and DAOE / FDAOE projectors
- Chebyshev moments, energy-window truncation, and spectral reconstruction
- entanglement, energy, and fidelity diagnostics

## Current Limits

- finite OBC `MPS` is the main supported state class
- periodic chains are not a general supported mode
- `TDVPEvolution` currently expects MPO generators
- DMT is currently implemented for operator-space workflows

## Documentation And Examples

- [Development documentation](https://jayren3996.github.io/MPSToolkit/dev/)
- [Getting started guide](https://jayren3996.github.io/MPSToolkit/dev/getting-started/)
- [ScarFinder manual](https://jayren3996.github.io/MPSToolkit/dev/manual/scarfinder/)
- [API reference](https://jayren3996.github.io/MPSToolkit/dev/api/)
- [Examples index](https://jayren3996.github.io/MPSToolkit/dev/examples/)
- [TEBD scheduler notebook](examples/tebd/scheduler_patterns.ipynb)
- [DMT scheduler notebook](examples/operator_space/dmt_scheduler.ipynb)
- [Chebyshev energy-cutoff notebook](examples/chebyshev/energy_cutoff_comparison.ipynb)

The hosted Documenter site is the main human-facing reference. This README is intentionally more operational and codebase-oriented so that AI agents and contributors can orient themselves quickly before opening the manual pages.
