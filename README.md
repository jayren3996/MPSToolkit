# MPSToolkit.jl

`MPSToolkit.jl` is a finite-`MPS` toolkit built on top of [`ITensors.jl`](https://itensor.github.io/ITensors.jl/stable/) and `ITensorMPS.jl`. It is organized around explicit tensor-network building blocks rather than one opaque driver, so the low-level evolution, truncation, observable, operator-space, and spectral pieces all stay directly callable.

[![Documentation - Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jayren3996.github.io/MPSToolkit/dev/)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jayren3996/MPSToolkit.git")
```

## Functionality

- dense-gate TEBD and helper-driven schedule construction from local Hamiltonians
- MPO-based TDVP for finite OBC `MPS`
- ScarFinder workflows with explicit projection, energy targeting, and selector-driven refinement
- Pauli-basis operator-space tools for coherent and open-system evolution
- operator-space DMT and DAOE / FDAOE projectors
- Chebyshev moments, energy-window truncation, and spectral reconstruction
- entanglement, energy, and fidelity diagnostics

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

## Documentation And Examples

- [Development documentation](https://jayren3996.github.io/MPSToolkit/dev/)
- [Examples index](https://jayren3996.github.io/MPSToolkit/dev/examples/)
- [TEBD scheduler notebook](examples/tebd/scheduler_patterns.ipynb)
- [DMT scheduler notebook](examples/operator_space/dmt_scheduler.ipynb)
- [Chebyshev energy-cutoff notebook](examples/chebyshev/energy_cutoff_comparison.ipynb)

The hosted docs contain the fuller manual, API reference, and workflow-specific examples. GitHub is best treated as the landing page, while the Documenter site is the main reference.
