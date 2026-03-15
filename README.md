# MPSToolkit.jl

`MPSToolkit.jl` is a low-level finite-`MPS` toolkit built on top of [`ITensors.jl`](https://itensor.github.io/ITensors.jl/stable/) and `ITensorMPS.jl`.

[![Documentation - Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jayren3996.github.io/MPSToolkit/stable/)
[![Documentation - Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jayren3996.github.io/MPSToolkit/dev/)

The package is organized around explicit tensor-network building blocks instead of a single opaque workflow. You can use the high-level loops when they help, but the lower-level evolution, truncation, observable, and operator-space pieces are all exposed directly.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jayren3996/MPSToolkit.git")
```

## What It Covers

Current public functionality includes:

- finite OBC `MPS` evolution from dense local gates
- helper-driven TEBD setup from dense local Hamiltonians
- odd-even-odd Strang TEBD helpers
- scheduled multi-site gate application, including custom per-bond and mixed-span schedules
- MPO-based TDVP
- explicit bond-dimension projection for ScarFinder workflows
- selector-driven trajectory refinement with entropy and fidelity selectors
- entanglement and energy diagnostics
- spin-1/2 Pauli-basis helpers for operator-space calculations
- operator-space Lindblad and local superoperator helpers
- operator-space DMT with low-level `dmt_step!` and scheduled `DMTGateEvolution`
- DAOE/FDAOE projector MPO construction
- Chebyshev moment generation and spectral reconstruction helpers

## Current Limits

- finite OBC `MPS` is the main supported state class
- finite-ring PBC is not a general supported mode
- `TDVPEvolution` currently expects MPO generators
- DMT is currently implemented for operator-space workflows, not as a generic physical-state truncation backend

## Package Layout

Feature namespaces:

- `MPSToolkit.Evolution`
- `MPSToolkit.Observables`
- `MPSToolkit.ScarFinder`
- `MPSToolkit.Bases`
- `MPSToolkit.OperatorSpace`
- `MPSToolkit.Models`
- `MPSToolkit.Chebyshev`

The design is intentionally explicit:

1. `evolve!(psi, evolution)`
2. `project!(psi, truncation)`
3. optional `match_energy!`
4. optional trajectory refinement

That split keeps TEBD, TDVP, DMT, and ScarFinder-style workflows composable rather than fused into one control path.

## Evolution APIs

TEBD and scheduled local-gate evolution:

- `LocalGateEvolution`
- `tebd_evolve!`
- `local_gates_from_hamiltonians`
- `tebd_evolution_from_hamiltonians`
- `tebd_strang_schedule`
- `tebd_strang_evolution`

Operator-space DMT:

- `DMTGateEvolution`
- `dmt_step!`
- `dmt_evolve!`

TDVP:

- `TDVPEvolution`
- `tdvp_evolve!`

ScarFinder and explicit projection:

- `BondDimTruncation`
- `EnergyTarget`
- `EntropySelector`
- `FidelitySelector`
- `scarfinder_step!`
- `scarfinder!`

Observables:

- `energy_density`
- `bond_entropy`
- `entanglement_spectrum`
- `fidelity_distance`

Operator-space helpers:

- `pauli_siteinds`
- `pauli_basis_state`
- `pauli_total_sz_state`
- `pauli_gate`
- `pauli_gate_from_hamiltonian`
- `pauli_lindblad_generator`
- `pauli_gate_from_lindbladian`
- `pauli_daoe_projector`
- `fdaoe_projector`

Model helpers:

- `spinhalf_matrices`
- `spinhalf_tfim_bond_hamiltonian`
- `spinhalf_xyz_bond_hamiltonian`

Chebyshev helpers:

- `ChebyshevRescaling`
- `SpectralFunction`
- `chebyshev_moments`
- `jackson_damping`
- `jackson_kernel`
- `reconstruct_chebyshev`
- `spectral_function`

## Quick Start

### Physical-Spin TEBD From Hamiltonians

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

for _ in 1:4
  evolve!(psi, evolution)
end

println(expect(psi, "Sz"))
println(maxlinkdim(psi))
```

### Operator-Space TEBD

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 6
sites = pauli_siteinds(nsites)
state = pauli_basis_state(sites, ["Z", "I", "I", "I", "I", "I"])

evolution = tebd_strang_evolution(
  nsites,
  0.05;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.8),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=64,
  cutoff=1e-12,
)

evolve!(state, evolution)
println(inner(pauli_basis_state(sites, ["Z", "I", "I", "I", "I", "I"]), state))
```

### Scheduled Operator-Space DMT

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 6
sites = pauli_siteinds(nsites)
state = pauli_basis_state(sites, ["Z", "I", "I", "I", "I", "I"])

gate = pauli_gate_from_hamiltonian(spinhalf_tfim_bond_hamiltonian(nsites, 1; J=1.0, g=0.6), 0.05)
schedule = collect(1:(nsites - 1))

evolution = DMTGateEvolution(
  fill(gate, length(schedule)),
  0.05;
  schedule=schedule,
  reverse_schedule=collect(reverse(schedule)),
  maxdim=16,
  cutoff=1e-10,
)

dmt_evolve!(state, evolution)
println(maxlinkdim(state))
```

## Example Notebooks And Scripts

Examples are grouped by workflow:

- `examples/tebd/`
- `examples/operator_space/`
- `examples/open_systems/`
- `examples/chebyshev/`
- `examples/scarfinder/`
- `examples/daoe/`
- `examples/tdvp/`

Good starting points:

- [examples/tebd/tebd_helper_apis.ipynb](examples/tebd/tebd_helper_apis.ipynb)
- [examples/tebd/scheduler_patterns.ipynb](examples/tebd/scheduler_patterns.ipynb)
- [examples/operator_space/operator_tebd_helper_apis.ipynb](examples/operator_space/operator_tebd_helper_apis.ipynb)
- [examples/operator_space/dmt_scheduler.ipynb](examples/operator_space/dmt_scheduler.ipynb)
- [examples/open_systems/pauli_lindblad_tebd.ipynb](examples/open_systems/pauli_lindblad_tebd.ipynb)
- [examples/tebd/xxz_tebd_vs_ed.ipynb](examples/tebd/xxz_tebd_vs_ed.ipynb)
- [examples/chebyshev/tfim_local_spectrum.ipynb](examples/chebyshev/tfim_local_spectrum.ipynb)

## Documentation

Additional docs:

- [Stable documentation](https://jayren3996.github.io/MPSToolkit/stable/)
- [Development documentation](https://jayren3996.github.io/MPSToolkit/dev/)
- [docs/index.md](docs/index.md)
- [docs/architecture.md](docs/architecture.md)
- [docs/api.md](docs/api.md)
- [docs/examples.md](docs/examples.md)
