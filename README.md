# MPSToolkit.jl

Low-level finite-`MPS` tooling built on `ITensors.jl` and `ITensorMPS.jl`.

`MPSToolkit.jl` is a finite-chain package organized around explicit tensor-network building blocks. The current supported workflows are:

- finite OBC `MPS` time evolution from local gates or MPO Hamiltonians
- finite OBC `MPS` ScarFinder loops with explicit projection
- dense local-gate TEBD, including multi-site and per-bond gate schedules
- MPO-driven TDVP
- entanglement entropy and entanglement spectrum diagnostics
- entropy-based and fidelity-based trajectory refinement
- spin-1/2 Pauli-basis helpers for operator-space work
- DAOE/FDAOE projector MPO construction for operator-space truncation

The package does not currently support:

- finite-ring PBC `MPS`
- automatic `OpSum` conversion inside `TDVPEvolution`
- DAOE/FDAOE-specific operator evolution workflows

## Installation

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Design

The package keeps the algorithm intentionally explicit:

1. `evolve!(psi, evolution)`
2. `project!(psi, truncation)`
3. optional `match_energy!`
4. optional trajectory refinement with a selector such as `EntropySelector` or `FidelitySelector`

That split makes the backend tunable. TEBD and TDVP are both first-class evolution engines, while truncation stays visible and replaceable.

## Public API

Feature namespaces:

- `MPSToolkit.Evolution`
- `MPSToolkit.Observables`
- `MPSToolkit.ScarFinder`
- `MPSToolkit.Bases`
- `MPSToolkit.OperatorSpace`
- `MPSToolkit.Models`

Core entry points:

- `scarfinder_step!(psi, evolution, truncation; target_energy=nothing, selector=nothing, kwargs...)`
- `scarfinder!(psi, evolution, truncation, niter; refine=false, selector=nothing, kwargs...)`

Backend primitives:

- `tebd_evolve!(psi, gate, bond; maxdim, cutoff)`
- `tebd_evolution_from_hamiltonians(hamiltonians, dt; kwargs...)`
- `tebd_strang_evolution(nsites, dt; local_hamiltonian, kwargs...)`
- `tdvp_evolve!(psi::MPS, evo::TDVPEvolution)`
- `project!(psi, truncation)`
- `energy_density(psi, op)`
- `bond_entropy(psi, bond)`
- `entanglement_spectrum(psi, bond)`

Configuration types:

- `LocalGateEvolution`
- `TDVPEvolution`
- `BondDimTruncation`
- `EnergyTarget`
- `EntropySelector`
- `FidelitySelector`

Model helpers:

- `spinhalf_tfim_bond_hamiltonian`
- `spinhalf_xyz_bond_hamiltonian`

## Quick Start

### Finite OBC TEBD

```julia
using MPSToolkit, ITensors, ITensorMPS, LinearAlgebra

sites = siteinds("S=1/2", 6)
psi = MPS(sites, n -> "Up")

sx = ComplexF64[0 1; 1 0] / 2
hbond = kron(sx, sx)
gate = exp(-0.05im * hbond)

evo = LocalGateEvolution(gate, 0.05; schedule=1:5, nstep=1, maxdim=16, cutoff=1e-12)
trunc = BondDimTruncation(8; cutoff=1e-12)

scarfinder!(psi, evo, trunc, 2; refine=true, selector=EntropySelector())
println(expect(psi, "Sz"))
```

### Finite OBC TDVP

```julia
using MPSToolkit, ITensors, ITensorMPS

sites = siteinds("S=1/2", 6)
psi = MPS(sites, n -> "Up")

opsum = OpSum()
for j in 1:length(sites)
  opsum += 1.0, "Sx", j
end
hamiltonian = MPO(opsum, sites)

evo = TDVPEvolution(
  hamiltonian,
  -0.1im;
  time_step=-0.05im,
  nsteps=2,
  normalize=true,
  solver_kwargs=(maxdim=16, cutoff=1e-12),
)
trunc = BondDimTruncation(8; cutoff=1e-12)

scarfinder!(psi, evo, trunc, 1; refine=true, selector=EntropySelector())
println(expect(psi, "Sz"))
println(energy_density(psi, hamiltonian))
```

## Example Suite

The example tree is grouped by workflow:

- `examples/benchmarks/`
- `examples/daoe/`
- `examples/open_systems/`
- `examples/operator_space/`
- `examples/scarfinder/`
- `examples/tebd/`

Operator-space examples:

- [examples/daoe/projectors.jl](/Users/ren/Codex/MPSToolkit/examples/daoe/projectors.jl)
- [examples/operator_space/tfim_local_z.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_local_z.jl)
- [examples/operator_space/tfim_total_sz.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_total_sz.jl)
- [examples/operator_space/tfim_string.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_string.jl)
- [examples/operator_space/tfim_autocorrelator.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_autocorrelator.jl)
- [examples/operator_space/tfim_entanglement.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_entanglement.jl)
- [examples/operator_space/custom_hamiltonians.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/custom_hamiltonians.jl)
- [examples/operator_space/xyz_local_z.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/xyz_local_z.jl)
- [examples/open_systems/pauli_lindblad_tebd.ipynb](/Users/ren/Codex/MPSToolkit/examples/open_systems/pauli_lindblad_tebd.ipynb)
- [examples/open_systems/boundary_driven_xxz_steady_state.ipynb](/Users/ren/Codex/MPSToolkit/examples/open_systems/boundary_driven_xxz_steady_state.ipynb)

Additional workflows:

- [examples/chebyshev/two_peak_spectrum.jl](/Users/ren/Codex/MPSToolkit/examples/chebyshev/two_peak_spectrum.jl)
- [examples/chebyshev/two_peak_spectrum.ipynb](/Users/ren/Codex/MPSToolkit/examples/chebyshev/two_peak_spectrum.ipynb)
- [examples/chebyshev/tfim_local_spectrum.ipynb](/Users/ren/Codex/MPSToolkit/examples/chebyshev/tfim_local_spectrum.ipynb)
- [examples/tebd/xxz_tebd_vs_ed.ipynb](/Users/ren/Codex/MPSToolkit/examples/tebd/xxz_tebd_vs_ed.ipynb)
- [examples/benchmarks/pbc_tdvp_vs_tebd.ipynb](/Users/ren/Codex/MPSToolkit/examples/benchmarks/pbc_tdvp_vs_tebd.ipynb)
- [examples/tebd/disordered_xxz_mbl_tebd.ipynb](/Users/ren/Codex/MPSToolkit/examples/tebd/disordered_xxz_mbl_tebd.ipynb)
- [examples/scarfinder/xyz_spiral.jl](/Users/ren/Codex/MPSToolkit/examples/scarfinder/xyz_spiral.jl)

## Documentation

Additional docs:

- [docs/index.md](/Users/ren/Codex/MPSToolkit/docs/index.md)
- [docs/architecture.md](/Users/ren/Codex/MPSToolkit/docs/architecture.md)
- [docs/api.md](/Users/ren/Codex/MPSToolkit/docs/api.md)
- [docs/examples.md](/Users/ren/Codex/MPSToolkit/docs/examples.md)
