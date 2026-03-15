# API Reference

This page focuses on the documented public surface of `MPSToolkit.jl`.

## Evolution Types

```@docs
LocalGateEvolution
DMTGateEvolution
TDVPEvolution
BondDimTruncation
EnergyTarget
SelectionContext
EntropySelector
FidelitySelector
```

## Evolution Entry Points

```@docs
evolve!
tebd_evolve!
dmt_step!
dmt_evolve!
local_gates_from_hamiltonians
tebd_evolution_from_hamiltonians
tebd_strang_schedule
tebd_strang_evolution
tdvp_evolve!
project!
```

## ScarFinder And Selection

```@docs
scarfinder_step!
scarfinder!
```

## Observables

```@docs
energy_density
bond_entropy
entanglement_spectrum
fidelity_distance
```

## Basis And Operator-Space Helpers

```@docs
pauli_matrices
pauli_basis
pauli_components
pauli_siteinds
pauli_basis_state
pauli_total_sz_state
pauli_gate
pauli_gate_from_hamiltonian
pauli_lindblad_generator
pauli_gate_from_lindbladian
pauli_daoe_projector
fdaoe_projector
```

## Model Helpers

```@docs
spinhalf_matrices
spinhalf_xyz_bond_hamiltonian
spinhalf_tfim_bond_hamiltonian
```

## Chebyshev

```@docs
ChebyshevRescaling
SpectralFunction
chebyshev_moments
energy_cutoff!
jackson_damping
jackson_kernel
reconstruct_chebyshev
spectral_function
```
