# API Reference

This page is the index for API-level documentation.

For most users, the best reading order is:

1. start with the relevant manual page
2. come here when you want the exact callable names and subsystem grouping

The workflow-specific manual pages carry most expanded source-backed docstrings:

- [TEBD And TDVP](manual/tebd-tdvp.md)
- [ScarFinder](manual/scarfinder.md)
- [Operator Space](manual/operator-space.md)
- [DAOE](manual/daoe.md)
- [DMT](manual/dmt.md)
- [Chebyshev](manual/chebyshev.md)

## Dispatch And Shared Entry Points

```@docs
evolve!
```

Also see `project!`, `scarfinder_step!`, `scarfinder!`, `MPSToolkit.ScarFinder.trajectory_refine!`, and `MPSToolkit.ScarFinder.match_energy!`.

## Evolution

- `LocalGateEvolution`
- `DMTGateEvolution`
- `TDVPEvolution`
- `tebd_evolve!`
- `dmt_step!`
- `dmt_evolve!`
- `tdvp_evolve!`
- `local_gates_from_hamiltonians`
- `tebd_evolution_from_hamiltonians`
- `tebd_strang_schedule`
- `tebd_strang_evolution`

## ScarFinder

- `BondDimTruncation`
- `EnergyTarget`
- `SelectionContext`
- `EntropySelector`
- `FidelitySelector`
- `MPSToolkit.score`

## Observables

```@docs
energy_density
bond_entropy
entanglement_spectrum
fidelity_distance
```

## Bases

- `pauli_matrices`
- `pauli_basis`
- `pauli_components`

## Operator Space

- `pauli_siteinds`
- `pauli_basis_state`
- `pauli_total_sz_state`
- `pauli_gate`
- `pauli_gate_from_hamiltonian`
- `pauli_lindblad_generator`
- `pauli_gate_from_lindbladian`
- `DMTOptions`
- `pauli_daoe_projector`
- `fdaoe_projector`

## Model Helpers

```@docs
spinhalf_matrices
spinhalf_xyz_bond_hamiltonian
spinhalf_tfim_bond_hamiltonian
```

## Chebyshev

- `ChebyshevRescaling`
- `SpectralFunction`
- `chebyshev_moments`
- `energy_cutoff!`
- `jackson_damping`
- `jackson_kernel`
- `reconstruct_chebyshev`
- `spectral_function`
