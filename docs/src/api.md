# API Reference

This page is the index for API-level documentation.

For most users, the best reading order is:

1. start with the relevant manual page
2. come here when you want the exact callable names and subsystem grouping

The workflow-specific manual pages carry most expanded source-backed docstrings:

- [TEBD](manual/tebd.md)
- [TDVP](manual/tdvp.md)
- [ScarFinder](manual/scarfinder.md)
- [Operator Space](manual/operator-space.md)
- [DMT](manual/dmt.md)
- [DAOE](manual/daoe.md)
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
- `operator_basis_matrices` — generic-`d` generalized Gell-Mann basis; `pauli_matrices`/`pauli_basis` are the `d = 2` case (see [Operator Space](manual/operator-space.md)).
- `local_dimension` — recovers `d` from an operator-space site or its dimension `d^2`.

## Operator Space

Generic in the local Hilbert space dimension `d` (`operator_*`), with every `pauli_*` name as
the `d = 2` case (see [Operator Space](manual/operator-space.md) for the basis and conventions,
and [DMT](manual/dmt.md) for a `d = 3` worked example):

- `operator_siteinds` / `pauli_siteinds`
- `operator_basis_state` / `pauli_basis_state`
- `operator_product_state` — literal tensor-product state builder; no `pauli_*` wrapper needed (not `d`-specific in convention).
- `operator_local_sum_state` — literal local-density sum; no `pauli_*` wrapper needed.
- `pauli_total_sz_state`, `pauli_domain_wall_state` — `d = 2`-only conveniences with their own bond-dimension-2 normalization; use `operator_local_sum_state` directly at other `d`.
- `operator_gate` / `pauli_gate`
- `operator_gate_from_hamiltonian` / `pauli_gate_from_hamiltonian`
- `operator_lindblad_generator` / `pauli_lindblad_generator`
- `operator_gate_from_lindbladian` / `pauli_gate_from_lindbladian`
- `operator_gate_from_imaginary_time` / `pauli_gate_from_imaginary_time`
- `operator_state_from_mpo` / `pauli_state_from_mpo`
- `operator_superoperator_mpo` / `pauli_superoperator_mpo`
- `operator_gibbs_state` / `pauli_gibbs_state`
- `operator_trace` / `pauli_trace`
- `operator_expectation` / `pauli_expectation`
- `operator_expectation_profile` / `pauli_expectation_profile`
- `pauli_pxp_constraint_state` — **spin-1/2 only**, raises `ArgumentError` at other `d`.
- `pauli_pxp_constraint_projector` — **spin-1/2 only**, raises `ArgumentError` at other `d`.
- `constrained_dmt_evolve!`
- `DMTOptions` — generic in `d`; see [DMT](manual/dmt.md) for `preserve_diameter` and the budget semantics (`connector_buffer` was removed).
- `pauli_daoe_projector` — **spin-1/2 only**, raises `ArgumentError` at other `d`.
- `pauli_fdaoe_projector` (alias: `fdaoe_projector`) — **spin-1/2 only**, raises `ArgumentError` at other `d`.

## Model Helpers

```@docs
spinhalf_matrices
spinhalf_xyz_bond_hamiltonian
spinhalf_tfim_bond_hamiltonian
```

PXP helpers (documented on the [DMT](manual/dmt.md) page):

- `pxp_term_hamiltonian`
- `pxp_term_support`
- `pxp_constraint_mpo`

## Chebyshev

- `ChebyshevRescaling`
- `SpectralFunction`
- `chebyshev_moments`
- `energy_cutoff!`
- `jackson_damping`
- `jackson_kernel`
- `reconstruct_chebyshev`
- `spectral_function`
