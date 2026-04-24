# API

## Entry Points

### `scarfinder_step!(psi, evolution, truncation; target_energy=nothing, selector=nothing, kwargs...)`

Runs one ScarFinder step:

1. `evolve!`
2. `project!`
3. optional `match_energy!`

### `scarfinder!(psi, evolution, truncation, niter; refine=false, selector=nothing, kwargs...)`

Runs `niter` ScarFinder steps, then optionally performs trajectory refinement with the provided selector and selector context.

## Evolution Configurations

### `LocalGateEvolution`

Fields:

- `gate`
- `dt`
- `schedule`
- `nstep`
- `maxdim`
- `cutoff`

Used by:

- `evolve!(::MPS, ::LocalGateEvolution)`
- `tebd_evolution_from_hamiltonians`

Notes:

- `gate` may be a single dense local matrix, a per-schedule vector of dense gates, or a callable gate provider
- single-site gates are allowed as long as their matrix dimension matches the local site dimension
- `schedule` is the ordered list of left-bond indices to update
- `nstep` repeats one full pass over `schedule`
- `maxdim` and `cutoff` are passed to the gate application/compression routine

### `TDVPEvolution`

Fields:

- `generator`
- `t`
- `time_step`
- `nsteps`
- `nsweeps`
- `reverse_step`
- `updater_backend`
- `updater`
- `normalize`
- `solver_kwargs`
- `schedule`

Used by:

- `evolve!(::MPS, ::TDVPEvolution)`
- `tdvp_evolve!(::MPS, ::TDVPEvolution)`

Notes:

- `generator` is currently expected to be an MPO
- `t` is the total evolution time passed to `ITensorMPS.tdvp`
- `time_step` controls the per-step TDVP increment
- if `nsteps` is `nothing`, the wrapper falls back to `nsweeps`
- `solver_kwargs` is passed through to `ITensorMPS.tdvp`

## Truncation

### `BondDimTruncation`

Fields:

- `maxdim`
- `cutoff`

Used by:

- `project!(::MPS, ::BondDimTruncation)`

Notes:

- this is the explicit post-evolution compression stage used by the ScarFinder loop
- for TDVP workflows it lets you keep the evolution backend and the ScarFinder truncation policy separate

## Targeting And Selection

### `EnergyTarget`

Fields:

- `target`
- `operator`
- `tol`
- `alpha`
- `maxstep`

Notes:

- `operator` may be a dense local operator or an MPO
- MPO targets are enforced with an internal TDVP correction step followed by the usual ScarFinder projection

### `SelectionContext`

Fields:

- `reference_state`

Notes:

- carries optional external data into selector scoring
- `FidelitySelector` uses `reference_state`

### `EntropySelector`

Fields:

- `bond`

Scoring is implemented in:

- `src/scarfinder/selectors.jl`

### `FidelitySelector`

Scores a state by its fidelity distance to `SelectionContext(reference_state=...)`.

## Backend Primitives

Exported helpers:

- `tebd_evolve!`
- `local_gates_from_hamiltonians`
- `tebd_evolution_from_hamiltonians`
- `tebd_strang_schedule`
- `tebd_strang_evolution`
- `tdvp_evolve!`
- `project!`
- `energy_density`
- `bond_entropy`
- `entanglement_spectrum`
- `fidelity_distance`
- `dmt_step!`
- `dmt_evolve!`

## Basis And Operator-Space Helpers

Exported helpers:

- `pauli_matrices`
- `pauli_basis`
- `pauli_components`
- `pauli_siteinds`
- `pauli_basis_state`
- `pauli_total_sz_state`
- `pauli_gate`
- `pauli_gate_from_hamiltonian`
- `pauli_lindblad_generator`
- `pauli_gate_from_lindbladian`
- `pauli_daoe_projector`
- `fdaoe_projector`

## Model Helpers

Exported helpers:

- `spinhalf_matrices`
- `spinhalf_tfim_bond_hamiltonian`
- `spinhalf_xyz_bond_hamiltonian`

## Chebyshev Helpers

Exported helpers:

- `ChebyshevRescaling`
- `SpectralFunction`
- `chebyshev_moments`
- `jackson_damping`
- `jackson_kernel`
- `reconstruct_chebyshev`
- `spectral_function`

Notes:

- `chebyshev_moments` assumes the input Hamiltonian has already been rescaled into the Chebyshev window `[-1, 1]`
- `spectral_function` stores physical frequency metadata separately from the raw moment sequence
- the first implementation keeps moment generation and spectrum reconstruction separate on purpose, so later post-processing methods can plug in cleanly

These remain public on purpose, so the package can be tuned below the level of the high-level ScarFinder loop.
