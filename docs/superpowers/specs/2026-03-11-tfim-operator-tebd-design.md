# TFIM Operator TEBD Design

## Goal

Add a development notebook that evolves the total spin operator `S^z = (1/2) \sum_j \sigma_j^z` under the open-boundary TFIM using TEBD in the operator-space basis.

## Scope

The notebook should:

- create Pauli-basis sites with `pauli_siteinds`
- build the initial operator state for total `S^z`
- construct two-site operator-space TEBD gates for the TFIM
- run time evolution with `LocalGateEvolution` and `evolve!`
- expose the core dynamics data so additional analysis can be added later

The notebook should not add reusable package APIs, plotting, or exact-diagonalization benchmarking.

## Approach

Keep the implementation notebook-local in `dev/`:

- use the local operator basis ordering `(I, X, Y, Z)`
- represent the initial operator as an operator-space `MPS`
- build each two-site gate by:
  - defining the physical two-site TFIM bond Hamiltonian
  - exponentiating it to a physical gate
  - converting the induced superoperator `O -> U O U^\dagger` into the Pauli operator basis
- assemble a Strang-style odd-even-odd TEBD schedule with bond-dependent field weights so the open-chain field term is counted once per site

## Verification

Verification should include:

- parsing the notebook JSON successfully
- loading the package in Julia
- running a Julia smoke test that builds the operator-space evolution object and advances the total-`S^z` operator by one step
