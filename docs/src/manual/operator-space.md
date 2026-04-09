# Operator Space

`MPSToolkit.jl` includes a Pauli-basis operator-space layer for dynamics and open-system workflows.

## Basis Helpers

The core basis and state helpers are:

- `pauli_siteinds`
- `pauli_basis_state`
- `pauli_total_sz_state`
- `pauli_matrices`
- `pauli_basis`
- `pauli_components`

The Pauli basis is the shared foundation for operator-space TEBD, DMT, Lindblad evolution, and DAOE-style projectors. In practice, a typical workflow is:

1. build Pauli-basis site indices with `pauli_siteinds`
2. prepare an operator-state `MPS`
3. map local Hamiltonians or Lindbladians into Pauli-space gates
4. run scheduled evolution or projection

## Minimal Example

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
```

The helper notebooks [operator_tebd_helper_apis.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/operator_tebd_helper_apis.ipynb) and [dmt_scheduler.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/dmt_scheduler.ipynb) show the fuller scheduler-driven workflow.

## Basis And Mapping API

```@docs
pauli_siteinds
pauli_basis_state
pauli_total_sz_state
pauli_matrices
pauli_basis
pauli_components
pauli_gate
pauli_gate_from_hamiltonian
pauli_lindblad_generator
pauli_gate_from_lindbladian
```

## Local Operator And Lindblad Maps

The package exposes helpers for converting local Hamiltonians and Lindbladians into operator-space gates:

- `pauli_gate`
- `pauli_gate_from_hamiltonian`
- `pauli_lindblad_generator`
- `pauli_gate_from_lindbladian`

## DMT And Projectors

Operator-space specific algorithms include:

- `dmt_step!`
- `dmt_evolve!`
- `DMTGateEvolution`
- `pauli_daoe_projector`
- `fdaoe_projector`

!!! note "DMT is transport-specific"
    The DMT functions (`dmt_step!`, `dmt_evolve!`, `DMTGateEvolution`) implement a
    truncation scheme designed specifically for **transport** (e.g. spin or energy diffusion).
    Their truncation rule preserves the identity component at every bond — an assumption that
    holds for transport but not for arbitrary operator-space tasks.  For general operator-space
    TEBD evolution, use [`LocalGateEvolution`](@ref) with the Pauli-basis helpers above.

These tools are intended for explicit operator-space workflows rather than as hidden internals.

## Examples

For runnable examples, see [Examples](../examples.md), especially the operator-space and open-system notebooks.
