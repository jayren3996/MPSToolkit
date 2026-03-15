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

These tools are intended for explicit operator-space workflows rather than as hidden internals.

## Examples

For runnable examples, see [Examples](../examples.md), especially the operator-space and open-system notebooks.
