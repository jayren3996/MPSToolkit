# Operator-Space Helper API Design

## Goal

Move the reusable operator-space construction code out of notebooks/examples and into `MPSToolkit` so examples can stay short and focused.

## Scope

Add a small public helper layer for spin-1/2 operator-space work:

- `pauli_basis_state`
- `pauli_total_sz_state`
- `pauli_gate`
- `pauli_gate_from_hamiltonian`

Then rewrite the development TFIM notebook to use those helpers instead of reimplementing them locally.

## Approach

Keep the helpers narrowly scoped to the Pauli operator basis, with `(I, X, Y, Z)` as the default local ordering. The helpers should support the common case directly, with optional basis overrides available but off the main path.

The notebook should no longer define reusable basis transforms, Pauli-string product-state builders, or total-`S^z` constructors. It should only keep TFIM-specific bond Hamiltonian and schedule setup plus the short evolution runner.

## Verification

Verification should include:

- targeted operator-space tests for the new public helpers
- successful package loading in Julia
- successful JSON parsing of the updated notebook
