# TEBD From Hamiltonians Design

## Goal

Add a small convenience layer that builds `LocalGateEvolution` directly from dense local Hamiltonians so examples do not have to manually exponentiate gate lists.

## Scope

Add two public helpers:

- `local_gates_from_hamiltonians`
- `tebd_evolution_from_hamiltonians`

The helpers should support:

- one dense Hamiltonian reused across the schedule
- a vector of dense Hamiltonians
- a callable Hamiltonian provider
- an optional `map_hamiltonian` hook for cases such as operator-space TEBD

Then simplify the TFIM operator notebook to use the new evolution helper.

## Approach

Keep the helpers in the evolution layer. By default they exponentiate dense local Hamiltonians with `exp(-im * step_dt * h)`, but callers may override that conversion through `map_hamiltonian`.

The operator-space notebook should keep only TFIM local Hamiltonian construction and schedule setup. It should no longer manually build the gate list comprehension.

## Verification

Verification should include:

- targeted finite-TEBD tests for the new helper API
- package load verification
- notebook JSON parsing
