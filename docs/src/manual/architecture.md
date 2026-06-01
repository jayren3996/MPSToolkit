# Architecture

*How MPSToolkit is organized, and why the algorithms share so much without becoming one monolithic driver.*

MPSToolkit is built on a single design decision: **time evolution and bond-dimension control are kept separate, and both are public.** Almost every algorithm in the package — TEBD, TDVP, DMT, ScarFinder — is a small composition of the same handful of explicit building blocks. This page is the map that shows how those blocks fit together and where to find each one.

## Public namespaces

The root module re-exports everything, but it is also split into feature namespaces so the surface of each subsystem is easy to see:

| Namespace | Responsibility | Manual page | Source |
|---|---|---|---|
| `MPSToolkit.Evolution` | TEBD / TDVP / DMT configuration objects and the `evolve!` drivers | [TEBD](tebd.md), [TDVP](tdvp.md), [DMT](dmt.md) | `src/evolution/` |
| `MPSToolkit.Observables` | energy density, entanglement entropy and spectrum, fidelity distance | [API Reference](../api.md) | `src/observables/` |
| `MPSToolkit.ScarFinder` | projection, energy targeting, selectors, the ScarFinder loop | [ScarFinder](scarfinder.md) | `src/scarfinder/` |
| `MPSToolkit.OperatorSpace` | Pauli-basis states/gates, DMT, DAOE/FDAOE projectors | [Operator Space](operator-space.md), [DMT](dmt.md), [DAOE](daoe.md) | `src/operator_space/` |
| `MPSToolkit.Bases` | Pauli matrices, basis, and component helpers | [Operator Space](operator-space.md) | `src/bases/` |
| `MPSToolkit.Models` | dense spin-½ bond Hamiltonians for examples and constructors | [API Reference](../api.md) | `src/models/` |
| `MPSToolkit.Chebyshev` | Chebyshev moments, kernels, reconstruction, energy-window cutoff | [Chebyshev](chebyshev.md) | `src/chebyshev/` |

## The evolution / projection split

The package never fuses "advance the state in time" with "control the bond dimension" into one opaque routine. Instead it exposes a few **dispatch points** that different algorithms specialize:

- `evolve!(state, evolution)` advances a state under an evolution object. The object — `LocalGateEvolution` (TEBD), `TDVPEvolution` (TDVP), or `DMTGateEvolution` (DMT) — carries the schedule, gates or generator, and truncation budget. The method chosen depends only on the object's type.
- `project!(state, truncation)` applies an explicit truncation/projection, separate from any time step.
- `energy_density(state, operator)` and the other observables read out diagnostics without mutating the state.

Because these are ordinary, separately-callable functions, a higher-level workflow is free to interleave them in whatever order it needs.

## How the algorithms reuse the blocks

```
              ┌──────────────────────── evolve! ───────────────────────┐
              │                          │                              │
       LocalGateEvolution        TDVPEvolution                 DMTGateEvolution
       (TEBD: dense gates        (TDVP: MPO generator,         (DMT: dense gates,
        + Strang schedule         variational stepping)          identity-preserving
        + SVD truncation)                                         truncation)
              │                                                         │
              └───────────── same scheduling abstraction ──────────────┘

  ScarFinder  =  evolve!  →  project!  →  match_energy!  →  trajectory_refine!
                (any of the evolution objects above, composed explicitly)
```

- **[TEBD](tebd.md)** and **[DMT](dmt.md)** share the *same* scheduling idea — a list of bonds at which local gates act. They differ only in the post-gate truncation rule: TEBD uses ordinary SVD truncation, DMT uses a transport-preserving truncation. That is why they live side by side in `src/evolution/` and `src/operator_space/` rather than as unrelated engines.
- **[TDVP](tdvp.md)** is the variational alternative: instead of dense local gates it integrates the projected Schrödinger equation using an MPO generator. It plugs into the same `evolve!` entry point.
- **[ScarFinder](scarfinder.md)** is deliberately *not* a backend with its own internal state machine. It is a short loop — evolve, project, optionally match a target energy, optionally refine along a trajectory — built entirely from the public functions above. You can call every one of those pieces yourself.
- **[Operator Space](operator-space.md)** is the foundation under DMT, DAOE, and open-system work: the same TEBD scheduler runs in the vectorized Pauli basis once local Hamiltonians/Lindbladians are mapped to operator-space gates.
- **[Chebyshev](chebyshev.md)** is the one largely independent subsystem: it expands spectral functions in Chebyshev polynomials of a rescaled Hamiltonian, reusing the MPS representation but not the gate scheduler.

This is the payoff of the split: the evolution backends are replaceable, projection stays explicit and tunable, targeting and refinement are separate from the evolution engine, and every lower-level stage is directly callable for experiments and debugging.

## Conventions

- `evolve!`, `project!`, and `energy_density` are dispatch points — check the concrete method before changing behavior.
- Dense local operators are interpreted through the **local site dimension** of the input MPS.
- Operator-space code assumes the local Pauli ordering **`(I, X, Y, Z)`** unless a docstring says otherwise.
- Most user-facing routines **mutate** their `MPS`/state argument in place and return it.

## Current state support

Supported now:

- finite open-boundary-condition (OBC) `MPS`

Not supported as a first-class mode:

- general periodic-chain workflows (only narrow helper cases exist)
