# MPSToolkit Design

**Goal:** Rename the finite-only package from `ScarFinder` to `MPSToolkit` and broaden it into a finite-chain tensor-network toolkit with evolution, observables, generic selection-driven ScarFinder, basis helpers, and a future operator-space extension point.

## Scope

The package will remain finite-`MPS` only. The immediate supported feature areas are:

- time evolution from local gates and from Hamiltonians
- entanglement entropy and entanglement spectrum
- ScarFinder
- Pauli/spin-1/2 basis helpers
- an `OperatorSpace` namespace reserved for later DAOE work

## Package Structure

The package root becomes `MPSToolkit`, with submodules:

- `MPSToolkit.Evolution`
- `MPSToolkit.Observables`
- `MPSToolkit.ScarFinder`
- `MPSToolkit.Bases`
- `MPSToolkit.OperatorSpace`

The root re-exports the high-use API, but the submodules provide stable organizational boundaries for future growth.

## Selector Design

The old selector model is generalized into a scoring interface:

- `score(selector, psi, context)`

The refinement loop keeps the state with the best score. The default comparison is lower-is-better.

Initial selector types:

- `EntropySelector`
- `FidelitySelector`

The selector `context` carries any required external information such as the initial state or reference observables.

## Migration

The rename is a hard rename:

- package name becomes `MPSToolkit`
- module name becomes `MPSToolkit`
- examples, tests, notebooks, and docs switch from `using ScarFinder` to `using MPSToolkit`

No compatibility alias for `ScarFinder` is retained.

## Implementation Strategy

Use the current finite backend logic as the implementation base. Avoid unnecessary algorithm rewrites while restructuring source layout and public naming. Add the new selector API first, then refactor exports/modules around it, then update docs/examples.
