# MPSToolkit DAOE/FDAOE Design

**Goal:** Add operator-space projection support for both DAOE and FDAOE by constructing projector MPOs that can be applied through the existing MPO machinery.

## Scope

This feature adds:

- DAOE projector MPO construction for spin-1/2 operator-space work
- FDAOE projector MPO construction for fermionic operator-space work
- operator-weight diagnostics for the corresponding operator bases

This feature does not add:

- a new monolithic DAOE execution function
- DAOE-specific evolution drivers
- Lindbladian/open-system evolution wrappers

## API Design

The public API should stay minimal and compositional:

- `pauli_daoe_projector(...)`
- `fdaoe_projector(...)`
- operator-space diagnostics such as weight-profile helpers

Usage is intentionally based on existing MPO application paths:

```julia
P = pauli_daoe_projector(...)
O = apply(P, O; maxdim=..., cutoff=...)
```

The same pattern applies to FDAOE.

## Module Placement

The functionality belongs under `OperatorSpace`, not under `Evolution` or `ScarFinder`.

Responsibilities:

- `OperatorSpace`: projector MPO builders and Pauli-basis diagnostics
- `Evolution`: physical time evolution of states or operators
- `ScarFinder`: projected state-search workflow for MPS states

## Configuration

There is no need for an extra runtime wrapper type if the projector is directly represented as an MPO. The required controls are the projector construction parameters and the standard MPO application compression controls already supported by ITensors/ITensorMPS.

Avoid exposing irrelevant parameters such as:

- `dt`
- `normalize`
- numerical “strategy” flags

Those belong either to evolution or to future internal implementation details.

## Testing Strategy

Tests should cover:

- projector MPO construction returns well-formed MPOs
- applying the projector MPO to simple test operators works without shape/index errors
- small analytic operator examples show expected suppression or retention patterns
- operator-weight diagnostics match known decompositions for simple inputs

Where possible, compare small cases against the EDKit reference implementation or exact small-system calculations.
