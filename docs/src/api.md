# API Reference

This page is the index for API-level documentation.

For most users, the best reading order is:

1. start with the relevant manual page
2. come here when you want the exact callable names, signatures, and source-backed docstrings

The workflow-specific manual pages carry the main narrative API docs:

- [TEBD And TDVP](manual/tebd-tdvp.md)
- [ScarFinder](manual/scarfinder.md)
- [Operator Space](manual/operator-space.md)
- [DAOE](manual/daoe.md)
- [DMT](manual/dmt.md)
- [Chebyshev](manual/chebyshev.md)

This page keeps the cross-cutting entry points and reusable helpers in one place.

## Dispatch And Shared Entry Points

```@docs
evolve!
```

## Observables

```@docs
energy_density
bond_entropy
entanglement_spectrum
fidelity_distance
```

## Model Helpers

```@docs
spinhalf_matrices
spinhalf_xyz_bond_hamiltonian
spinhalf_tfim_bond_hamiltonian
```
