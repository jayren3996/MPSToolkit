# Operator-Space Examples Design

## Goal

Add a dedicated DAOE/FDAOE example pair that teaches the operator-space API with one runnable Julia script and one explanatory notebook.

## Scope

The new example pair should:

- build Pauli-basis sites with `pauli_siteinds`
- construct `pauli_daoe_projector` and `fdaoe_projector`
- build Pauli-basis product states in the local `(I, X, Y, Z)` basis
- apply the projector MPOs with the standard `apply` path
- show preserved versus damped strings with explicit expected overlaps
- explain the meaning of `lstar`, `wstar`, and `gamma`

The example pair should not add a DAOE-specific evolution workflow. It should stay focused on projector construction and application.

## Approach

Use one shared tutorial for both DAOE and FDAOE:

- `examples/operator_space_daoe.jl`
- `examples/notebooks/06_operator_space_daoe.ipynb`

This keeps the operator-space API in one place, minimizes duplicated helper code, and makes the difference between Pauli-weight damping and fermionic-weight damping easy to compare.

## Example Content

The script and notebook should both cover:

1. `pauli_siteinds`
2. `pauli_daoe_projector`
3. `fdaoe_projector`
4. basis-product operator states in the `(I, X, Y, Z)` basis
5. projector application through `apply`
6. overlap checks against known damping factors

Recommended example strings:

- DAOE with `lstar = 1`
  - `X I I I` is preserved
  - `X Y Z I` is damped by `exp(-2gamma)`
- FDAOE with `wstar = 2`
  - `Z I I I` is preserved
  - `X X I I` is preserved
  - `Z Z I I` is damped by `exp(-2gamma)`
  - `X Z Z X` is preserved because the middle `Z` operators are inside the Jordan-Wigner tail

## Documentation Updates

Update:

- `README.md`
- `docs/examples.md`

The docs should add the new script and notebook to the curated example index and remove outdated placeholder language around operator-space functionality.

## Verification

Verification should include:

- running the new script with `julia --project=. examples/operator_space_daoe.jl`
- parsing the new notebook JSON
- running the full package test suite with `julia --project=. -e 'using Pkg; Pkg.test()'`
