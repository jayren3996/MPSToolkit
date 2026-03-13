# Script-First Example Suite Design

## Goal

Replace the current verbose example set with a compact script-first suite centered on operator-space workflows, while pushing repeated helper code into the package.

## Scope

The new suite should:

- prioritize short `.jl` examples over notebooks
- go deeper on operator-space usage
- include one generic random-circuit TEBD example
- include one XYZ ScarFinder example using the user-specified couplings and target energy reference

The refactor should also add any reusable helpers exposed by the example pressure, especially around TEBD schedules, spin-1/2 model construction, and operator-space dynamics.

## Approach

Add a small helper layer to the package:

- Strang-schedule TEBD helpers
- spin-1/2 TFIM / XYZ model helpers
- a compact operator-space dynamics runner

Then replace the old examples with a new script set that keeps each file focused on one workflow and limits local scaffolding to problem-specific details only.

## Verification

Verification should include:

- updated unit tests for the new helper API
- successful execution of representative new examples
- successful package load
- documentation updates matching the new script set
