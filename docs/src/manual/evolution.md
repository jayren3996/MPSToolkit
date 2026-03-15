# Evolution

## TEBD

The TEBD layer is built around scheduled local-gate application:

- `LocalGateEvolution`
- `tebd_evolve!`
- `local_gates_from_hamiltonians`
- `tebd_evolution_from_hamiltonians`
- `tebd_strang_schedule`
- `tebd_strang_evolution`

The key design split is:

- the schedule decides where gates act and in what order
- the backend decides how each local update is applied and truncated

## TDVP

TDVP support is currently MPO-based:

- `TDVPEvolution`
- `tdvp_evolve!`

This lets you keep TDVP as a separate evolution backend while still plugging it into the same higher-level loops and diagnostics.

## DMT

DMT is currently implemented for operator-space workflows:

- `DMTGateEvolution`
- `dmt_step!`
- `dmt_evolve!`

Like TEBD, the scheduling layer is separate from the local update rule. The difference is that DMT applies an operator-space truncation strategy designed to preserve selected local information across each scheduled update.

## Projection And Refinement

Projection-oriented workflows are exposed directly:

- `BondDimTruncation`
- `project!`
- `EnergyTarget`
- `SelectionContext`
- `EntropySelector`
- `FidelitySelector`
- `scarfinder_step!`
- `scarfinder!`

This keeps ScarFinder-style iteration explicit rather than hard-wired into any one evolution engine.
