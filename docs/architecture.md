# Architecture

## Public Namespaces

The root module is split into feature namespaces:

- `MPSToolkit.Evolution`
- `MPSToolkit.Observables`
- `MPSToolkit.ScarFinder`
- `MPSToolkit.Bases`
- `MPSToolkit.OperatorSpace`

## Shared Kernel

The core loop lives in:

- [src/scarfinder/algorithm.jl](/Users/ren/Codex/MPSToolkit/src/scarfinder/algorithm.jl)

Its contract is intentionally small:

1. evolve the state
2. project/truncate the state
3. optionally apply target-energy correction
4. optionally refine by scanning a short trajectory with a selector

The kernel does not know which finite evolution backend it is driving. It only dispatches on `evolution` and `truncation`.

## Backend Surface

Evolution and observable files:

- [src/evolution/tebd.jl](/Users/ren/Codex/MPSToolkit/src/evolution/tebd.jl)
- [src/evolution/tdvp.jl](/Users/ren/Codex/MPSToolkit/src/evolution/tdvp.jl)
- [src/observables/entanglement.jl](/Users/ren/Codex/MPSToolkit/src/observables/entanglement.jl)
- [src/observables/energy.jl](/Users/ren/Codex/MPSToolkit/src/observables/energy.jl)

Shared projection dispatch:

- [src/scarfinder/dispatch.jl](/Users/ren/Codex/MPSToolkit/src/scarfinder/dispatch.jl)

## Why The Evolution/Projection Split Matters

This package does not fuse time evolution and bond-dimension control into one opaque routine.

That split is deliberate:

- TEBD and TDVP can drive the same ScarFinder loop
- projection remains explicit and tunable
- energy targeting and refinement stay separate from the evolution engine
- internal stages remain directly callable for experiments

## State Support

Supported now:

- finite OBC `MPS`

Not supported now:

- finite-ring PBC `MPS`
