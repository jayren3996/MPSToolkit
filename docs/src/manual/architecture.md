# Architecture

## Public Namespaces

The root module is split into feature namespaces:

- `MPSToolkit.Evolution`
- `MPSToolkit.Observables`
- `MPSToolkit.ScarFinder`
- `MPSToolkit.Bases`
- `MPSToolkit.OperatorSpace`
- `MPSToolkit.Models`
- `MPSToolkit.Chebyshev`

## Shared Kernel

The ScarFinder loop is intentionally small:

1. evolve the state
2. project or truncate it
3. optionally apply target-energy correction
4. optionally refine by scanning a short trajectory with a selector

That split is why TEBD, TDVP, and projection-based workflows can share code without becoming one monolithic backend.

## Why The Evolution / Projection Split Matters

This package does not fuse time evolution and bond-dimension control into a single opaque routine.

That design keeps:

- evolution backends replaceable
- projection explicit and tunable
- targeting and refinement separate from the evolution engine
- lower-level stages directly callable for experiments

## Current State Support

Supported now:

- finite OBC `MPS`

Not supported now:

- general periodic-chain workflows as a first-class supported mode
