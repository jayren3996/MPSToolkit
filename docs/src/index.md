# MPSToolkit.jl

*Composable finite matrix-product-state algorithms — evolution, projection, operator space, and spectra — in Julia.*

`MPSToolkit.jl` is a finite matrix-product-state (MPS) toolkit built on
[`ITensors.jl`](https://itensor.github.io/ITensors.jl/stable/) and `ITensorMPS.jl`. It is designed
for users who want explicit control over tensor-network workflows instead of one monolithic driver.
The important moving parts stay visible:

- TEBD is dense local gates plus an explicit schedule
- TDVP is an MPO generator plus explicit solver settings
- ScarFinder is a short loop over public evolution, projection, energy-matching, and selector routines
- operator-space workflows keep the Pauli-basis representation explicit
- Chebyshev tools expose both moment generation and spectral reconstruction

## 60-Second Demo

Build a Strang-split TEBD schedule from a bond Hamiltonian, evolve a Néel state, and read off two
diagnostics — every ingredient is a separate, public object:

```julia
using MPSToolkit, ITensors, ITensorMPS

# A 12-site transverse-field Ising chain, started in a Néel product state.
nsites = 12
sites  = siteinds("S=1/2", nsites)
psi    = MPS(sites, n -> isodd(n) ? "Up" : "Dn")

# TEBD is an explicit object: an odd–even–odd Strang schedule of local gates built from the
# bond Hamiltonian. Here, 200 sweeps of dt = 0.05 (total evolution time t = 10).
evolution = tebd_strang_evolution(
    nsites, 0.05;
    local_hamiltonian = (bond, weight) ->
        weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J = 1.0, g = 0.8),
    nstep  = 200,
    maxdim = 32,
    cutoff = 1e-12,
)

evolve!(psi, evolution)         # mutates psi in place and returns it

expect(psi, "Sz")               # staggered order has relaxed toward 0
bond_entropy(psi, nothing)      # half-chain entanglement has grown from 0 to ≈ 1.47
```

## What MPSToolkit Is Good For

- build finite-chain TEBD workflows from dense local Hamiltonians
- run MPO-based TDVP on finite open-boundary `MPS`
- study explicit projection-and-refinement loops such as ScarFinder
- work in Pauli operator space for open-system or projector-based calculations
- compute Chebyshev moments and reconstruct spectral functions

## Choose Your Path

| If you want to … | Start with |
|---|---|
| Build and run a finite-chain TEBD workflow | [Getting Started](getting-started.md) → [TEBD and TDVP](manual/tebd-tdvp.md) |
| Run MPO-based TDVP on a finite OBC `MPS` | [TEBD and TDVP](manual/tebd-tdvp.md) |
| Study explicit projection-and-refinement loops | [ScarFinder](manual/scarfinder.md) |
| Work in Pauli operator space | [Operator Space](manual/operator-space.md) |
| Truncate or dissipate operator space | [DMT](manual/dmt.md) · [DAOE](manual/daoe.md) |
| Reconstruct spectral functions | [Chebyshev](manual/chebyshev.md) |

## Design Philosophy

The package is intentionally modular. Higher-level workflows are thin wrappers around reusable
lower-level functions, so the same pieces can be recombined for custom projects, experiments, and
debugging. The manual is therefore organized by subsystem rather than by one giant "run everything"
interface — if you already know what kind of computation you want, jump straight to the matching
manual page.

## Current Limits

- finite OBC `MPS` is the main supported state class
- periodic chains are not a general supported mode
- `TDVPEvolution` currently expects MPO generators
- DMT is currently implemented for operator-space workflows

## Navigation

- [Getting Started](getting-started.md) — minimal setup, installation, and first examples
- [Architecture](manual/architecture.md) — how the codebase is organized and how the building blocks relate
- [TEBD and TDVP](manual/tebd-tdvp.md) — dense-gate TEBD helpers and MPO-based TDVP
- [ScarFinder](manual/scarfinder.md) — explicit projection workflows, selector refinement, and step-count guidance
- [Operator Space](manual/operator-space.md) — Pauli-basis states, gates, and operator-space evolution
- [DAOE](manual/daoe.md) — DAOE and FDAOE projector construction
- [DMT](manual/dmt.md) — operator-space DMT truncation and scheduling
- [Chebyshev](manual/chebyshev.md) — moment generation, energy-window cutoff, and spectral reconstruction
- [Examples](examples.md) — script and notebook entry points
- [API Reference](api.md) — source-backed API docs grouped by subsystem
