<div align="center">

<img src="docs/src/assets/logo.svg" alt="MPSToolkit.jl logo" width="170"/>

# MPSToolkit.jl

**Composable finite matrix-product-state algorithms in Julia.**

Time evolution, projection, operator space, and Chebyshev spectra — explicit building blocks, not one opaque driver.

[![Docs](https://img.shields.io/badge/docs-dev-9558B2.svg)](https://jayren3996.github.io/MPSToolkit/dev/) [![CI](https://github.com/jayren3996/MPSToolkit/actions/workflows/documentation.yml/badge.svg)](https://github.com/jayren3996/MPSToolkit/actions/workflows/documentation.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Julia](https://img.shields.io/badge/Julia-1.10%2B-389826.svg)](https://julialang.org)

</div>

---

MPSToolkit.jl is a finite matrix-product-state (MPS) toolkit built on
[`ITensors.jl`](https://itensor.github.io/ITensors.jl/stable/) and
[`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl). It is organized around explicit,
reusable tensor-network building blocks, so time evolution, projection, observables,
operator-space tools, and spectral routines all stay directly callable and recombinable.

## ✨ Features

|  |  |
| --- | --- |
| 🧩 **Composable building blocks** | TEBD, TDVP, projection, operator space, and spectra are separate, public objects you recombine — swap any one without rewriting the others. |
| ⏱ **Dense-gate TEBD** | Explicit odd–even–odd Strang schedules from local Hamiltonians via `tebd_strang_evolution` / `tebd_strang_schedule`. |
| 🌀 **MPO-based TDVP** | Finite open-boundary `MPS` evolution through `TDVPEvolution`. |
| 🎯 **ScarFinder workflows** | Explicit projection (`project!`), energy targeting (`match_energy!`), and selector-driven refinement (`scarfinder!`, `trajectory_refine!`). |
| ⚛️ **Operator space, any local dimension** | States and gates (`operator_siteinds`, `operator_basis_state`, `operator_gate_from_hamiltonian`) for coherent and open-system evolution; `pauli_*` names are the spin-1/2 case. |
| 🧪 **Operator-space truncation** | Density-matrix truncation (DMT) and DAOE / FDAOE dissipative projectors. |
| 📈 **Chebyshev spectra** | Moments, energy-window cutoff, Jackson damping, and spectral-function reconstruction. |
| 📐 **Diagnostics** | `energy_density`, `bond_entropy`, `entanglement_spectrum`, and `fidelity_distance`. |

## 📦 Installation

```julia
pkg> add https://github.com/jayren3996/MPSToolkit.git
```

Then load it with `using MPSToolkit`.

## 🚀 Quick Start

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

The schedule, the gates, the evolution object, and the diagnostics are all separate pieces — swap
any one without rewriting the others.

## 🧭 Choose Your Path

| If you want to … | Start with |
| --- | --- |
| Build and run a finite-chain TEBD workflow | [Getting Started](https://jayren3996.github.io/MPSToolkit/dev/getting-started/) → [TEBD](https://jayren3996.github.io/MPSToolkit/dev/manual/tebd/) |
| Run MPO-based TDVP on a finite OBC `MPS` | [TDVP](https://jayren3996.github.io/MPSToolkit/dev/manual/tdvp/) |
| Study explicit projection-and-refinement loops | [ScarFinder](https://jayren3996.github.io/MPSToolkit/dev/manual/scarfinder/) |
| Work in Pauli operator space | [Operator Space](https://jayren3996.github.io/MPSToolkit/dev/manual/operator-space/) |
| Truncate or dissipate operator space | [DMT](https://jayren3996.github.io/MPSToolkit/dev/manual/dmt/) · [DAOE](https://jayren3996.github.io/MPSToolkit/dev/manual/daoe/) |
| Reconstruct spectral functions | [Chebyshev](https://jayren3996.github.io/MPSToolkit/dev/manual/chebyshev/) |

## 🛠 Usage

<details>
<summary><b>Explicit ScarFinder projection loop</b></summary>

<br>

ScarFinder is a thin wrapper around public building blocks: it alternates a TEBD step with an
explicit bond-truncation projection. Nothing is hidden — you build the schedule, the truncation,
and the loop yourself.

```julia
using MPSToolkit, ITensors, ITensorMPS

sites = siteinds("S=1/2", 10)
psi   = MPS(sites, n -> isodd(n) ? "Up" : "Dn")

# A periodic odd–even–odd schedule on the 10-site ring; entry `10` is the boundary bond (10, 1).
schedule = [1, 3, 5, 7, 9, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9]
weights  = [fill(0.5, 5); fill(1.0, 5); fill(0.5, 5)]   # Strang prefactors
bond_hamiltonians = [w * spinhalf_xyz_bond_hamiltonian(; Jx = 0.88, Jy = 1.0, Jz = 0.29)
                     for w in weights]

evolution  = tebd_evolution_from_hamiltonians(bond_hamiltonians, 0.1;
                                              schedule = schedule, maxdim = 1, cutoff = 1e-12)
truncation = BondDimTruncation(1; cutoff = 1e-12)       # rank-1 projection on every bond

# 150 projected TEBD sweeps: evolve, then project, repeat.
scarfinder!(psi, evolution, truncation, 150; refine = false)

bond_entropy(psi, 5)
```

</details>

<details>
<summary><b>Operator-space autocorrelator</b></summary>

<br>

TEBD runs just as well in Pauli operator space: build the site indices with `pauli_siteinds`,
map bond Hamiltonians to gates with `pauli_gate_from_hamiltonian`, and record the full
autocorrelation trace `⟨O(0), O(t)⟩` along the trajectory.

```julia
using MPSToolkit, ITensors, ITensorMPS

nsites = 8
sites  = pauli_siteinds(nsites)

# TEBD in Pauli operator space: gates come from the TFIM bond Hamiltonian.
evolution = tebd_strang_evolution(
    nsites, 0.04;
    local_hamiltonian = (bond, weight) ->
        weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J = 1.0, g = 1.05),
    map_hamiltonian   = pauli_gate_from_hamiltonian,
    maxdim = 96, cutoff = 1e-12,
)

local_z = pauli_basis_state(sites, ["I", "I", "I", "Z", "I", "I", "I", "I"])
op      = copy(local_z)

trace = ComplexF64[inner(local_z, op)]
for _ in 1:10
    evolve!(op, evolution)
    push!(trace, inner(local_z, op))   # ⟨O(0), O(t)⟩ at each step
end
trace
```

</details>

## 🗂 Repository Map

MPSToolkit prefers explicit workflows over hidden orchestration. High-level routines like
ScarFinder are deliberately thin wrappers around public building blocks (`evolve!`, `project!`,
`match_energy!`, and the selector scoring helpers), and most user-facing routines mutate their
`MPS`/state argument in place and return it.

| Directory | What lives there |
| --- | --- |
| [`src/evolution/`](src/evolution) | TEBD and TDVP configuration types (`LocalGateEvolution`, `DMTGateEvolution`, `TDVPEvolution`) and the concrete `evolve!` methods. |
| [`src/scarfinder/`](src/scarfinder) | Projection settings, selector types, the explicit ScarFinder loop, and post-step energy matching. |
| [`src/observables/`](src/observables) | Energy density, entanglement entropy and spectrum, and fidelity-style diagnostics. |
| [`src/operator_space/`](src/operator_space) | Pauli-basis helpers, MPO vectorizers, thermal/domain-wall preparation, expectation sweeps, DMT (plus constraint-checkpointed DMT), and DAOE / FDAOE projectors. |
| [`src/chebyshev/`](src/chebyshev) | Chebyshev moments, energy-window projection, damping kernels, and spectral reconstruction. |
| [`src/models/`](src/models), [`src/bases/`](src/bases) | Small dense helper matrices used by examples and constructors. |

The public surface is grouped into feature namespaces — `MPSToolkit.Evolution`, `.Observables`,
`.ScarFinder`, `.OperatorSpace`, `.Bases`, `.Models`, and `.Chebyshev` — and every name is also
exported from the root module.

## 📐 Conventions Worth Knowing

- `evolve!`, `project!`, and `energy_density` are **dispatch points** — check the concrete methods before changing behavior.
- Dense local operators are interpreted through the **local site dimension** of the input `MPS`.
- Operator-space code vectorizes onto a generalized Gell-Mann basis, identity first, generic in the local Hilbert space dimension `d` (`operator_*`); at `d = 2` this is exactly the Pauli ordering **`(I, X, Y, Z)`** that every `pauli_*` name uses, unless a docstring says otherwise.
- **Finite OBC `MPS`** is the main supported state class; periodic behavior is limited to a few helper cases.
- **ScarFinder step-count guard.** `scarfinder!`, `scarfinder_step!`, and `trajectory_refine!` treat an effective evolution step count of `1` as a misconfiguration: they emit a warning and internally use `10` for that call. This rule is local to ScarFinder — global TEBD, DMT, and TDVP constructors keep their own defaults, and the energy-correction substeps inside `match_energy!` intentionally remain single-step.

### ⚠️ Current limits

- finite OBC `MPS` is the main supported state class
- periodic chains are not a general supported mode
- `TDVPEvolution` currently expects MPO generators
- DMT is currently implemented for operator-space workflows

## 📚 Documentation

Full documentation lives at **[jayren3996.github.io/MPSToolkit](https://jayren3996.github.io/MPSToolkit/dev/)**.

- [Getting Started](https://jayren3996.github.io/MPSToolkit/dev/getting-started/) — the first finite-chain TEBD workflow.
- [Architecture](https://jayren3996.github.io/MPSToolkit/dev/manual/architecture/) — how the building blocks fit together.
- [API Reference](https://jayren3996.github.io/MPSToolkit/dev/api/) — signatures and mutation behavior.

## 🗂 Examples

Runnable notebooks and scripts live in [`examples/`](examples):

- [`examples/tebd/`](examples/tebd) — XXZ TEBD vs ED, disordered-XXZ MBL dynamics, scheduler patterns, helper APIs
- [`examples/scarfinder/`](examples/scarfinder) — PXP ScarFinder, XYZ spiral
- [`examples/operator_space/`](examples/operator_space) — TFIM autocorrelators, operator strings, entanglement, custom Hamiltonians, DMT scheduling, XXZ spin-transport exponents, PXP constrained energy transport
- [`examples/tdvp/`](examples/tdvp) — PBC TDVP vs TEBD
- [`examples/chebyshev/`](examples/chebyshev) — energy-cutoff comparison, two-peak spectra, local spectral functions
- [`examples/open_systems/`](examples/open_systems) — boundary-driven XXZ steady state, Pauli–Lindblad TEBD

## 📝 Citation

If MPSToolkit.jl is useful in your research, please cite it:

```bibtex
@misc{MPSToolkit_jl,
  author       = {Ren, Jie},
  title        = {{MPSToolkit.jl}: composable finite matrix-product-state algorithms in Julia},
  year         = {2026},
  howpublished = {\url{https://github.com/jayren3996/MPSToolkit}},
}
```

## 🙏 Acknowledgements

MPSToolkit.jl builds directly on the [ITensor](https://itensor.org/) ecosystem —
[`ITensors.jl`](https://github.com/ITensor/ITensors.jl) and
[`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl) — and on
[`KrylovKit.jl`](https://github.com/Jutho/KrylovKit.jl).

## 📄 License

[MIT](LICENSE) © 2026 Jie Ren and contributors.
