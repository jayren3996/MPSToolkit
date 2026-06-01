# TEBD

*Time-Evolving Block Decimation: real-time evolution of a finite open-boundary MPS by sweeping dense local gates and truncating the bond dimension after each one.*

## Background

TEBD evolves a matrix-product state under a Hamiltonian that is a sum of local terms,

```math
H = \sum_j h_{j,j+1},
```

where each ``h_{j,j+1}`` acts only on a nearest-neighbor pair of sites. The exact propagator ``\exp(-\mathrm{i} H t)`` couples every site at once and cannot be applied to an MPS directly. The idea behind TEBD is to approximate it by a product of small, *local* propagators that each act on a single bond, and to apply those one at a time.

### Trotter–Suzuki decomposition

Split the Hamiltonian into a part supported on odd bonds and a part supported on even bonds,

```math
H = H_\mathrm{odd} + H_\mathrm{even}, \qquad
H_\mathrm{odd} = \sum_{j\ \mathrm{odd}} h_{j,j+1}, \quad
H_\mathrm{even} = \sum_{j\ \mathrm{even}} h_{j,j+1}.
```

All terms inside ``H_\mathrm{odd}`` act on disjoint pairs of sites, so they commute with one another; the same holds for ``H_\mathrm{even}``. The two groups do *not* commute with each other, and that non-commutativity is where the algorithm's systematic error comes from. For a single time step ``\mathrm{d}t`` the first-order Trotter–Suzuki formula reads

```math
\exp\!\big(-\mathrm{i}(H_\mathrm{odd}+H_\mathrm{even})\,\mathrm{d}t\big)
= \exp(-\mathrm{i} H_\mathrm{odd}\,\mathrm{d}t)\,
  \exp(-\mathrm{i} H_\mathrm{even}\,\mathrm{d}t)
  + \mathcal{O}(\mathrm{d}t^2).
```

Because ``H_\mathrm{odd}`` and ``H_\mathrm{even}`` are each sums of commuting two-site terms, their exponentials factorize exactly into a product of independent two-site gates ``\exp(-\mathrm{i}\,h_{j,j+1}\,\mathrm{d}t)``. Each such gate is a dense ``d^2 \times d^2`` matrix (``d`` the local site dimension) that can be applied to an MPS one bond at a time.

### Strang (second-order) splitting

The first-order formula has a local error of order ``\mathrm{d}t^2`` per step, hence a global error of order ``\mathrm{d}t`` after a fixed total time. The symmetric **Strang splitting** halves the odd layer and brackets the even layer between the two halves,

```math
\exp\!\big(-\mathrm{i} H\,\mathrm{d}t\big)
= \exp\!\big(-\tfrac{\mathrm{i}}{2} H_\mathrm{odd}\,\mathrm{d}t\big)\,
  \exp\!\big(-\mathrm{i} H_\mathrm{even}\,\mathrm{d}t\big)\,
  \exp\!\big(-\tfrac{\mathrm{i}}{2} H_\mathrm{odd}\,\mathrm{d}t\big)
  + \mathcal{O}(\mathrm{d}t^3).
```

The odd-numbered terms in the expansion cancel by the symmetry of the bracketing, leaving a local error of order ``\mathrm{d}t^3`` and a global error of order ``\mathrm{d}t^2``. This is the default decomposition used by MPSToolkit's TEBD helpers, and the prefactors ``\tfrac{1}{2}`` (odd, twice) and ``1`` (even, once) are exactly what the Strang schedule encodes.

### The brick-wall schedule

Reading the Strang formula right to left as a sequence of gate applications produces the familiar **brick-wall** circuit: a half step on every odd bond, then a full step on every even bond, then a second half step on every odd bond. One step of the algorithm is therefore one *sweep* over the bond list

```math
(\underbrace{1, 3, 5, \dots}_{\text{odd}},\ \underbrace{2, 4, 6, \dots}_{\text{even}},\ \underbrace{1, 3, 5, \dots}_{\text{odd}}),
```

with the half/full weighting applied to the Hamiltonian before exponentiation. On a finite open chain the boundary bonds simply drop out of the lists where they do not exist; no periodic wraparound is involved.

### Truncation and the bond-dimension budget

Applying a two-site gate fuses two MPS tensors, acts with the gate, and then splits the result back into two tensors with a singular value decomposition. That SVD generically *increases* the bond dimension across the cut. To keep the simulation tractable the new bond is truncated: singular values below `cutoff` are discarded, and at most `maxdim` of them are retained. These two knobs define the bond-dimension budget. Truncating is what makes TEBD efficient, and it is also the second, independent source of error.

### Entanglement growth and the area law

The reason a finite bond dimension can ever be enough is the **area law**: ground states and low-energy states of gapped local Hamiltonians have half-chain entanglement entropy that saturates to a constant, so the Schmidt spectrum decays fast and a modest `maxdim` captures it faithfully. Real-time evolution is the hard case. After a quench, entanglement entropy across a cut typically grows *linearly* in time,

```math
S(t) \sim v\, t,
```

so the bond dimension needed to represent the state grows *exponentially*, ``\chi \sim \mathrm{e}^{S(t)}``. No fixed `maxdim` can keep up indefinitely. For most quenches this entanglement barrier, not CPU time, sets the longest time you can reach, which is why post-quench dynamics are reliable only out to some model-dependent horizon.

### Error sources

Two errors accumulate independently and should be balanced against each other:

- **Trotter error.** Order ``\mathrm{d}t^2`` per step for the Strang scheme, accumulating to ``\mathcal{O}(\mathrm{d}t^2)`` over a fixed total time at fixed `nstep \cdot dt`. Reducing `dt` (and increasing `nstep` to keep the total time fixed) makes it smaller.
- **Truncation error.** Set by the discarded weight at each SVD, controlled by `cutoff` and `maxdim`. It grows as entanglement grows, regardless of how small `dt` is.

Shrinking `dt` below the point where truncation dominates buys nothing, and an enormous `maxdim` is wasted if the Trotter step is coarse. The two should be tightened together.

## TEBD in MPSToolkit

`MPSToolkit.jl` keeps physical-state evolution explicit: TEBD is represented as scheduled local-gate application, exposed as a first-class public API rather than a hidden internal. The layer is built around three ideas:

1. **A local dense gate, or a gate provider.** A gate is just a dense matrix. You can supply one matrix reused at every schedule entry, a vector of matrices indexed by schedule position, or a callable `(bond, index) -> gate`. This is the `gate` field of `LocalGateEvolution`.
2. **An explicit schedule.** A schedule is an ordered list of bond labels saying where each gate is applied and in what order. The Strang brick wall is one such schedule, but it is data you can inspect and modify, not a hidden code path.
3. **A truncation budget.** `maxdim` and `cutoff` are applied at every gate application via the underlying ITensor SVD.

These are bundled in [`LocalGateEvolution`](@ref), constructed as `LocalGateEvolution(gate, dt; schedule, nstep, maxdim, cutoff)`. Here `dt` is the logical time increment associated with *one* `evolve!` call, `schedule` is the ordered bond list, and `nstep` is how many complete passes through `schedule` a single `evolve!` performs.

### Helper constructors

For most workflows you do not build gates by hand. The helpers convert dense local Hamiltonians into gates and assemble the schedule for you:

- [`tebd_strang_schedule`](@ref)`(nsites)` returns the nearest-neighbor odd–even–odd `(schedule, weights)` pair for a finite OBC chain, with the half/full Strang prefactors in `weights`.
- [`local_gates_from_hamiltonians`](@ref)`(hamiltonians, dt)` exponentiates dense local Hamiltonian data into gates via the default map ``h \mapsto \exp(-\mathrm{i}\,\mathrm{d}t\,h)``. It preserves the shape of its input: a matrix in, a matrix out; a vector in, a vector out; a callable in, a callable out.
- [`tebd_evolution_from_hamiltonians`](@ref)`(hamiltonians, dt; schedule, ...)` wraps the previous two into a ready-to-run `LocalGateEvolution` (an explicit `schedule` is required).
- [`tebd_strang_evolution`](@ref)`(nsites, dt; local_hamiltonian, ...)` is the highest-level entry point. You pass a builder `local_hamiltonian(bond, weight) -> h_local`; the helper builds the Strang schedule, evaluates your builder at each scheduled `(bond, weight)`, exponentiates, and returns the `LocalGateEvolution`. The `weight` argument is the Strang prefactor and should multiply the bond Hamiltonian so the half/full bracketing is correct.

### How `evolve!` runs

Running the evolution is a single call to the generic [`evolve!`](@ref)`(psi, evolution)` dispatch, which mutates `psi` in place and returns it. For a `LocalGateEvolution`, one call performs `nstep` complete traversals of `schedule`; at each schedule entry the gate for that `(bond, index)` is resolved and applied with [`tebd_evolve!`](@ref) under the configured `maxdim`/`cutoff`. So the total simulated time of one `evolve!` call is `nstep` Strang sweeps, each of logical duration `dt`.

### Conventions

- **Dense local operators.** Gates are plain dense matrices. Their support is inferred from the chain's local site dimension ``d``: a ``d^2 \times d^2`` matrix is a two-site gate, ``d \times d`` is a one-site gate, and so on. A uniform local dimension across the chain is assumed.
- **Finite open boundaries.** Schedules and gates target a finite OBC chain. The Strang schedule contains no wraparound bond. (`tebd_evolve!` does support a single two-site gate across the last–first boundary when a schedule explicitly uses `bond == length(psi)`, but the standard Strang workflow never does this.)
- **Model builders.** [`spinhalf_tfim_bond_hamiltonian`](@ref) and [`spinhalf_xyz_bond_hamiltonian`](@ref) return dense ``4 \times 4`` two-site Hamiltonians ready to drop into the helpers. The TFIM builder splits its transverse field evenly between neighboring bonds in the bulk and fully on the edges, so that the sum over bond Hamiltonians reproduces the standard open-chain TFIM.

## Worked example

A 12-site transverse-field Ising chain quenched from a Néel product state. The state is evolved for a total time ``t = 10`` as 200 Strang sweeps of ``\mathrm{d}t = 0.05``, and we read off the staggered magnetization and the half-chain entanglement afterwards.

```julia
using MPSToolkit, ITensors, ITensorMPS

# A 12-site transverse-field Ising chain, started in a Néel product state.
nsites = 12
sites  = siteinds("S=1/2", nsites)
psi    = MPS(sites, n -> isodd(n) ? "Up" : "Dn")

# Strang (odd–even–odd) schedule of local gates from the bond Hamiltonian.
# 200 sweeps of dt = 0.05 (total time t = 10).
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

What to expect:

- `expect(psi, "Sz")` (from `ITensorMPS`) returns the per-site ``\langle S^z_j \rangle``. The initial Néel state has these pinned at ``\pm\tfrac{1}{2}`` in a perfect ``+,-,+,-`` pattern. Under the TFIM quench the staggered order is not conserved and relaxes toward zero, so the alternating pattern is washed out.
- [`bond_entropy`](@ref)`(psi, nothing)` returns the von Neumann entropy across the middle cut (passing `nothing` selects the central bond — here the cut between sites 6 and 7). The product initial state is unentangled, ``S = 0``; by ``t = 10`` it has grown to roughly ``1.47`` as the quench spreads correlations across the cut.

With `maxdim = 32` and `cutoff = 1e-12`, the budget is comfortable at these times: a half-chain entropy near ``1.47`` corresponds to an effective Schmidt rank well under 32, so the truncation error stays negligible and the Strang step at ``\mathrm{d}t = 0.05`` is the dominant error. Pushing to much longer times would eventually saturate the budget and demand a larger `maxdim`.

## Building schedules by hand

The helpers are conveniences over public primitives. When you need a non-standard circuit (a plain left-to-right sweep, a custom brick pattern, mixed-span shallow layers, or a position-dependent gate) you can assemble the pieces yourself.

Start from the raw Strang schedule to see what the helper produces:

```julia
using MPSToolkit

nsites = 8
schedule, weights = tebd_strang_schedule(nsites)
# schedule == [1, 3, 5, 7,  2, 4, 6,  1, 3, 5, 7]
# weights  == [0.5, 0.5, 0.5, 0.5,  1.0, 1.0, 1.0,  0.5, 0.5, 0.5, 0.5]
```

Convert a vector of per-entry dense Hamiltonians into per-entry gates, then bundle them with an explicit schedule. Building the gate vector aligned with the schedule lets each entry carry its own Strang weight:

```julia
using MPSToolkit

nsites = 8
dt = 0.05
schedule, weights = tebd_strang_schedule(nsites)

# One dense bond Hamiltonian per schedule entry, weighted for the Strang bracketing.
hamiltonians = [
    weight * spinhalf_xyz_bond_hamiltonian(; Jx = 1.0, Jy = 1.0, Jz = 0.5)
    for (bond, weight) in zip(schedule, weights)
]

# Exponentiate h -> exp(-im * dt * h) entrywise, preserving the vector shape.
gates = local_gates_from_hamiltonians(hamiltonians, dt)

evolution = LocalGateEvolution(
    gates, dt;
    schedule = schedule,
    nstep    = 100,
    maxdim   = 64,
    cutoff   = 1e-10,
)
```

`tebd_evolution_from_hamiltonians` collapses the last two steps — exponentiation and wrapping — into one call, and accepts the same matrix / vector / callable shapes:

```julia
using MPSToolkit

nsites = 8
dt = 0.05
schedule, weights = tebd_strang_schedule(nsites)

hamiltonians = [
    weight * spinhalf_xyz_bond_hamiltonian(; Jx = 1.0, Jy = 1.0, Jz = 0.5)
    for (bond, weight) in zip(schedule, weights)
]

evolution = tebd_evolution_from_hamiltonians(
    hamiltonians, dt;
    schedule = schedule,
    nstep    = 100,
    maxdim   = 64,
    cutoff   = 1e-10,
)
```

For full control over a single gate application — for instance to drive a hand-rolled sweep loop — `tebd_evolve!` applies one dense gate to the block starting at `bond`:

```julia
using MPSToolkit, ITensors, ITensorMPS

nsites = 8
sites  = siteinds("S=1/2", nsites)
psi    = MPS(sites, n -> isodd(n) ? "Up" : "Dn")

gate = local_gates_from_hamiltonians(
    spinhalf_xyz_bond_hamiltonian(; Jx = 1.0, Jy = 1.0, Jz = 0.5),
    0.05,
)

# One left-to-right sweep of the same two-site gate across every bond.
for bond in 1:(nsites - 1)
    tebd_evolve!(psi, gate, bond; maxdim = 64, cutoff = 1e-10)
end
```

A gate provider adds position dependence: passing a callable `(bond, index) -> gate` to `LocalGateEvolution` lets each scheduled application use a different gate, which is how disordered or spatially modulated Hamiltonians are handled.

## Tips and pitfalls

- **Balance `dt` against the truncation budget.** The Strang error is ``\mathcal{O}(\mathrm{d}t^2)`` per unit time; the truncation error is set by `cutoff`/`maxdim`. Shrinking `dt` past the point where truncation dominates wastes work, and a large `maxdim` is wasted on a coarse `dt`. Tighten them together, and converge results in *both* knobs separately.
- **`nstep` versus `dt` set total time.** A single `evolve!` advances by `nstep` Strang sweeps of duration `dt`, i.e. total time `nstep * dt`. To reach a target time more accurately, keep `nstep * dt` fixed while increasing `nstep` (smaller `dt`).
- **Entanglement caps the reachable time.** After a quench the half-chain entropy grows roughly linearly, so the required bond dimension grows exponentially. Watch [`bond_entropy`](@ref): once the effective Schmidt rank approaches `maxdim`, the discarded weight is no longer negligible and results past that time are not trustworthy regardless of how small `dt` is. Raising `maxdim` buys more time, but only logarithmically.
- **Weight the Hamiltonian, not the gate, for Strang.** In `tebd_strang_evolution` and the hand-built schedules, multiply the bond Hamiltonian by `weight` *before* exponentiating. The half/full bracketing is what makes the scheme second order; applying the weight to an already-exponentiated gate is wrong.
- **Match gate support to the local dimension.** A two-site gate on spin-½ sites is ``4 \times 4``. If a gate's size is not a power of the local dimension `tebd_evolve!` raises an error rather than silently misinterpreting the support.
- **ScarFinder wants more than one step.** Plain TEBD keeps its own constructor default `nstep = 1` (one complete schedule pass). ScarFinder, however, expects a meaningful evolution segment between refinement steps. If you intend to pass a `LocalGateEvolution` into `scarfinder_step!`, `scarfinder!`, or `trajectory_refine!`, set `nstep = 10`. Passing a single-step evolution triggers a warning and ScarFinder internally upgrades it to `10` for that call only. See [ScarFinder](scarfinder.md).

## API

```@docs
LocalGateEvolution
local_gates_from_hamiltonians
tebd_evolution_from_hamiltonians
tebd_strang_schedule
tebd_strang_evolution
tebd_evolve!
```

## Examples

- [`examples/tebd/tebd_helper_apis.ipynb`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/tebd/tebd_helper_apis.ipynb) — the helper constructors end to end.
- [`examples/tebd/scheduler_patterns.ipynb`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/tebd/scheduler_patterns.ipynb) — standard sweeps, brick schedules, and mixed-span shallow circuits.
- [`examples/tebd/xxz_tebd_vs_ed.ipynb`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/tebd/xxz_tebd_vs_ed.ipynb) — TEBD benchmarked against exact diagonalization on the XXZ chain.
- [`examples/tebd/disordered_xxz_mbl_tebd.ipynb`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/tebd/disordered_xxz_mbl_tebd.ipynb) — disordered XXZ and many-body localization dynamics via a position-dependent gate provider.

## References

- G. Vidal, [Efficient Simulation of One-Dimensional Quantum Many-Body Systems](https://arxiv.org/abs/quant-ph/0310089) (arXiv:quant-ph/0310089).
- S. R. White and A. E. Feiguin, [Real-Time Evolution Using the Density Matrix Renormalization Group](https://arxiv.org/abs/cond-mat/0403313) (arXiv:cond-mat/0403313).
- S. Paeckel, T. Köhler, A. Swoboda, S. R. Manmana, U. Schollwöck, and C. Hubig, [Time-evolution methods for matrix-product states](https://arxiv.org/abs/1901.05824) (arXiv:1901.05824).
- M. Suzuki, [Generalized Trotter's formula and systematic approximants of exponential operators and inner derivations with applications to many-body problems](https://doi.org/10.1007/BF01609348), Commun. Math. Phys. **51**, 183 (1976).
