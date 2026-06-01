# TDVP

*MPO-driven variational time evolution on the matrix-product-state manifold.*

## Background

The time-dependent variational principle (TDVP) integrates the Schrödinger equation directly on the manifold of matrix-product states of fixed bond dimension, rather than discretizing the propagator with local gates. This makes it the natural evolution method whenever the generator is supplied as a matrix-product operator (MPO), including long-range and otherwise non-nearest-neighbour Hamiltonians that TEBD cannot factor into a product of two-site gates.

### The variational principle on the MPS manifold

Fix a bond dimension and let ``\mathcal{M}`` be the set of all MPS that can be written with bonds no larger than that budget. ``\mathcal{M}`` is not a linear subspace of the full Hilbert space: a sum of two MPS of bond dimension ``D`` generally has bond dimension ``2D``. It is instead a smooth (nonlinear) manifold embedded in Hilbert space. Exact time evolution ``|\psi(t)\rangle = e^{-iHt}|\psi(0)\rangle`` almost always leaves this manifold, because applying ``H`` grows entanglement and hence the bond dimension.

TDVP asks for the best approximation *that stays on the manifold*. At each instant it chooses the velocity ``|\dot\psi\rangle`` tangent to ``\mathcal{M}`` that is as close as possible to the exact velocity ``-iH|\psi\rangle``. Formally this is the Dirac–Frenkel variational principle,

```math
\big\| \, |\dot\psi\rangle + i H |\psi\rangle \, \big\| \;=\; \min ,
\qquad |\dot\psi\rangle \in T_{|\psi\rangle}\mathcal{M},
```

where ``T_{|\psi\rangle}\mathcal{M}`` is the tangent space to the manifold at the current state.

### The tangent-space projector

The minimization above is an orthogonal projection. Let ``\hat P_{T_{|\psi\rangle}}`` be the orthogonal projector onto the tangent space ``T_{|\psi\rangle}\mathcal{M}``. The TDVP equation of motion is

```math
i \, \frac{d}{dt} |\psi\rangle \;=\; \hat P_{T_{|\psi\rangle}} \, H \, |\psi\rangle .
```

This is exact time evolution with the Hamiltonian replaced by its tangent-space projection. For MPS the tangent-space projector decomposes into a sum of local terms, one per site and per bond, built from the left- and right-orthogonal environments of the canonical form. Crucially these local terms can be applied one at a time by sweeping through the chain, and each one generates an effective *local* Schrödinger equation that is integrated by exponentiating a small effective Hamiltonian (a Krylov/Lanczos exponential of an MPO-contracted local operator). That sweep-and-exponentiate structure is exactly what `tdvp` does internally.

### One-site versus two-site TDVP

There are two standard variants, distinguished by how much of the wavefunction is updated at each local step.

- **One-site TDVP** updates a single MPS tensor at a time. The tangent-space projection is taken at *fixed* bond dimension, so the bond dimensions never change during the sweep. This makes one-site TDVP strictly norm- and (for time-independent ``H``) energy-conserving, very cheap, and completely stable, at the cost of being unable to grow entanglement: if you start with a low-bond-dimension state it stays low, even when the true dynamics demands more. This is the right tool once the bond dimension has saturated, or for a state that is genuinely well-described at fixed ``D``.

- **Two-site TDVP** updates two neighbouring tensors jointly, then re-splits them with a truncated SVD. Because the joint two-site object lives in a larger space, the SVD can *increase* the shared bond up to the truncation budget. Two-site TDVP can therefore grow entanglement adaptively, which is essential early in a quench when bonds are still small, but it reintroduces a truncation (projection onto the larger-but-still-finite manifold) and a slightly higher cost per step. A common strategy is to start two-site to let bonds grow and then switch to one-site once the budget is reached.

### Error and conservation properties

TDVP has no Trotter error: the time integrator is not built from a product-formula splitting of ``H``, so there is no error term from non-commuting Hamiltonian pieces. The only approximation is the **projection error** — the part of ``-iH|\psi\rangle`` that points off the manifold and is discarded by ``\hat P_{T_{|\psi\rangle}}``. Its size is set entirely by how well the chosen bond dimension captures the true state. When the bond dimension is large enough to represent the exact dynamics, TDVP is exact (up to the local exponential's Krylov tolerance and finite step size).

For a time-independent Hamiltonian, the projected flow conserves energy ``\langle\psi|H|\psi\rangle`` exactly in one-site TDVP, because the tangent-space projector is Hermitian and contains ``H|\psi\rangle``'s tangential part — the energy's time derivative is a real expectation of a commutator that vanishes under projection. Two-site TDVP conserves energy up to the truncation performed in the SVD re-split. The norm is preserved by the local unitary exponentials; the optional renormalization step exists mainly to absorb the tiny drift from truncation and from imaginary-time runs.

### When TDVP beats TEBD

TEBD is excellent for short-range Hamiltonians that split cleanly into a small set of two-site gates, and it is conceptually transparent. TDVP is preferable when:

- the Hamiltonian is **long-range or otherwise not a sum of nearest-neighbour two-site terms** — anything you would naturally write as an `OpSum`/`MPO` (power-law couplings, three-site terms, dressed constraints) is handled with no extra work, whereas TEBD would require swap networks or large fused gates;
- you want **no Trotter error**, trading it for a controllable projection error governed solely by the bond dimension;
- you are already carrying an MPO for the same Hamiltonian (e.g. for energy measurements or DMRG), so reusing it as the TDVP generator is free.

Conversely, TEBD is often simpler and faster for strictly local models with small bonds, and two-site TDVP shares TEBD's need for a truncation budget while exact-bond one-site TDVP cannot grow entanglement on its own.

## TDVP in MPSToolkit

`MPSToolkit.jl` exposes TDVP as a first-class evolution backend through the [`TDVPEvolution`](@ref) configuration object and the [`tdvp_evolve!`](@ref) driver. The driver is a thin, in-place wrapper around ITensorMPS's `tdvp`: the configuration object stores the generator MPO and all solver settings, and `tdvp_evolve!` forwards them to `tdvp` and copies the returned state back into the input MPS so that the package-wide in-place `evolve!` convention is preserved.

### `TDVPEvolution`

The constructor takes the generator and the total evolution interval as **positional** arguments, with everything else keyword:

```julia
TDVPEvolution(generator, t;
    time_step=nothing,
    nsteps=nothing,
    nsweeps=nothing,
    reverse_step=true,
    updater_backend="exponentiate",
    updater=nothing,
    normalize=false,
    solver_kwargs=(;),
    schedule=nothing,
)
```

- `generator`: the MPO (or any object `tdvp` accepts) that generates the evolution.
- `t`: the total evolution interval, passed positionally to `tdvp`. The underlying `tdvp` evolves by ``e^{t\,\cdot\,\text{generator}}`` applied to the MPO, so **real-time** Schrödinger evolution under a Hamiltonian ``H`` uses an *imaginary* `t` of the form `-im*T` (the generator is then ``H``), and **imaginary-time** evolution / cooling uses a *real, negative* `t`.
- `time_step`: the internal TDVP time step. When left at `nothing`, it is forwarded as `time_step=nothing` and `tdvp` derives the stepping from `t` and the step count.
- `nsteps`: preferred number of internal TDVP steps. This is the effective step count when it is set.
- `nsweeps`: legacy fallback step-count field, used only when `nsteps` is `nothing`.
- `reverse_step`: whether `tdvp` alternates left-to-right and right-to-left sweep directions (a symmetric, second-order sweep when `true`).
- `updater_backend`: the local-exponential backend forwarded to `tdvp` (default `"exponentiate"`, the Krylov exponential).
- `updater`: an optional custom local updater callback; when `nothing`, `tdvp`'s default updater is used.
- `normalize`: whether `tdvp` renormalizes the state. Defaults to `false`, which preserves the exact (unit) norm of real-time evolution; set it `true` for imaginary-time runs where the norm decays.
- `solver_kwargs`: a `NamedTuple` of extra keyword arguments forwarded **last** to `tdvp` (e.g. `maxdim`, `cutoff`, `nsite`). It must not contain any of the reserved keys `time_step`, `nsteps`, `reverse_step`, `updater_backend`, `updater`, `normalize` — passing one throws an `ArgumentError`, because being forwarded last it would silently override the dedicated keyword. Set those through the dedicated keywords instead.
- `schedule`: optional metadata kept only for interface symmetry with the TEBD-style evolution objects; TDVP does not use a bond schedule.

The constructor validates that `nsteps`/`nsweeps`, when given, are at least `1`, and enforces the reserved-key rule on `solver_kwargs`.

### `tdvp_evolve!` and `evolve!`

[`tdvp_evolve!`](@ref)`(psi, evo)` performs one TDVP evolution call in place:

1. It resolves the effective step count as `evo.nsteps` if present, otherwise `evo.nsweeps`; when both are `nothing`, the step count is left to `tdvp`, which derives it from `evo.t` and `evo.time_step`.
2. It assembles the keyword arguments (`time_step`, `nsteps`, `reverse_step`, `updater_backend`, `normalize`, then `solver_kwargs...`) and calls `tdvp(evo.generator, evo.t, psi; ...)`, threading the optional `updater` when one is set.
3. Because `tdvp` returns a *new* MPS, the wrapper copies the result back into the original storage (`psi[:] = evolved`) so that callers keep the same object.

The generic [`evolve!`](@ref) dispatches `evolve!(psi::MPS, evo::TDVPEvolution)` straight to `tdvp_evolve!`, so the same `evolve!` entry point used for TEBD and DMT also drives TDVP. This is what lets a `TDVPEvolution` be handed to higher-level loops such as `scarfinder_step!`.

### The MPO-generator requirement

Unlike the dense-gate backends, TDVP needs its generator as an MPO. In practice you build one with ITensorMPS's `OpSum`/`MPO` machinery on a set of physical indices (`siteinds`), exactly as you would for a DMRG Hamiltonian. The worked example below constructs such an MPO from scratch so it is fully self-contained.

## Worked example

The following builds a transverse-field Ising MPO on a small open chain, prepares a polarized initial state, and runs real-time TDVP. Everything used here is real API; there are no undefined variables.

```julia
using MPSToolkit
using ITensors
using ITensorMPS

# 1. Physical indices and the Hamiltonian MPO (open-boundary TFIM).
N = 16
sites = siteinds("S=1/2", N)

J = 1.0   # Ising coupling
g = 0.8   # transverse field

os = OpSum()
for j in 1:(N - 1)
    os += -J, "Sz", j, "Sz", j + 1       # -J Sz_j Sz_{j+1}
end
for j in 1:N
    os += -g, "Sx", j                     # -g Sx_j
end
H = MPO(os, sites)

# 2. Initial state: all spins up along z (a domain-wall-free product state).
psi = MPS(sites, n -> "Up")

# 3. Real-time TDVP. The generator is H, so a Schrodinger step over total
#    time T uses t = -im*T. The default two-site sweep lets the bond grow.
T_total = 1.0
nsteps = 20                     # 20 internal steps -> internal dt = T_total/20
evolution = TDVPEvolution(
    H, -im * T_total;
    nsteps=nsteps,
    normalize=false,            # real-time evolution preserves the norm
    solver_kwargs=(maxdim=64, cutoff=1e-12),
)

# 4. Evolve in place and inspect observables.
energy0 = real(inner(psi', H, psi))
evolve!(psi, evolution)          # equivalently: tdvp_evolve!(psi, evolution)
energy1 = real(inner(psi', H, psi))

@show energy0 energy1            # energy is conserved up to truncation
@show maxlinkdim(psi)            # bond dimension after two-site growth
@show abs(norm(psi) - 1)         # norm stays ~ 1 in real time
```

Expected behaviour: starting from the polarized product state (bond dimension 1), the default two-site sweep grows the bond dimension as entanglement builds, up to the `maxdim=64` budget. The energy ``\langle\psi|H|\psi\rangle`` is conserved up to the SVD truncation in the two-site re-split (it would be conserved to Krylov precision in one-site TDVP), and the norm stays essentially `1` because real-time evolution is unitary and `normalize=false`. Increasing `T_total` or `N` makes the dynamics richer; selecting a one-site sweep via `tdvp`'s `nsite` option freezes the bond dimension at its current value.

## Tips and pitfalls

- **Time argument and direction.** Remember that `t` is the *total* interval and that `tdvp` exponentiates `t * generator`. With the Hamiltonian as the generator, real-time evolution needs `t = -im*T`; imaginary-time cooling needs a real, negative `t` (e.g. `t = -beta`). Mixing these up silently turns a quench into a cooling run or vice versa.

- **`time_step` versus `nsteps`.** You can either set `nsteps` (the internal step count, with the internal step size derived as `t / nsteps`) or set `time_step` directly and let the step count follow. Leaving both `nsteps`/`nsweeps` and `time_step` at `nothing` hands the choice to `tdvp`. Pick a step small enough that the local Krylov exponential is well-conditioned; for real time, a smaller step reduces the projection error accumulated per sweep.

- **One-site vs two-site, `maxdim`, `cutoff`.** Use the two-site sweep (the `tdvp` default) with a sensible `maxdim`/`cutoff` while entanglement is still growing, then switch to a one-site sweep via `tdvp`'s `nsite` option (forwarded through `solver_kwargs`) once the bond saturates to lock the cost and make the evolution strictly conserving. `maxdim` and `cutoff` only bite in the two-site re-split; under one-site TDVP the bond dimension is fixed and these have no effect.

- **Real vs imaginary time and `normalize`.** Keep `normalize=false` for real-time evolution so you can use the norm as a truncation-error diagnostic. For imaginary-time evolution the norm decays by construction, so set `normalize=true` (or renormalize yourself) to track the evolving state.

- **Energy conservation as a check.** Because TDVP has no Trotter error, a drift in ``\langle\psi|H|\psi\rangle`` during a real-time run is a direct readout of truncation/projection error. Watching it is the cheapest way to tell whether your bond dimension is adequate.

- **ScarFinder step count.** When you intend to pass a `TDVPEvolution` into `scarfinder_step!`, `scarfinder!`, or `trajectory_refine!`, prefer setting `nsteps=10` explicitly. ScarFinder treats an effective step count of `1` as a bad main-loop setting: it will emit a warning and internally use `10` for that call only. See [ScarFinder](scarfinder.md) for the full rule.

## API

```@docs
TDVPEvolution
tdvp_evolve!
```

## Examples

- [examples/tdvp/pbc_tdvp_vs_tebd.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/tdvp/pbc_tdvp_vs_tebd.ipynb)

## References

- Jutho Haegeman, J. Ignacio Cirac, Tobias J. Osborne, Henri Verschelde, and Frank Verstraete, [Time-Dependent Variational Principle for Quantum Lattices](https://arxiv.org/abs/1103.0936)
- Jutho Haegeman, Christian Lubich, Ivan Oseledets, Bart Vandereycken, and Frank Verstraete, [Unifying time evolution and optimization with matrix product states](https://arxiv.org/abs/1408.5056)
- Sebastian Paeckel, Thomas Köhler, Andreas Swoboda, Salvatore R. Manmana, Ulrich Schollwöck, and Claudius Hubig, [Time-evolution methods for matrix-product states](https://arxiv.org/abs/1901.05824)
