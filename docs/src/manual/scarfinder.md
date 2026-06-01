# ScarFinder

*An explicit evolve → project → energy-match → refine loop for locating optimal quantum-many-body-scar trajectories in finite open-boundary MPS dynamics.*

## Background

Generic isolated quantum many-body systems are expected to thermalize. The eigenstate thermalization hypothesis predicts that local observables relax to their thermal values and that an arbitrary initial state loses memory of its preparation under unitary dynamics. **Quantum many-body scars** are a striking exception to this picture: a small, structured set of nonthermal eigenstates sits embedded in an otherwise thermal spectrum, and a special class of initial states couples almost exclusively to that subspace. Such states evade thermalization for long times, producing persistent coherent oscillations of local observables. Because only a measure-zero set of states behaves this way while the bulk of the spectrum thermalizes normally, the phenomenon is called **weak ergodicity breaking** — ergodicity fails on a thin, dynamically isolated submanifold rather than across the whole Hilbert space.

The canonical setting is the **PXP model**, the effective description of a one-dimensional chain of Rydberg atoms in the **Rydberg-blockade** regime. When the blockade forbids two neighboring atoms from being simultaneously excited, the hopping term that flips a single atom survives only if both neighbors are in the ground state. Quench experiments on these arrays, prepared in the period-two Néel (``\mathbb{Z}_2``) product state, showed long-lived revivals of the local density rather than the rapid relaxation expected of a thermalizing chain. The revivals are the dynamical fingerprint of an approximate, equally spaced tower of scar eigenstates, and the ``\mathbb{Z}_2`` state has large overlap with that tower.

A useful way to think about scar dynamics is geometric. The exact trajectory ``|\psi(t)\rangle`` of a scarring initial state stays close to a low-dimensional, low-entanglement **variational manifold**: a family of states (for PXP, essentially the time-dependent-variational-principle manifold of low bond dimension) on which the dynamics nearly closes. The revivals occur because the trajectory traces an almost-periodic orbit confined to that manifold, leaking only slowly into the surrounding thermal bulk. "Finding a scar trajectory" then means identifying the initial condition and the projected orbit on which this leakage is minimized — the trajectory that revives best.

ScarFinder operationalizes that geometric idea. Rather than diagonalizing the Hamiltonian or assuming prior knowledge of the scar tower, it repeatedly evolves a state for a short time under the physical generator and **projects** it back onto a chosen low-entanglement manifold. The projection is realized concretely as bond-dimension truncation of the MPS. Iterating evolve-then-project drives the state toward the trajectory that the manifold can carry coherently, which is precisely the optimal scar orbit. Optional energy targeting fixes the trajectory to a chosen energy density, and an optional refinement pass scans a short projected orbit and keeps the state a selector ranks best.

## The ScarFinder loop

ScarFinder is deliberately not a hidden backend with its own internal state machine. It is a small, composable loop assembled from public building blocks, each of which you can also call on its own:

1. **Evolve.** Advance the state for a short interval under a physical or effective generator using the standard evolution objects (`LocalGateEvolution` / dense TEBD gates, `DMTGateEvolution`, or `TDVPEvolution`) through [`evolve!`](@ref).
2. **Project.** Truncate the evolved state back into the chosen variational manifold with [`project!`](@ref). For a finite `MPS`, `project!(psi, ::BondDimTruncation)` delegates to ITensor's `truncate!`, so the manifold is exactly the set of MPS at or below a target bond dimension. Truncation *is* the projection — there is no separate manifold-fitting routine.
3. **Energy-match (optional).** Nudge the projected state toward a target energy density with [`MPSToolkit.ScarFinder.match_energy!`](@ref), driven by an `EnergyTarget`. This keeps the trajectory pinned to the spectral region where the scar tower lives instead of drifting after each truncation.
4. **Refine (optional).** Scan a few extra projected steps with [`MPSToolkit.ScarFinder.trajectory_refine!`](@ref) and keep the state that a selector (`EntropySelector` or `FidelitySelector`) scores best, ranking lowest-score-first.

[`scarfinder_step!`](@ref) bundles steps 1–3 into one bare update: evolve, project, then optionally match energy. [`scarfinder!`](@ref) runs that step `niter` times and, if `refine=true`, finishes with step 4. The decomposition is the point of the design. Because projection is explicit rather than fused into the evolution engine, you can swap the projector, the energy-correction rule, or the selector independently, and you can reproduce a full ScarFinder step by hand from `evolve!`, `project!`, and `match_energy!`.

```julia
using MPSToolkit

truncation = BondDimTruncation(16; cutoff=1e-10)
target = EnergyTarget(0.0; operator=hamiltonian_mpo, tol=1e-8, alpha=0.05, maxstep=8)
selector = EntropySelector(; bond=3)

scarfinder!(
  psi,
  evolution,
  truncation,
  10;
  target_energy=target,
  refine=true,
  selector=selector,
  refine_steps=6,
)
```

## Step-count rule

ScarFinder treats an effective main-loop evolution step count of `1` as a misconfiguration. A single Trotter pass per projection generally does not advance the trajectory far enough for the evolve-then-project cycle to behave as intended, so the entry points guard against it.

When [`scarfinder_step!`](@ref), [`scarfinder!`](@ref), or [`MPSToolkit.ScarFinder.trajectory_refine!`](@ref) receive an evolution whose effective step count is `1` —

- a `LocalGateEvolution` with `nstep == 1`,
- a `DMTGateEvolution` with `nstep == 1`, or
- a `TDVPEvolution` with effective steps `1` (via `nsteps`, or `nsweeps` when `nsteps` is `nothing`)

— ScarFinder emits a warning and internally rebuilds an equivalent evolution object with step count `10` for that call only. The original object is untouched.

This normalization is local to ScarFinder. The global TEBD, DMT, and TDVP constructors keep their own defaults; only the ScarFinder entry points apply the rule. In `scarfinder!` it is applied once up front, so a long loop does not emit repeated warnings. If you are deliberately building an evolution object for ScarFinder, set `nstep=10` (or `nsteps=10`) explicitly to advance more per projection and silence the warning.

The internal one-step correction moves inside `match_energy!` are *not* affected. Those use `nstep=1` (dense path) or `nsteps=1` (MPO/TDVP path) on purpose: each is a single narrow post-step energy-correction or rollback move, not the main ScarFinder trajectory evolution.

## PXP example

The PXP model is the standard kinetically constrained spin-chain model for Rydberg-blockaded atom arrays, and one of the best-studied settings for scar dynamics. In the ScarFinder paper, this projection loop is run on PXP to uncover a trajectory with nearly perfect ``\mathbb{Z}_2`` revivals in the thermodynamic limit, without assuming prior knowledge of the scar tower.

For an open chain, the Hamiltonian used in the example is

```math
H_{\mathrm{PXP}} = \sum_{j=2}^{L-1} 2\, P^{\downarrow}_{j-1}\, X_j\, P^{\downarrow}_{j+1},
```

where ``P^{\downarrow} = |{\downarrow}\rangle\langle{\downarrow}|`` projects a site onto its ground state and ``X_j`` flips the center spin. The sandwiching projectors encode the blockade: the central flip on site ``j`` acts only when both neighbors ``j-1`` and ``j+1`` are in ``|{\downarrow}\rangle``. The factor of ``2`` follows from writing the central flip as ``X = 2 S^x`` with spin-``1/2`` operators, which matches the `"Sx"` convention used when building the MPO below.

The workflow on PXP follows the four-step loop directly:

1. start from the low-entanglement ``\mathbb{Z}_2`` product state, which has large overlap with the scar tower,
2. evolve it for a short interval under the constrained Hamiltonian,
3. explicitly project back onto the low-bond-dimension manifold, optionally re-pinning the energy,
4. iterate until the trajectory converges to a stable revival orbit.

The shipped notebook uses the public API directly. The Hamiltonian appears in two forms: a dense three-site local gate ``2\,P^{\downarrow} X P^{\downarrow}`` that drives the TEBD evolution, and an equivalent `MPO` that measures the energy density and supplies the operator for energy targeting. The schedule applies the three-site gate on every site triple along the open chain.

```julia
using MPSToolkit
using ITensors
using ITensorMPS
using LinearAlgebra

function pxp_local_hamiltonian()
    projector_dn = ComplexF64[0 0; 0 1]
    pauli_x = ComplexF64[0 1; 1 0]
    return kron(projector_dn, kron(pauli_x, projector_dn))
end

function pxp_schedule(nsites::Int)
    starts = Int[]
    for offset in 1:3
        append!(starts, offset:3:(nsites - 2))
    end
    return starts
end

function pxp_mpo(sites)
    opsum = OpSum()
    for j in 2:(length(sites) - 1)
        opsum += 2.0, "ProjDn", j - 1, "Sx", j, "ProjDn", j + 1
    end
    return MPO(opsum, sites)
end

nsites = 32
sites = siteinds("S=1/2", nsites)
schedule = pxp_schedule(nsites)
local_hamiltonian = pxp_local_hamiltonian()
hamiltonian_mpo = pxp_mpo(sites)

scar_evolution = tebd_evolution_from_hamiltonians(
    fill(local_hamiltonian, length(schedule)),
    0.1;
    schedule=schedule,
    nstep=10,
    maxdim=64,
    cutoff=1e-10,
)

truncation = BondDimTruncation(1; cutoff=1e-10)
energy_target = EnergyTarget(0.0; operator=hamiltonian_mpo, tol=1e-8, alpha=0.1, maxstep=4)
z2 = MPS(sites, n -> isodd(n) ? "Up" : "Dn")
psi = deepcopy(z2)

for _ in 1:200
    scarfinder_step!(psi, scar_evolution, truncation; target_energy=energy_target)
end
```

A few details worth noting. The dense gate is built from `kron(projector_dn, kron(pauli_x, projector_dn))`, so it acts on three consecutive sites, and `pxp_schedule` lays those three-site windows across the chain. The MPO term `2.0, "ProjDn", j-1, "Sx", j, "ProjDn", j+1` mirrors the gate (the custom `"ProjDn"` op is the down-projector ``P^{\downarrow}``), so the energy measured for targeting is consistent with the generator that drives the dynamics. Here the projection budget is `BondDimTruncation(1)`, the most aggressive manifold: each projection compresses the state back to bond dimension one, which is exactly the product-state manifold that the ``\mathbb{Z}_2`` revival orbit is built on. The example passes `nstep=10` explicitly, the recommended ScarFinder setting, so no step-count warning is emitted.

If you drive ScarFinder with TDVP instead of dense-gate TEBD, the same recommendation applies: prefer an explicit `nsteps=10` configuration when the evolution object is intended for ScarFinder.

## Selectors and targeting

The search logic is exposed as small, independently replaceable pieces:

- `BondDimTruncation` sets the explicit projection budget — the `maxdim` (and `cutoff`) that define the variational manifold. Smaller `maxdim` is a stronger projection. This is the object [`project!`](@ref) consumes for finite-`MPS` truncation.
- `EnergyTarget` configures the post-step energy-correction loop used by `match_energy!`. Its fields are the desired `target` energy, the `operator` used to measure it (a dense matrix or an `MPO`), an early-stop `tol` on the residual, a proportional step size `alpha`, and a cap `maxstep` on correction iterations per ScarFinder step. A dense `operator` routes through a local-gate correction path; an `MPO` routes through a one-step TDVP correction path. Both apply a small imaginary-time-like nudge ``\propto \alpha\,(\langle H\rangle - \text{target})`` and roll the move back if it overshoots or fails to improve, so the energy lands near the target rather than oscillating.
- `EntropySelector` scores a candidate by the bipartite entanglement entropy at a chosen `bond` (`EntropySelector(; bond=3)`); with `bond=nothing` the backend default bond is used. Lower entropy scores better, so this selector prefers the least-entangled candidate along the refinement scan — the most scar-like.
- `FidelitySelector` scores a candidate by its fidelity distance to a reference state. The reference is supplied separately through `SelectionContext` rather than stored in the selector itself, so the same selector can be reused with different references.
- `SelectionContext` carries optional auxiliary data shared across selector calls — currently the `reference_state` that `FidelitySelector` requires. It is threaded into refinement through the `selector_context` keyword.

All selectors are ranked **lowest-score-first** by `MPSToolkit.ScarFinder.trajectory_refine!`: on a tie the incumbent (earlier) state is kept, so the refinement scan never displaces an equally good starting state, and `NaN` scores are ignored (they never win) and trigger a warning. The common scoring entry point is `MPSToolkit.score(selector, psi, context)`, which both `match_energy!`-driven workflows and direct callers can use to evaluate a candidate.

Because each piece is a plain value rather than a buried code path, you can replace the projector, the energy-correction rule, or the selector without rewriting the loop.

## API

```@docs
BondDimTruncation
EnergyTarget
SelectionContext
EntropySelector
FidelitySelector
project!
MPSToolkit.score
MPSToolkit.ScarFinder.match_energy!
MPSToolkit.ScarFinder.trajectory_refine!
scarfinder_step!
scarfinder!
```

## Examples

- [examples/scarfinder/pxp_scarfinder.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/scarfinder/pxp_scarfinder.ipynb) — an ``L = 32`` open-chain PXP notebook using three-site TEBD gates, explicit bond-dimension projection, energy targeting, and overlap/entropy/energy diagnostics for the ``\mathbb{Z}_2`` revival orbit.
- [examples/scarfinder/xyz_spiral.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/scarfinder/xyz_spiral.jl) — a spin-``1/2`` XYZ-chain script tracking a projected spiral trajectory.

## References

- Jie Ren, Andrew Hallam, Lei Ying, and Zlatko Papić, *ScarFinder: A detector of optimal scar trajectories in quantum many-body dynamics*, [PRX Quantum (accepted)](https://journals.aps.org/prxquantum/accepted/10.1103/8g8w-nkwx).
- C. J. Turner, A. A. Michailidis, D. A. Abanin, M. Serbyn, and Z. Papić, *Weak ergodicity breaking from quantum many-body scars*, [arXiv:1711.03528](https://arxiv.org/abs/1711.03528).
- M. Serbyn, D. A. Abanin, and Z. Papić, *Quantum many-body scars and weak breaking of ergodicity*, [arXiv:2011.09486](https://arxiv.org/abs/2011.09486).
