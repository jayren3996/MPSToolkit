# ScarFinder

The ScarFinder layer packages a small but flexible projection workflow:

1. evolve the state
2. truncate or project it
3. optionally match a target energy
4. optionally refine along a short projected trajectory

That split makes ScarFinder usable with both local-gate and MPO-based evolution backends.

## Minimal Workflow

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

The important point is that projection is explicit. `scarfinder_step!` does not hide a bespoke evolution backend; it composes the same public `evolve!`, `project!`, and `match_energy!` style building blocks that you can also call yourself.

## PXP Background

The PXP model is the standard kinetically constrained spin-chain model for Rydberg-blockaded atom arrays, and it is one of the best-known settings for quantum many-body scar dynamics. In the ScarFinder paper, this projection loop is used on the PXP model to uncover a trajectory with nearly perfect revivals in the thermodynamic limit without assuming prior knowledge of the scar tower.

For an open chain, the Hamiltonian used in the example is

```math
H_{\mathrm{PXP}} = \sum_{j=2}^{L-1} 2 P^{\downarrow}_{j-1} X_j P^{\downarrow}_{j+1},
```

where `P^\downarrow = |\downarrow\rangle\langle\downarrow|` enforces the blockade constraint on the neighboring sites and `X_j` flips the center spin.

Conceptually, the PXP use case looks like:

1. choose a low-entanglement initial product-state ansatz
2. evolve it for a short time under the constrained Hamiltonian
3. explicitly project back to the chosen variational manifold
4. iterate until the trajectory converges toward a stable revival orbit

That makes PXP a good mental model for what ScarFinder is for: not just finding low-entanglement states, but isolating coherent trajectories embedded inside otherwise thermalizing many-body dynamics.

## Example

- [examples/scarfinder/pxp_scarfinder.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/scarfinder/pxp_scarfinder.ipynb)
  A `L = 32` open-chain PXP notebook using 3-site TEBD gates, `dt = 0.01` diagnostics, `Delta t = 0.1` projected steps, and a 200-step ScarFinder loop.

## Core PXP Setup

The notebook uses the public API directly. The minimal ScarFinder loop is:

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
    nstep=1,
    maxdim=64,
    cutoff=1e-10,
)

diagnostic_evolution = tebd_evolution_from_hamiltonians(
    fill(local_hamiltonian, length(schedule)),
    0.01;
    schedule=schedule,
    nstep=1,
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

The notebook adds warmup, convergence histories, and a short fine-TEBD diagnostic against the `|Z2>` / `|Z2bar>` revival family, but the loop above is the core ScarFinder setup.

## Selectors And Targeting

- `BondDimTruncation` controls the explicit projection budget.
- `EnergyTarget` adds a post-step correction loop toward a desired expectation value.
- `EntropySelector` and `FidelitySelector` provide small scoring rules for trajectory refinement.

This is useful when you want to search for low-entanglement or approximately recurrent trajectories without hard-wiring the search logic into the evolution engine itself.

## API

```@docs
BondDimTruncation
EnergyTarget
SelectionContext
EntropySelector
FidelitySelector
project!
MPSToolkit.ScarFinder.match_energy!
MPSToolkit.ScarFinder.trajectory_refine!
scarfinder_step!
scarfinder!
```

## References

- Jie Ren, Andrew Hallam, Lei Ying, and Zlatko Papic, [ScarFinder: A detector of optimal scar trajectories in quantum many-body dynamics](https://journals.aps.org/prxquantum/accepted/10.1103/8g8w-nkwx)
