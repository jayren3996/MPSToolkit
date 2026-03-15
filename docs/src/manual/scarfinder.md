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

Conceptually, the PXP use case looks like:

1. choose a low-entanglement initial product-state ansatz
2. evolve it for a short time under the constrained Hamiltonian
3. explicitly project back to the chosen variational manifold
4. iterate until the trajectory converges toward a stable revival orbit

That makes PXP a good mental model for what ScarFinder is for: not just finding low-entanglement states, but isolating coherent trajectories embedded inside otherwise thermalizing many-body dynamics.

## Example

The current repository example is the periodic XYZ spiral benchmark:

- [examples/scarfinder/xyz_spiral.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/scarfinder/xyz_spiral.jl)

It uses the same core projected-evolution loop as the PXP application, but on a model where the target scarred trajectory is known analytically. That makes it a good first benchmark before moving on to the less structured PXP setting discussed in the paper.

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
