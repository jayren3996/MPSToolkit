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
