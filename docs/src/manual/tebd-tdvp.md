# TEBD And TDVP

`MPSToolkit.jl` keeps physical-state evolution explicit. TEBD is represented as scheduled local-gate application, while TDVP is represented as MPO-driven variational evolution. Both are first-class public APIs rather than hidden internals.

## TEBD

The TEBD layer is designed around three ideas:

1. a local dense gate or gate provider
2. an explicit schedule describing where those gates act
3. a truncation budget applied at each scheduled step

For many workflows, the helper constructors are the easiest entry point:

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 6
sites = siteinds("S=1/2", nsites)
psi = MPS(sites, n -> isodd(n) ? "Up" : "Dn")

evolution = tebd_strang_evolution(
  nsites,
  0.05;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.8),
  maxdim=32,
  cutoff=1e-12,
)

evolve!(psi, evolution)
```

If you need more control, you can build `LocalGateEvolution` manually or generate gates with `local_gates_from_hamiltonians` and `tebd_evolution_from_hamiltonians`. The scheduler notebooks under [Examples](../examples.md) show standard sweeps, brick schedules, and mixed-span shallow circuits.

## TDVP

`TDVPEvolution` is the MPO-based evolution path. It is a good fit when the generator is naturally represented as an MPO and you want variational time stepping instead of explicit dense local gates.

```julia
using MPSToolkit

evolution = TDVPEvolution(
  generator_mpo,
  -0.1im;
  time_step=-0.1im,
  nsteps=1,
  normalize=true,
  solver_kwargs=(; maxdim=64, cutoff=1e-10),
)

evolve!(psi, evolution)
```

## Choosing Step Counts

For plain TEBD or TDVP usage, the evolution objects keep their own constructor defaults. In other words:

- `LocalGateEvolution(...; nstep=1)` still means one complete schedule pass
- `DMTGateEvolution(...; nstep=1)` still means one forward-and-reverse DMT sweep
- `TDVPEvolution(...; nsteps=1)` still means one TDVP step

ScarFinder is the only place with extra policy on top of those defaults. If you intend to pass one of these evolution objects into `scarfinder_step!`, `scarfinder!`, or `trajectory_refine!`, prefer setting:

- `nstep=10` for TEBD- or DMT-style evolutions
- `nsteps=10` for TDVP evolutions

If you pass a single-step evolution object into ScarFinder, ScarFinder will warn and internally upgrade it to `10` for that ScarFinder call only.

## API

```@docs
LocalGateEvolution
local_gates_from_hamiltonians
tebd_evolution_from_hamiltonians
tebd_strang_schedule
tebd_strang_evolution
tebd_evolve!
TDVPEvolution
tdvp_evolve!
```

## References

- G. Vidal, [Efficient Simulation of One-Dimensional Quantum Many-Body Systems](https://arxiv.org/abs/quant-ph/0310089)
- Steven R. White and Adrian E. Feiguin, [Real-Time Evolution Using the Density Matrix Renormalization Group](https://arxiv.org/abs/cond-mat/0403313)
- Jutho Haegeman, J. Ignacio Cirac, Tobias J. Osborne, Henri Verschelde, and Frank Verstraete, [Time-Dependent Variational Principle for Quantum Lattices](https://arxiv.org/abs/1103.0936)
