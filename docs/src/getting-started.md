# Getting Started

## Namespaces

The root module is split into feature namespaces:

- `MPSToolkit.Evolution`
- `MPSToolkit.Observables`
- `MPSToolkit.ScarFinder`
- `MPSToolkit.Bases`
- `MPSToolkit.OperatorSpace`
- `MPSToolkit.Models`
- `MPSToolkit.Chebyshev`

This layout is deliberate: evolution, projection, observables, and spectral tools stay composable instead of being fused into one control path.

## Minimal TEBD Example

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
println(expect(psi, "Sz"))
```

## Minimal Operator-Space Example

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 6
sites = pauli_siteinds(nsites)
state = pauli_basis_state(sites, ["Z", "I", "I", "I", "I", "I"])

evolution = tebd_strang_evolution(
  nsites,
  0.05;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.8),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=64,
  cutoff=1e-12,
)

evolve!(state, evolution)
```

## Next Steps

- For the package layout and shared algorithm split, see [Architecture](manual/architecture.md).
- For TEBD / TDVP / DMT concepts, see [Evolution](manual/evolution.md).
- For operator-space workflows, see [Operator Space](manual/operator-space.md).
- For spectral workflows, see [Chebyshev](manual/chebyshev.md).
- For runnable notebooks and scripts, see [Examples](examples.md).
