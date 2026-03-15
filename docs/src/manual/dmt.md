# DMT

Density matrix truncation (DMT) is implemented here as an operator-space truncation backend that shares the same gate-scheduling idea as TEBD while changing the local truncation rule. The scheduler decides where gates act; the DMT kernel decides how the post-gate truncation is performed.

## Minimal Workflow

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 6
sites = pauli_siteinds(nsites)
state = pauli_basis_state(sites, ["Z", "I", "I", "I", "I", "I"])

gate = pauli_gate_from_hamiltonian(spinhalf_tfim_bond_hamiltonian(nsites, 1; J=1.0, g=0.6), 0.05)
schedule = collect(1:(nsites - 1))

evolution = DMTGateEvolution(
  fill(gate, length(schedule)),
  0.05;
  schedule=schedule,
  reverse_schedule=collect(reverse(schedule)),
  maxdim=16,
  cutoff=1e-10,
)

dmt_evolve!(state, evolution)
```

If you want to work at the lowest level, `dmt_step!` applies one local gate and then performs the associated DMT truncation. The scheduled driver simply repeats that step over forward and reverse schedules.

## Relation To TEBD

- `LocalGateEvolution` and `DMTGateEvolution` both represent scheduled local updates.
- `tebd_evolve!` uses ordinary gate application and SVD truncation.
- `dmt_step!` uses gate application followed by DMT-preserving truncation.

This keeps TEBD and DMT on the same conceptual footing instead of treating DMT as an unrelated evolution family.

## API

```@docs
DMTOptions
DMTGateEvolution
dmt_step!
dmt_evolve!
```

## Examples

- [examples/operator_space/dmt_scheduler.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/dmt_scheduler.ipynb)
- [examples/open_systems/pauli_lindblad_tebd.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/open_systems/pauli_lindblad_tebd.ipynb)

## References

- Stuart Yi-Thomas, Brayden Ware, Jay D. Sau, and Christopher David White, [Comparing numerical methods for hydrodynamics in a one-dimensional lattice spin model](https://arxiv.org/abs/2310.06886)
- En-Jui Kuo, Brayden Ware, Peter Lunts, Mohammad Hafezi, and Christopher David White, [Energy diffusion in weakly interacting chains with fermionic dissipation-assisted operator evolution](https://arxiv.org/abs/2311.17148)
