# XYZ local-operator evolution
#
# What this shows:
# - `spinhalf_xyz_bond_hamiltonian(; Jx, Jy, Jz)` provides a dense two-site XYZ bond term.
# - `Jx`, `Jy`, and `Jz` are the `SxSx`, `SySy`, and `SzSz` couplings.
# - Swapping the model helper changes the physics while keeping the operator-space workflow
#   identical to the TFIM scripts.

using ITensors
using ITensorMPS
using MPSToolkit

nsites = 6
sites = pauli_siteinds(nsites)
# `Jx`, `Jy`, and `Jz` are the XYZ couplings multiplying the `SxSx`, `SySy`, and `SzSz`
# bond terms returned by `spinhalf_xyz_bond_hamiltonian`.
# `maxdim` and `cutoff` remain the TEBD compression controls.
evolution = tebd_strang_evolution(
  nsites,
  0.04;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=0.8, Jz=0.4),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=64,
  cutoff=1e-12,
)

initial = pauli_basis_state(sites, ["I", "I", "Z", "I", "I", "I"])
operator_state = copy(initial)
nsteps = 8
autocorrelation = ComplexF64[inner(initial, operator_state)]

for _ in 1:nsteps
  evolve!(operator_state, evolution)
  push!(autocorrelation, inner(initial, operator_state))
end

println("XYZ local-Z autocorrelation = ", autocorrelation)
