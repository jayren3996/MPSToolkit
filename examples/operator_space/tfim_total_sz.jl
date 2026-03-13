# TFIM extensive-operator evolution
#
# What this shows:
# - `pauli_total_sz_state(sites)` builds the Pauli-basis `MPS` for `∑_j S_j^z`.
# - `spinhalf_tfim_bond_hamiltonian(nsites, bond; J, g)` returns the dense two-site TFIM
#   term with open-boundary field splitting built in.
# - `J` controls the Ising `Sz Sz` coupling and `g` the transverse-field `Sx` strength.
# - `measure_bond` selects the bond where the explicit loop below reports operator entanglement.

using ITensors
using ITensorMPS
using MPSToolkit

nsites = 6
sites = pauli_siteinds(nsites)
# `J` is the nearest-neighbor Ising coupling and `g` is the transverse field.
# `maxdim` and `cutoff` control TEBD compression during the evolution.
evolution = tebd_strang_evolution(
  nsites,
  0.05;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.8),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=64,
  cutoff=1e-12,
)

# `pauli_total_sz_state(sites)` builds the extensive initial operator.
initial = pauli_total_sz_state(sites)
operator_state = copy(initial)
measure_bond = 3
nsteps = 8

autocorrelation = ComplexF64[inner(initial, operator_state)]
operator_entanglement = Float64[bond_entropy(operator_state, measure_bond)]

for _ in 1:nsteps
  evolve!(operator_state, evolution)
  push!(autocorrelation, inner(initial, operator_state))
  push!(operator_entanglement, bond_entropy(operator_state, measure_bond))
end

println("final autocorrelation = ", autocorrelation[end])
println("entanglement trace = ", operator_entanglement)
