# TFIM operator-entanglement growth
#
# What this shows:
# - The same explicit TEBD loop can track operator entanglement instead of only the final state.
# - `measure_bond=4` asks for the bipartite operator entanglement across the center bond.

using ITensors
using ITensorMPS
using MPSToolkit

nsites = 8
sites = pauli_siteinds(nsites)
# `J` is the Ising coupling, `g` is the transverse field, and the remaining keywords control
# TEBD compression exactly as in the other TFIM scripts.
evolution = tebd_strang_evolution(
  nsites,
  0.04;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=1.05),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=96,
  cutoff=1e-12,
)

initial = pauli_basis_state(sites, ["I", "I", "I", "Z", "I", "I", "I", "I"])
operator_state = copy(initial)
measure_bond = 4
nsteps = 10
operator_entanglement = Float64[bond_entropy(operator_state, measure_bond)]

for _ in 1:nsteps
  evolve!(operator_state, evolution)
  push!(operator_entanglement, bond_entropy(operator_state, measure_bond))
end

println("operator entanglement = ", operator_entanglement)
