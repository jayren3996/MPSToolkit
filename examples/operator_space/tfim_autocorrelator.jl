# TFIM autocorrelator trace
#
# What this shows:
# - The explicit TEBD loop below records the whole autocorrelation trace `⟨O(0), O(t)⟩`
#   rather than only the final operator.
# - `nsteps=10` means ten Strang sweeps are applied.
# - `measure_bond=4` selects the middle bond of the 8-site chain for the accompanying
#   operator-entanglement measurement stored in the returned data.

using ITensors
using ITensorMPS
using MPSToolkit

nsites = 8
sites = pauli_siteinds(nsites)
# `J` is the Ising coupling and `g` the transverse-field strength inside the TFIM bond term.
# `maxdim` and `cutoff` are the TEBD compression controls.
evolution = tebd_strang_evolution(
  nsites,
  0.04;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=1.05),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=96,
  cutoff=1e-12,
)

local_z = pauli_basis_state(sites, ["I", "I", "I", "Z", "I", "I", "I", "I"])
operator_state = copy(local_z)
measure_bond = 4
nsteps = 10
times = Float64[0.0]
autocorrelation = ComplexF64[inner(local_z, operator_state)]
operator_entanglement = Float64[bond_entropy(operator_state, measure_bond)]

for step in 1:nsteps
  evolve!(operator_state, evolution)
  push!(times, step * evolution.dt)
  push!(autocorrelation, inner(local_z, operator_state))
  push!(operator_entanglement, bond_entropy(operator_state, measure_bond))
end

for (time, value) in zip(times, autocorrelation)
  println(time, " ", value)
end
