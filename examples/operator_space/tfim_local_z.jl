# TFIM local-operator evolution
#
# What this shows:
# - `tebd_strang_evolution(nsites, dt; local_hamiltonian, ...)` builds an odd-even-odd
#   nearest-neighbor TEBD schedule automatically.
# - `dt` is the physical time step for one full Strang sweep.
# - `local_hamiltonian=(bond, weight) -> ...` must return the dense local Hamiltonian for
#   one scheduled update. `bond` is the left site of the two-site term, and `weight` is the
#   Strang prefactor (`0.5` on half steps, `1.0` on full steps).
# - `map_hamiltonian=pauli_gate_from_hamiltonian` converts each dense local Hamiltonian into
#   the corresponding Pauli-basis superoperator gate.
# - `maxdim` and `cutoff` are the compression controls passed through each TEBD update.
# - The measurement loop is written out explicitly below: each step calls `evolve!`, then
#   records `inner(reference, operator_state)` and `bond_entropy(operator_state, measure_bond)`.

using ITensors
using ITensorMPS
using MPSToolkit

nsites = 6
sites = pauli_siteinds(nsites)
# `nsites` is the chain length and `0.05` is the physical time step of one Strang sweep.
# `J=1.0` sets the `Sz Sz` coupling, `g=0.8` sets the transverse field.
# `maxdim` bounds the compressed bond dimension and `cutoff` is the truncation threshold.
evolution = tebd_strang_evolution(
  nsites,
  0.05;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.8),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=64,
  cutoff=1e-12,
)

# `pauli_basis_state(..., ["Z", "I", ...])` initializes a single-site `Z` operator.
initial = pauli_basis_state(sites, ["Z", "I", "I", "I", "I", "I"])
operator_state = copy(initial)
measure_bond = 3
nsteps = 8

times = Float64[0.0]
autocorrelation = ComplexF64[inner(initial, operator_state)]
operator_entanglement = Float64[bond_entropy(operator_state, measure_bond)]

for step in 1:nsteps
  evolve!(operator_state, evolution)
  push!(times, step * evolution.dt)
  push!(autocorrelation, inner(initial, operator_state))
  push!(operator_entanglement, bond_entropy(operator_state, measure_bond))
end

println("times = ", times)
println("autocorrelation = ", autocorrelation)
println("final maxlinkdim = ", maxlinkdim(operator_state))
