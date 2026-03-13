# TFIM Pauli-string evolution
#
# What this shows:
# - `pauli_basis_state(sites, labels)` builds a single Pauli string in operator space.
#   Here `["X", "Y", "Z", "I", "I", "I"]` means `X ⊗ Y ⊗ Z ⊗ I ⊗ I ⊗ I`.
# - The TEBD arguments have the same meaning as in `operator_space_tfim_local_z.jl`:
#   `dt` is the sweep time step, `maxdim` is the bond cap, and `cutoff` is the compression
#   threshold.

using ITensors
using ITensorMPS
using MPSToolkit

nsites = 6
sites = pauli_siteinds(nsites)
# `nsites`, `dt`, `J`, `g`, `maxdim`, and `cutoff` play the same roles as in the local-`Z`
# TFIM example: chain length, sweep time step, model couplings, and compression controls.
evolution = tebd_strang_evolution(
  nsites,
  0.05;
  local_hamiltonian=(bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.8),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=64,
  cutoff=1e-12,
)

initial = pauli_basis_state(sites, ["X", "Y", "Z", "I", "I", "I"])
operator_state = copy(initial)
nsteps = 8
autocorrelation = ComplexF64[inner(initial, operator_state)]

for _ in 1:nsteps
  evolve!(operator_state, evolution)
  push!(autocorrelation, inner(initial, operator_state))
end

println("initial string autocorrelation trace = ", autocorrelation)
