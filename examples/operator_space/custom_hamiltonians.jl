# Custom local-Hamiltonian operator dynamics
#
# What this shows:
# - You can bypass model helpers entirely and provide dense local Hamiltonian matrices
#   directly through `local_hamiltonian=(bond, weight) -> ...`.
# - `spinhalf_matrices()` returns dense `I, Sx, Sy, Sz` matrices in the spin convention
#   used by the package, so you can assemble custom terms without redefining the basis.
# - The example alternates `Sx ⊗ Sx` and `Sz ⊗ Sz` by bond index and still uses the same
#   TEBD and operator-dynamics helpers as the model-specific scripts.

using ITensors
using ITensorMPS
using MPSToolkit

spins = spinhalf_matrices()
nsites = 6
sites = pauli_siteinds(nsites)

# `nsites` is the chain length.
# `0.03` is the time step of one full odd-even-odd Strang sweep.
# `local_hamiltonian` returns the dense local term for one scheduled update.
# `map_hamiltonian` converts each dense local term into the corresponding Pauli-basis superoperator gate.
# `maxdim` is the maximum intermediate bond dimension during TEBD compression.
# `cutoff` discards singular values below the given threshold during compression.
custom_evolution = tebd_strang_evolution(
  nsites,
  0.03;
  # `bond` selects which local term to use; `weight` applies the Strang half-step/full-step factor.
  local_hamiltonian=(bond, weight) -> weight * (isodd(bond) ? kron(spins.Sx, spins.Sx) : kron(spins.Sz, spins.Sz)),
  map_hamiltonian=pauli_gate_from_hamiltonian,
  maxdim=64,
  cutoff=1e-12,
)

initial = pauli_basis_state(sites, ["X", "I", "I", "I", "I", "I"])
operator_state = copy(initial)
nsteps = 8
autocorrelation = ComplexF64[inner(initial, operator_state)]

for _ in 1:nsteps
  evolve!(operator_state, custom_evolution)
  push!(autocorrelation, inner(initial, operator_state))
end

println("custom dynamics autocorrelation = ", autocorrelation)
