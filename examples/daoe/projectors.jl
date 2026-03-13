# Operator-space projector example
#
# What this shows:
# - `pauli_siteinds(nsites)` creates `nsites` Pauli-basis sites in the default local
#   ordering `(I, X, Y, Z)`.
# - `pauli_daoe_projector(...; lstar, gamma)` preserves strings up to Pauli weight `lstar` and
#   damps higher-weight strings by factors controlled by `gamma`.
# - `fdaoe_projector(...; wstar, gamma)` uses the fermionic / Jordan-Wigner weight cutoff
#   `wstar` instead of the raw Pauli weight.
# - `apply(projector, state; maxdim, cutoff)` applies the MPO and compresses the result:
#   `maxdim` bounds bond growth and `cutoff` drops singular values below the threshold.
#
# The script prints Hilbert-Schmidt overlaps `⟨O, P O⟩` for a few representative strings.

using ITensors
using ITensorMPS
using MPSToolkit

sites = pauli_siteinds(4)
gamma = 0.3

# `lstar` is the preserved Pauli-weight cutoff for DAOE.
# `gamma` sets the damping strength for each unit of weight above `lstar`.
# `wstar` is the fermionic-weight cutoff for FDAOE.
cases = [
  ("DAOE XIII", pauli_daoe_projector(sites; lstar=1, gamma=gamma), ["X", "I", "I", "I"]),
  ("DAOE XYZI", pauli_daoe_projector(sites; lstar=1, gamma=gamma), ["X", "Y", "Z", "I"]),
  ("FDAOE ZIII", fdaoe_projector(sites; wstar=2, gamma=gamma), ["Z", "I", "I", "I"]),
  ("FDAOE XZZX", fdaoe_projector(sites; wstar=2, gamma=gamma), ["X", "Z", "Z", "X"]),
]

for (name, projector, labels) in cases
  state = pauli_basis_state(sites, labels)
  # `maxdim` and `cutoff` control the MPO-application compression used by `apply`.
  value = inner(state, apply(projector, state; maxdim=32, cutoff=0.0))
  println(rpad(name, 14), " => ", value)
end
