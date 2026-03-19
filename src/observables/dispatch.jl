"""
    energy_density(psi, op; kwargs...)

Backend-dispatched energy-density estimator.

# Arguments
- `psi`: State whose energy density should be estimated.
- `op`: Observable representation understood by a concrete backend. In this package that
  is typically either a dense local operator matrix or an `MPO`.

# Keyword Arguments
- `kwargs...`: Backend-specific options forwarded to the concrete method.

# Returns
- A real-valued energy density.

# Notes
- The fallback method throws a `MethodError`.
"""
function energy_density(psi, op; kwargs...)
  throw(MethodError(energy_density, (psi, op)))
end

"""
    bond_entropy(psi, bond)

Backend-dispatched bond entanglement entropy estimator.

# Arguments
- `psi`: State whose bipartite entanglement should be measured.
- `bond`: Backend-specific bond selector. For finite `MPS` states, `nothing` means
  "use the default central bond".

# Returns
- The von Neumann entropy associated with the selected cut.
"""
function bond_entropy(psi, bond)
  throw(MethodError(bond_entropy, (psi, bond)))
end

"""
    entanglement_spectrum(psi, bond)

Backend-dispatched entanglement-spectrum estimator.

# Arguments
- `psi`: State whose Schmidt spectrum should be extracted.
- `bond`: Backend-specific bond selector.

# Returns
- A vector of normalized Schmidt probabilities.
"""
function entanglement_spectrum(psi, bond)
  throw(MethodError(entanglement_spectrum, (psi, bond)))
end
