"""
    energy_density(psi, op; kwargs...)

Backend-dispatched local energy-density estimator for a dense local operator.
"""
function energy_density(psi, op; kwargs...)
  throw(MethodError(energy_density, (psi, op)))
end

"""
    bond_entropy(psi, bond)

Return the entanglement entropy associated with a backend-specific bond choice.
"""
function bond_entropy(psi, bond)
  throw(MethodError(bond_entropy, (psi, bond)))
end

"""
    entanglement_spectrum(psi, bond)

Return the normalized Schmidt probabilities associated with a backend-specific bond choice.
"""
function entanglement_spectrum(psi, bond)
  throw(MethodError(entanglement_spectrum, (psi, bond)))
end
