"""
    _diagonal_values(tensor)

Extract the diagonal entries of a two-index diagonal `ITensor`.
"""
function _diagonal_values(tensor::ITensor)
  tensor_inds = inds(tensor)
  length(tensor_inds) == 2 || throw(ArgumentError("expected a two-index diagonal tensor"))
  n = min(dim(tensor_inds[1]), dim(tensor_inds[2]))
  return [tensor[tensor_inds[1] => i, tensor_inds[2] => i] for i in 1:n]
end

"""
    _entropy_from_values(values)

Compute the von Neumann entropy from a vector of Schmidt values or equivalent amplitudes.
"""
function _entropy_from_values(values)
  probabilities = abs2.(values)
  total = sum(probabilities)
  iszero(total) && return 0.0
  normalized = probabilities ./ total
  return -sum(p -> iszero(p) ? 0.0 : p * log(p), normalized)
end

"""
    bond_entropy(psi::MPS, bond)

Return the bond entropy for a finite MPS by orthogonalizing around the selected bond and extracting singular values.
"""
function bond_entropy(psi::MPS, bond::Union{Nothing,Int})
  spectrum = entanglement_spectrum(psi, bond)
  return _entropy_from_values(sqrt.(spectrum))
end

"""
    entanglement_spectrum(psi::MPS, bond)

Return the Schmidt probabilities for a finite MPS across the selected bond.
"""
function entanglement_spectrum(psi::MPS, bond::Union{Nothing,Int})
  bond_index = isnothing(bond) ? max(1, length(psi) ÷ 2) : bond
  (bond_index < 1 || bond_index >= length(psi)) && return Float64[]

  centered = orthogonalize(copy(psi), bond_index)
  left_inds = uniqueinds(centered[bond_index], centered[bond_index + 1])
  _, singular_values, _ = svd(centered[bond_index], left_inds)
  values = abs2.(_diagonal_values(singular_values))
  total = sum(values)
  iszero(total) && return values
  return values ./ total
end

"""
    fidelity_distance(psi, reference_state)

Return a selector-friendly fidelity distance, with lower values corresponding to higher fidelity.
"""
function fidelity_distance(psi, reference_state)
  isnothing(reference_state) && throw(ArgumentError("fidelity selection requires a reference state"))
  return 1 - abs(inner(psi, reference_state))
end
