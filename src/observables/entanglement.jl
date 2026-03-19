"""
    _diagonal_values(tensor)

Extract the diagonal entries of a two-index diagonal `ITensor`.

# Arguments
- `tensor`: Diagonal singular-value tensor returned by `svd`.

# Returns
- A dense Julia vector containing the diagonal entries in index order.
"""
function _diagonal_values(tensor::ITensor)
  tensor_inds = inds(tensor)
  length(tensor_inds) == 2 || throw(ArgumentError("expected a two-index diagonal tensor"))
  n = min(dim(tensor_inds[1]), dim(tensor_inds[2]))
  return [tensor[tensor_inds[1] => i, tensor_inds[2] => i] for i in 1:n]
end

"""
    _entropy_from_values(values)

Compute the von Neumann entropy from Schmidt values or equivalent amplitudes.

# Arguments
- `values`: Schmidt coefficients, amplitudes, or any vector whose squared magnitudes should
  be interpreted as probabilities.

# Returns
- The Shannon/von-Neumann entropy `-∑ p log(p)` after normalization.

# Notes
- Zero-probability entries are skipped safely.
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

Return the bond entropy for a finite `MPS`.

# Arguments
- `psi`: Finite matrix-product state.
- `bond`: Bond index at which to cut the chain. If `nothing`, the middle bond is used.

# Returns
- The von Neumann entropy associated with the Schmidt values across the selected cut.

# Notes
- This method delegates spectrum extraction to [`entanglement_spectrum`](@ref) and then
  converts the normalized Schmidt probabilities back into amplitudes before evaluating the
  entropy.
"""
function bond_entropy(psi::MPS, bond::Union{Nothing,Int})
  spectrum = entanglement_spectrum(psi, bond)
  return _entropy_from_values(sqrt.(spectrum))
end

"""
    entanglement_spectrum(psi::MPS, bond)

Return the normalized Schmidt probabilities for a finite `MPS`.

# Arguments
- `psi`: Finite matrix-product state.
- `bond`: Bond index at which to cut. `nothing` selects the central bond.

# Returns
- A vector of normalized Schmidt probabilities.

# Notes
- Out-of-range bonds return an empty vector rather than throwing.
- The state is copied before orthogonalization, so the input `psi` is not mutated.
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

Return a selector-friendly fidelity distance.

# Arguments
- `psi`: Candidate state.
- `reference_state`: Reference state against which fidelity should be measured.

# Returns
- `1 - |⟨reference_state|psi⟩|`, so lower values correspond to higher fidelity.

# Notes
- This helper is intentionally phrased as a distance because ScarFinder selectors minimize
  their score.
"""
function fidelity_distance(psi, reference_state)
  isnothing(reference_state) && throw(ArgumentError("fidelity selection requires a reference state"))
  return 1 - abs(inner(psi, reference_state))
end
