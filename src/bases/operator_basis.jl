const _OPERATOR_BASIS_CACHE = Dict{Int,Vector{Matrix{ComplexF64}}}()

"""
    operator_basis_matrices(d)

Return the normalized generalized Gell-Mann basis for a local Hilbert space of dimension `d`.

# Arguments
- `d`: Local Hilbert space dimension (`2S + 1` for spin `S`).

# Returns
- A cached `Vector` of `d^2` dense `d x d` matrices that is Hermitian, traceless apart from
  element 1, and orthonormal under the Hilbert-Schmidt inner product `tr(A' * B)`. Element 1 is
  the normalized identity `I / sqrt(d)`.

# Notes
- The ordering is: the identity; then, for each pair `j < k`, the real and imaginary
  off-diagonal generators; then the `d - 1` diagonal generators. For `d = 2` this reproduces
  `(I, X, Y, Z) / sqrt(2)` exactly, so the `pauli_*` helpers are the `d = 2` case of the
  `operator_*` ones.
- The DMT kernel relies on exactly two properties of this basis: orthonormality, and the
  identity sitting at index 1. Any basis with those properties can be substituted; Hermiticity
  is not required by the algorithm.

# Examples
```jldoctest
julia> length(operator_basis_matrices(3))
9
```
"""
function operator_basis_matrices(d::Integer)
  local_dim = Int(d)
  local_dim >= 2 || throw(ArgumentError("local Hilbert space dimension must be at least 2, got $(local_dim)"))
  return get!(_OPERATOR_BASIS_CACHE, local_dim) do
    _build_operator_basis(local_dim)
  end
end

function _build_operator_basis(d::Int)
  basis = Vector{Matrix{ComplexF64}}()
  push!(basis, Matrix{ComplexF64}(I, d, d) / sqrt(d))
  for j in 1:d, k in (j + 1):d
    symmetric = zeros(ComplexF64, d, d)
    symmetric[j, k] = 1 / sqrt(2)
    symmetric[k, j] = 1 / sqrt(2)
    push!(basis, symmetric)
    antisymmetric = zeros(ComplexF64, d, d)
    antisymmetric[j, k] = -im / sqrt(2)
    antisymmetric[k, j] = im / sqrt(2)
    push!(basis, antisymmetric)
  end
  for level in 1:(d - 1)
    entries = zeros(ComplexF64, d)
    for j in 1:level
      entries[j] = 1
    end
    entries[level + 1] = -level
    push!(basis, Matrix(Diagonal(entries)) / sqrt(level * (level + 1)))
  end
  return basis
end

"""
    local_dimension(site_dimension)
    local_dimension(site::Index)

Return the physical local dimension `d` encoded by an operator-space site of dimension `d^2`.

# Arguments
- `site_dimension`: Operator-space site dimension, or an `Index` carrying it.

# Returns
- The integer `d` with `d^2 == site_dimension`.

# Notes
- Operator-space sites always have a perfect-square dimension because they carry a basis for
  the `d x d` onsite operator space. This is how every `operator_*` helper recovers `d` without
  being told.
"""
function local_dimension(site_dimension::Integer)
  squared = Int(site_dimension)
  d = isqrt(squared)
  d * d == squared || throw(ArgumentError("operator-space site dimension $(squared) is not a perfect square"))
  d >= 2 || throw(ArgumentError("operator-space site dimension $(squared) implies local dimension $(d) < 2"))
  return d
end

local_dimension(site::Index) = local_dimension(dim(site))
