"""
    pauli_matrices(; include_identity=true)

Return the dense spin-1/2 Pauli matrices.

# Keyword Arguments
- `include_identity`: If `true`, include the identity matrix as field `I`.

# Returns
- A named tuple whose fields follow the ordering `(I, X, Y, Z)` when the identity is
  included, and `(X, Y, Z)` otherwise.

# Examples
```julia
julia> pauli_matrices().Z
2x2 Matrix{ComplexF64}:
 1+0im   0+0im
 0+0im  -1+0im
```
"""
function pauli_matrices(; include_identity::Bool=true)
  matrices = (
    I=ComplexF64[1 0; 0 1],
    X=ComplexF64[0 1; 1 0],
    Y=ComplexF64[0 -im; im 0],
    Z=ComplexF64[1 0; 0 -1],
  )
  return include_identity ? matrices : (X=matrices.X, Y=matrices.Y, Z=matrices.Z)
end

"""
    pauli_basis(; include_identity=true)

Return the spin-1/2 Pauli basis as labeled dense matrices.

# Keyword Arguments
- `include_identity`: Forwarded to [`pauli_matrices`](@ref).

# Returns
- A vector of `Pair{Symbol, Matrix}` entries in the same local ordering used throughout the
  operator-space code.
"""
function pauli_basis(; include_identity::Bool=true)
  matrices = pauli_matrices(; include_identity=include_identity)
  return collect(pairs(matrices))
end

"""
    pauli_components(operator; include_identity=true)

Decompose a `2 x 2` operator into Pauli-basis coefficients.

# Arguments
- `operator`: Dense single-site operator matrix.

# Keyword Arguments
- `include_identity`: If `true`, return the coefficient of the identity alongside the
  Pauli components.

# Returns
- A named tuple of Hilbert-Schmidt coefficients, normalized so that
  `operator == sum(coeffs[name] * pauli_matrices()[name] for name in keys(coeffs))`.
"""
function pauli_components(operator::AbstractMatrix; include_identity::Bool=true)
  size(operator) == (2, 2) || throw(ArgumentError("Pauli decomposition requires a 2 x 2 operator"))
  basis = pauli_matrices(; include_identity=include_identity)
  return NamedTuple{keys(basis)}(Tuple(tr(adjoint(matrix) * operator) / 2 for matrix in values(basis)))
end
