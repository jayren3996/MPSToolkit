"""
    pauli_matrices(; include_identity=true)

Return the spin-1/2 Pauli matrices as a named tuple of dense `ComplexF64` matrices.
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

Return the spin-1/2 Pauli basis as a vector of labeled matrix pairs.
"""
function pauli_basis(; include_identity::Bool=true)
  matrices = pauli_matrices(; include_identity=include_identity)
  return collect(pairs(matrices))
end

"""
    pauli_components(operator; include_identity=true)

Decompose a `2 x 2` operator into Pauli-basis coefficients using the Hilbert-Schmidt inner product.
"""
function pauli_components(operator::AbstractMatrix; include_identity::Bool=true)
  size(operator) == (2, 2) || throw(ArgumentError("Pauli decomposition requires a 2 x 2 operator"))
  basis = pauli_matrices(; include_identity=include_identity)
  return NamedTuple{keys(basis)}(Tuple(tr(adjoint(matrix) * operator) / 2 for matrix in values(basis)))
end
