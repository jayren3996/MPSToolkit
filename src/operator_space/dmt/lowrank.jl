"""
    _protected_basis(protected)

Return a thin orthonormal basis whose columns span the columns of `protected`.

# Notes
- Uses a QR factorization rather than an SVD: only the span matters, the input is
  `chi x d^(2n)` with `d^(2n)` tiny compared with `chi`, and a thin QR is orders of magnitude
  cheaper than `svd(...; full=true)` (measured 0.0013 s versus 0.505 s at `chi = 1800`).
- No rank detection is performed. A rank-deficient `protected` yields a few arbitrary extra
  orthonormal columns, which over-protects rather than under-protects, and is therefore safe.
"""
function _protected_basis(protected::AbstractMatrix)
  size(protected, 2) == 0 && return zeros(ComplexF64, size(protected, 1), 0)
  ncols = min(size(protected, 1), size(protected, 2))
  factorization = qr(Matrix{ComplexF64}(protected))
  return Matrix(factorization.Q)[:, 1:ncols]
end

"""
    _bond_mul(bond_matrix, x)
    _bond_adj(bond_matrix, x)

Apply the bond matrix (or its adjoint) to `x`, with an `O(chi)` fast path when the bond matrix
is `Diagonal`. The bond matrix is diagonal when the bond tensor is factorized by SVD and
triangular when it is factorized by QR; the DMT rule is identical either way.
"""
_bond_mul(bond_matrix::Diagonal, x) = diag(bond_matrix) .* x
_bond_mul(bond_matrix::AbstractMatrix, x) = bond_matrix * x
_bond_adj(bond_matrix::Diagonal, x) = conj(diag(bond_matrix)) .* x
_bond_adj(bond_matrix::AbstractMatrix, x) = bond_matrix' * x

"""
    _dmt_connector(bond_matrix, q0, r0)

Return `(a, b, has_connector)` describing the rank-one trace connector
`C = a * transpose(b) = (S r0)(q0' S) / (q0' S r0)` for the bond matrix `S`.

# Notes
- When the identity overlap `q0' S r0` is negligible relative to the operator scale — a
  traceless operator has no trace component to anchor on — the connector is disabled and `a`
  is returned as zeros, so the caller still truncates and still enforces `maxdim`.
- `B = S - C` annihilates the identity directions (`B r0 = 0`, `q0' B = 0`), so the connector
  costs no extra bond dimension: its range already lies inside the protected block.
"""
function _dmt_connector(bond_matrix::AbstractMatrix, q0::AbstractVector, r0::AbstractVector)
  scaled_right = _bond_mul(bond_matrix, r0)
  denominator = dot(q0, scaled_right)
  has_connector = abs(denominator) > sqrt(eps(Float64)) * norm(bond_matrix)
  a = has_connector ? ComplexF64.(scaled_right ./ denominator) : zeros(ComplexF64, length(r0))
  b = ComplexF64.(conj(_bond_adj(bond_matrix, q0)))    # b = (q0' S)^T = conj(S' q0)
  return a, b, has_connector
end

"""
    _dmt_complement_ops(bond_matrix, a, b, QL, QR)

Return matrix-free operators for `D = P_L^perp B P_R^perp` with `B = S - a b^T` and
`P^perp = I - Q Q'`.

# Returns
- A named tuple with `mul(X) = D * X`, `adj(X) = D' * X`, and the reusable intermediates
  `QLtB = QL' * B`, `BQR = B * QR`, `BQRc = B QR - QL (QL' B QR)`.

# Notes
- Nothing here is `chi x chi`. With a diagonal bond matrix every product costs
  `O(chi * ncols * d^(2n))`, which is what makes the randomized truncation in
  [`_truncated_svd`](@ref) cheap; with a triangular one it is `O(chi^2 * ncols)`.
"""
function _dmt_complement_ops(bond_matrix::AbstractMatrix, a::AbstractVector, b::AbstractVector,
                             QL::AbstractMatrix, QR::AbstractMatrix)
  bmul(x) = _bond_mul(bond_matrix, x) .- a * (transpose(b) * x)
  badj(x) = _bond_adj(bond_matrix, x) .- conj(b) * (a' * x)
  qltb = badj(QL)'
  bqr = bmul(QR)
  bqrc = bqr - QL * (qltb * QR)
  mul(x) = bmul(x) - QL * (qltb * x) - bqrc * (QR' * x)
  function adj(x)
    w = badj(x) - qltb' * (QL' * x)
    return w - QR * (QR' * w)
  end
  return (mul=mul, adj=adj, QLtB=qltb, BQR=bqr, BQRc=bqrc)
end

"""
    _truncated_svd(mul, adj, chi, rank; mode=:dense, oversample=10, power=2)

Return the leading `rank` singular triplet `(U, S, V)` of the `chi x chi` operator defined by
the matrix-free products `mul` and `adj`.

# Keyword Arguments
- `mode`: `:dense` materializes the operator and calls LAPACK; `:random` uses a block
  randomized range finder with `power` power iterations and `oversample` extra probe vectors.

# Notes
- The randomized branch cannot break DMT's preservation guarantee: its range lies inside the
  range of the operator, which is already doubly orthogonal to the protected subspaces, so the
  protected border and connector are still reinstated exactly. It only affects how optimal the
  discarded weight is (measured overlap 0.998-0.999 with the dense result).
"""
function _truncated_svd(mul, adj, chi::Integer, rank::Integer; mode::Symbol=:dense,
                        oversample::Integer=10, power::Integer=2)
  keep = max(min(Int(rank), Int(chi)), 0)
  keep == 0 && return (zeros(ComplexF64, chi, 0), Float64[], zeros(ComplexF64, chi, 0))
  if mode === :dense
    factorization = svd(mul(Matrix{ComplexF64}(I, chi, chi)))
    return (factorization.U[:, 1:keep], factorization.S[1:keep], factorization.V[:, 1:keep])
  end
  mode === :random || throw(ArgumentError("DMT truncated SVD mode must be :dense or :random"))
  probes = min(Int(chi), keep + Int(oversample))
  sketch = mul(randn(ComplexF64, Int(chi), probes))
  for _ in 1:Int(power)
    sketch = mul(adj(Matrix(qr(sketch).Q)))
  end
  range_basis = Matrix(qr(sketch).Q)
  projected = adj(range_basis)'                   # probes x chi
  factorization = svd(projected)
  keep = min(keep, length(factorization.S))
  return (range_basis * factorization.U[:, 1:keep], factorization.S[1:keep], factorization.V[:, 1:keep])
end

"""
    _dmt_refactor(left, right, maxdim, cutoff)

Return the SVD `(U, S, V)` of the low-rank product `left * right'` without forming it.

# Notes
- Costs `O(chi r^2 + r^3)` for `chi x r` inputs, versus `O(chi^3)` for a dense SVD (measured
  0.60 s versus 14.2 s at `chi = 3600, r = 400`). The DMT bond matrix is always available in
  this factored form, so the dense matrix is never needed.
"""
function _dmt_refactor(left::AbstractMatrix, right::AbstractMatrix, maxdim::Integer, cutoff::Real)
  size(left, 1) == size(right, 1) || throw(ArgumentError("DMT refactorization needs matching row counts"))
  size(left, 2) == size(right, 2) || throw(ArgumentError("DMT refactorization needs matching column counts"))
  left_qr = qr(Matrix{ComplexF64}(left))
  right_qr = qr(Matrix{ComplexF64}(right))
  core = svd(Matrix(left_qr.R) * Matrix(right_qr.R)')
  keep = length(core.S)
  if keep > 0 && cutoff > 0
    keep = max(count(>(cutoff * core.S[1]), core.S), 1)
  end
  keep = min(keep, Int(maxdim))
  return (Matrix(left_qr.Q) * core.U[:, 1:keep], core.S[1:keep], Matrix(right_qr.Q) * core.V[:, 1:keep])
end
