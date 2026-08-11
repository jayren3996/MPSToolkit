"""
    _operator_span(matrix_dim, d)

Infer how many sites of local dimension `d` a dense physical matrix of size
`matrix_dim x matrix_dim` acts on.
"""
function _operator_span(matrix_dim::Integer, d::Integer)
  local_dim = Int(d)
  local_dim >= 2 || throw(ArgumentError("local Hilbert space dimension must be at least 2"))
  size_value = Int(matrix_dim)
  size_value >= local_dim || throw(ArgumentError("operator dimension $(size_value) is smaller than the local dimension $(local_dim)"))
  span = 1
  product = local_dim
  while product < size_value
    product *= local_dim
    span += 1
  end
  product == size_value ||
    throw(ArgumentError("matrix dimension $(size_value) is not a power of the local dimension $(local_dim)"))
  return span
end

"""
    _operator_basis_operators(nsites, d)

Enumerate the normalized multi-site operator-basis matrices for `nsites` sites of local
dimension `d`, ordered so the site-`1` label varies slowest (matching `kron`).
"""
function _operator_basis_operators(nsites::Integer, d::Integer)
  local_basis = operator_basis_matrices(d)
  operators = local_basis
  for _ in 2:Int(nsites)
    operators = [kron(left, right) for left in operators for right in local_basis]
  end
  return operators
end

"""
    operator_gate(op; d=2)

Convert a dense physical map `A` on one or more sites of local dimension `d` into the
corresponding dense superoperator `rho -> A rho A'` in the normalized operator basis.

# Arguments
- `op`: Dense (or sparse) physical matrix of size `d^span x d^span`. Unitary for real-time
  evolution; Hermitian positive for an imaginary-time slice.

# Keyword Arguments
- `d`: Local Hilbert space dimension.

# Returns
- A dense `d^(2 span) x d^(2 span)` matrix `G` with `G[a, b] = tr(P_a' * A * P_b * A')`,
  where `P` are the multi-site operators of [`operator_basis_matrices`](@ref).

# Notes
- Built as `G = W' * kron(conj(A), A) * W`, where `W[:, b] = vec(P_b)`. This costs
  `O(d^(6 span))` rather than the `O(d^(7 span))` of a per-column loop, which matters once
  `d^span` exceeds a few: a three-site spin-1 gate is `19683 x 19683` in operator space, so
  keep `span` small at higher spin.
- Sparse input (for example from `EDKit.spin(...; D = d)`) is densified internally.
"""
function operator_gate(op::AbstractMatrix; d::Integer=2)
  size(op, 1) == size(op, 2) || throw(ArgumentError("operator-space gate requires a square matrix"))
  local_dim = Int(d)
  span = _operator_span(size(op, 1), local_dim)
  dense = Matrix{ComplexF64}(op)
  basis = _operator_basis_operators(span, local_dim)
  nbasis = length(basis)
  physical_dim = size(dense, 1)
  # W[:, b] = vec(P_b) with the same column-major vec used by kron(conj(A), A).
  w = Matrix{ComplexF64}(undef, physical_dim^2, nbasis)
  for b in 1:nbasis
    w[:, b] = vec(basis[b])
  end
  return _real_operator_data(adjoint(w) * (kron(conj(dense), dense) * w))
end

# Set to `true` to suppress the real downcast in `_real_operator_data`, forcing every gate and
# coefficient vector to stay `ComplexF64`. This exists so that "the complex path is unchanged" can be established
# mechanically -- flip it on and the whole suite must reproduce the numbers it produced before the
# element type was threaded -- rather than by reading diffs.
const _OPERATOR_GATE_FORCE_COMPLEX = Ref(false)

"""
    _real_operator_data(x)

Return `real(x)` when the imaginary part of `x` is numerical noise, and `x` unchanged otherwise.
Used for both operator-space gate matrices and operator-space coefficient vectors.

# Notes
- This recovers an exact property rather than approximating one. In a basis of *Hermitian*
  operators every Hermiticity-preserving superoperator has a real matrix:
  `conj(G[a, b]) = conj(tr(P_a' A P_b A')) = tr((P_a A P_b A')') = tr(A P_b A' P_a) = G[a, b]`,
  using `P' = P` and the cyclicity of the trace. Note what the derivation does *not* need: `A`
  need be neither unitary nor Hermitian, so the downcast is unconditional rather than a
  Hermiticity heuristic. Measured residues across every builder here are 1e-18 to 1.2e-16
  relative, at `d = 2, 3, 4`.
- The tolerance is `sqrt(eps(Float64))` *relative*, this repo's existing threshold for exactly
  this class of question (see `operator_gate_from_imaginary_time`'s Hermiticity check,
  `_dmt_connector`'s connector-negligibility test, and `_reject_unresolvable_trace`). It sits
  eight decades above the observed residue and far below any physical imaginary part. An absolute
  threshold would misfire on a large-norm generator.
- Falls back to the complex matrix instead of throwing, because
  [`operator_basis_matrices`](@ref) explicitly sanctions substituting a non-Hermitian basis, and
  for one of those a genuinely complex gate is the correct answer.
"""
function _real_operator_data(x::AbstractArray)
  _OPERATOR_GATE_FORCE_COMPLEX[] && return x
  norm(imag(x)) <= sqrt(eps(Float64)) * max(norm(x), one(Float64)) || return x
  return real(x)
end

"""
    operator_gate_from_hamiltonian(h, dt; d=2)

Build the operator-space superoperator induced by `exp(-im * dt * h)` for a dense local
Hamiltonian `h` on sites of local dimension `d`.

# Arguments
- `h`: Dense (or sparse) physical-space Hamiltonian.
- `dt`: Real-time step.

# Keyword Arguments
- `d`: Local Hilbert space dimension.

# Returns
- The operator-space gate corresponding to `exp(-im * dt * h)`, as returned by
  [`operator_gate`](@ref).
"""
function operator_gate_from_hamiltonian(h::AbstractMatrix, dt::Real; d::Integer=2)
  size(h, 1) == size(h, 2) || throw(ArgumentError("operator-space Hamiltonian must be square"))
  return operator_gate(exp(-im * dt * Matrix{ComplexF64}(h)); d=d)
end

"""
    operator_gate_from_imaginary_time(h, dbeta; d=2)

Build the operator-space superoperator for one two-sided imaginary-time slice
`rho -> e^{-(dbeta/2) h} rho e^{-(dbeta/2) h}` for a dense local Hermitian `h`.

# Arguments
- `h`: Dense local Hermitian operator on one or more sites of local dimension `d`.
- `dbeta`: Imaginary-time increment for this slice (each side picks up `dbeta / 2`).

# Keyword Arguments
- `d`: Local Hilbert space dimension.

# Returns
- A dense operator-space gate suitable for `tebd_evolve!`.
"""
function operator_gate_from_imaginary_time(h::AbstractMatrix, dbeta::Real; d::Integer=2)
  size(h, 1) == size(h, 2) || throw(ArgumentError("imaginary-time Hamiltonian must be square"))
  dense = Matrix{ComplexF64}(h)
  norm(dense - dense') <= sqrt(eps(Float64)) * max(norm(dense), one(Float64)) ||
    throw(ArgumentError("imaginary-time Hamiltonian must be Hermitian"))
  return operator_gate(exp(-(dbeta / 2) * dense); d=d)
end

"""
    operator_lindblad_generator(h, jumps; d=2)

Build the dense operator-space Lindbladian generator for a local Hamiltonian `h` and local
jump operators `jumps`, implementing `-im[H, rho] + sum_j (L_j rho L_j' - {L_j' L_j, rho}/2)`.

# Arguments
- `h`: Dense Hamiltonian acting on one or more sites of local dimension `d`.
- `jumps`: One jump operator or a collection of jump operators with the same dimension as
  `h`.

# Keyword Arguments
- `d`: Local Hilbert space dimension.

# Returns
- A dense generator matrix in the normalized operator basis of
  [`operator_basis_matrices`](@ref).
"""
function operator_lindblad_generator(h::AbstractMatrix, jumps; d::Integer=2)
  size(h, 1) == size(h, 2) || throw(ArgumentError("operator-space Hamiltonian must be square"))
  local_dim = Int(d)
  span = _operator_span(size(h, 1), local_dim)
  basis = _operator_basis_operators(span, local_dim)
  dense_h = Matrix{ComplexF64}(h)
  jump_list = _lindblad_jump_list(jumps)
  for jump in jump_list
    size(jump) == size(dense_h) || throw(ArgumentError("jump operators must match the Hamiltonian dimension"))
  end
  jump_data = [(Matrix{ComplexF64}(jump), Matrix{ComplexF64}(jump)' * Matrix{ComplexF64}(jump)) for jump in jump_list]

  generator = Matrix{ComplexF64}(undef, length(basis), length(basis))
  for column in eachindex(basis)
    evolved = -im * (dense_h * basis[column] - basis[column] * dense_h)
    for (jump, jump_dag_jump) in jump_data
      evolved += jump * basis[column] * jump'
      evolved -= 0.5 * (jump_dag_jump * basis[column] + basis[column] * jump_dag_jump)
    end
    for row in eachindex(basis)
      generator[row, column] = tr(basis[row]' * evolved)
    end
  end
  # Real by the same argument as in `_real_operator_data`: the Lindbladian preserves
  # Hermiticity, so its matrix in a Hermitian basis is real. `operator_gate_from_lindbladian`
  # inherits this for free, since `exp(::Matrix{Float64})` is a `Matrix{Float64}`.
  return _real_operator_data(generator)
end

"""
    operator_gate_from_lindbladian(h, jumps, dt; d=2)

Build the dense operator-space TEBD gate `exp(dt * operator_lindblad_generator(h, jumps; d))`.

# Arguments
- `h`: Dense local Hamiltonian.
- `jumps`: Jump operator or collection of jump operators.
- `dt`: Time increment.

# Keyword Arguments
- `d`: Local Hilbert space dimension.

# Returns
- `exp(dt * operator_lindblad_generator(h, jumps; d))`.
"""
function operator_gate_from_lindbladian(h::AbstractMatrix, jumps, dt::Real; d::Integer=2)
  return exp(dt * operator_lindblad_generator(h, jumps; d=d))
end

"""
    pauli_gate(unitary)

Spin-1/2 case of [`operator_gate`](@ref): convert a dense physical unitary acting on
spin-1/2 sites into the corresponding dense superoperator in the default local Pauli
ordering `(I, X, Y, Z)`.

# Arguments
- `unitary`: Dense physical-space unitary acting on one or more spin-1/2 sites.

# Returns
- A dense matrix representing the induced operator-space map in the Pauli basis.

# Notes
- The returned matrix entry `G[α_out, α_in]` is `tr(P_{α_out}' * U * P_{α_in} * U')`,
  where `P_α` are the normalized `n`-site Pauli basis operators produced by
  [`operator_basis_matrices`](@ref)`(2)`. This convention matches `pauli_basis_state`,
  `pauli_siteinds`, and `pauli_lindblad_generator`, so the resulting super-operator can be
  applied to a Pauli-basis MPS via `tebd_evolve!` consistently.
"""
pauli_gate(unitary::AbstractMatrix) = operator_gate(unitary; d=2)

"""
    pauli_gate_from_hamiltonian(h, dt)

Spin-1/2 case of [`operator_gate_from_hamiltonian`](@ref): build the Pauli-basis
superoperator induced by a dense Hamiltonian over one time step.

# Arguments
- `h`: Dense physical-space Hamiltonian.
- `dt`: Real-time step.

# Returns
- The operator-space gate corresponding to `exp(-im * dt * h)`.
"""
pauli_gate_from_hamiltonian(h::AbstractMatrix, dt::Real) = operator_gate_from_hamiltonian(h, dt; d=2)

"""
    pauli_lindblad_generator(h, jumps)

Spin-1/2 case of [`operator_lindblad_generator`](@ref): build the dense Pauli-basis
Lindbladian generator induced by the local Hamiltonian `h` and the local jump operators
`jumps`. The local operator basis is assumed to be ordered as `(I, X, Y, Z)` on each
spin-1/2 site.

# Arguments
- `h`: Dense Hamiltonian acting on one or more spin-1/2 sites.
- `jumps`: One jump operator or a collection of jump operators with the same dimension as
  `h`.

# Returns
- A dense generator matrix in the normalized Pauli basis.

# Notes
- The generator implements the standard Lindblad action
  `-im[H, ρ] + Σ_j (L_j ρ L_j† - 1/2 {L_j†L_j, ρ})`.
"""
pauli_lindblad_generator(h::AbstractMatrix, jumps) = operator_lindblad_generator(h, jumps; d=2)

"""
    pauli_gate_from_lindbladian(h, jumps, dt)

Spin-1/2 case of [`operator_gate_from_lindbladian`](@ref): build the dense Pauli-basis TEBD
gate generated by the local Lindbladian defined by `h` and `jumps` over one time step `dt`.

# Arguments
- `h`: Dense local Hamiltonian.
- `jumps`: Jump operator or collection of jump operators.
- `dt`: Time increment.

# Returns
- `exp(dt * pauli_lindblad_generator(h, jumps))`.
"""
pauli_gate_from_lindbladian(h::AbstractMatrix, jumps, dt::Real) = operator_gate_from_lindbladian(h, jumps, dt; d=2)
