"""
    pauli_basis_state(sites, labels; coefficient=1.0)

Build a Pauli-basis product `MPS` in the default local ordering `(I, X, Y, Z)`.
Each entry of `labels` may be an integer basis index or one of `"I"`, `"X"`, `"Y"`, `"Z"`
(or the corresponding symbols).

# Arguments
- `sites`: Pauli-space site indices, typically created by [`pauli_siteinds`](@ref).
- `labels`: One local basis label per site.

# Keyword Arguments
- `coefficient`: Overall scalar prefactor stored on the first tensor.

# Returns
- A product-state `MPS` in operator space.

# Examples
```julia
sites = pauli_siteinds(3)
rho = pauli_basis_state(sites, [:I, :Z, :I])
```
"""
function pauli_basis_state(sites, labels::AbstractVector; coefficient::Number=1.0)
  length(sites) == length(labels) || throw(ArgumentError("Pauli-basis labels must have one entry per site"))

  tensors = ITensor[]
  for (index, (site, label)) in enumerate(zip(sites, labels))
    tensor = ITensor(site)
    tensor[site => _pauli_basis_label(label)] = index == 1 ? coefficient : 1.0
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

"""
    pauli_total_sz_state(sites; coefficient=0.5)

Build the Pauli-basis `MPS` representing `coefficient * sum_j σ_j^z`.
With the default `coefficient=0.5`, this is the total spin operator `∑_j S_j^z`.

# Arguments
- `sites`: Pauli-space site indices.

# Keyword Arguments
- `coefficient`: Scalar multiplying each local `σ^z` contribution.

# Returns
- An `MPS` representation of the summed operator in Pauli space.

# Notes
- The returned state is not a simple product state; it uses bond dimension `2` to encode
  the operator sum compactly.
"""
function pauli_total_sz_state(sites; coefficient::Number=0.5)
  nsites = length(sites)
  nsites < 1 && throw(ArgumentError("pauli_total_sz_state requires at least one site"))
  nsites == 1 && return pauli_basis_state(sites, [4]; coefficient=coefficient)

  left_link = Index(2, "OperatorStateLink,n=1")
  first = ITensor(sites[1], left_link)
  first[sites[1] => 1, left_link => 1] = 1.0
  first[sites[1] => 4, left_link => 2] = coefficient

  tensors = ITensor[first]

  for j in 2:(nsites - 1)
    right_link = Index(2, "OperatorStateLink,n=$(j)")
    tensor = ITensor(dag(left_link), sites[j], right_link)
    tensor[dag(left_link) => 1, sites[j] => 1, right_link => 1] = 1.0
    tensor[dag(left_link) => 1, sites[j] => 4, right_link => 2] = coefficient
    tensor[dag(left_link) => 2, sites[j] => 1, right_link => 2] = 1.0
    push!(tensors, tensor)
    left_link = right_link
  end

  last = ITensor(dag(left_link), sites[nsites])
  last[dag(left_link) => 1, sites[nsites] => 4] = coefficient
  last[dag(left_link) => 2, sites[nsites] => 1] = 1.0
  push!(tensors, last)

  return MPS(tensors)
end

"""
    pauli_gate(unitary)

Convert a dense physical unitary acting on spin-1/2 sites into the corresponding dense
superoperator in the default local Pauli ordering `(I, X, Y, Z)`.

# Arguments
- `unitary`: Dense physical-space unitary acting on one or more spin-1/2 sites.

# Returns
- A dense matrix representing the induced operator-space map in the Pauli basis.

# Notes
- The returned matrix entry `G[α_out, α_in]` is `tr(P_{α_out}' * U * P_{α_in} * U')`,
  where `P_α` are the normalized `n`-site Pauli basis operators produced by
  `_pauli_basis_operators`. This convention matches `pauli_basis_state`,
  `pauli_siteinds`, and `pauli_lindblad_generator`, so the resulting
  super-operator can be applied to a Pauli-basis MPS via `tebd_evolve!`
  consistently.
"""
function pauli_gate(unitary::AbstractMatrix)
  size(unitary, 1) == size(unitary, 2) || throw(ArgumentError("Pauli-basis gate requires a square unitary"))
  nsites = _spinhalf_span(size(unitary, 1))
  basis_operators = _pauli_basis_operators(nsites)
  u = ComplexF64.(unitary)
  u_dag = adjoint(u)
  nbasis = length(basis_operators)
  gate = Matrix{ComplexF64}(undef, nbasis, nbasis)
  for column in eachindex(basis_operators)
    evolved = u * basis_operators[column] * u_dag
    for row in eachindex(basis_operators)
      gate[row, column] = tr(adjoint(basis_operators[row]) * evolved)
    end
  end
  return gate
end

"""
    pauli_gate_from_hamiltonian(h, dt)

Build the Pauli-basis superoperator induced by a dense Hamiltonian over one time step.

# Arguments
- `h`: Dense physical-space Hamiltonian.
- `dt`: Real-time step.

# Returns
- The operator-space gate corresponding to `exp(-im * dt * h)`.
"""
function pauli_gate_from_hamiltonian(h::AbstractMatrix, dt::Real)
  size(h, 1) == size(h, 2) || throw(ArgumentError("Pauli-basis Hamiltonian must be square"))
  return pauli_gate(exp(-im * dt * h))
end

"""
    pauli_lindblad_generator(h, jumps)

Build the dense Pauli-basis Lindbladian generator induced by the local Hamiltonian `h` and
the local jump operators `jumps`. The local operator basis is assumed to be ordered as
`(I, X, Y, Z)` on each spin-1/2 site.

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
function pauli_lindblad_generator(h::AbstractMatrix, jumps)
  size(h, 1) == size(h, 2) || throw(ArgumentError("Pauli-basis Hamiltonian must be square"))
  nsites = _spinhalf_span(size(h, 1))
  basis_operators = _pauli_basis_operators(nsites)
  local_dim = size(h, 1)
  jump_list = _lindblad_jump_list(jumps)
  for jump in jump_list
    size(jump, 1) == local_dim && size(jump, 2) == local_dim || throw(ArgumentError("jump operators must match the Hamiltonian dimension"))
  end

  lindblad_action = operator -> begin
    updated = -im * (h * operator - operator * h)
    for jump in jump_list
      jump_dag_jump = jump' * jump
      updated += jump * operator * jump'
      updated -= 0.5 * (jump_dag_jump * operator + operator * jump_dag_jump)
    end
    return updated
  end

  generator = Matrix{ComplexF64}(undef, length(basis_operators), length(basis_operators))
  for column in eachindex(basis_operators)
    evolved = lindblad_action(basis_operators[column])
    for row in eachindex(basis_operators)
      generator[row, column] = tr(basis_operators[row]' * evolved)
    end
  end
  return generator
end

"""
    pauli_gate_from_lindbladian(h, jumps, dt)

Build the dense Pauli-basis TEBD gate generated by the local Lindbladian defined by `h`
and `jumps` over one time step `dt`.

# Arguments
- `h`: Dense local Hamiltonian.
- `jumps`: Jump operator or collection of jump operators.
- `dt`: Time increment.

# Returns
- `exp(dt * pauli_lindblad_generator(h, jumps))`.
"""
function pauli_gate_from_lindbladian(h::AbstractMatrix, jumps, dt::Real)
  return exp(dt * pauli_lindblad_generator(h, jumps))
end

"""
    _pauli_basis_label(label)

Normalize one local Pauli-basis label into its integer index in the `(I, X, Y, Z)` ordering.
"""
function _pauli_basis_label(label::Integer)
  1 <= label <= 4 || throw(ArgumentError("Pauli-basis labels must lie in 1:4"))
  return Int(label)
end

"""
    _pauli_basis_label(label::Symbol)

Symbol overload of [`_pauli_basis_label`](@ref).
"""
function _pauli_basis_label(label::Symbol)
  return _pauli_basis_label(String(label))
end

"""
    _pauli_basis_label(label::AbstractString)

String overload of [`_pauli_basis_label`](@ref). Whitespace is stripped and case is ignored.
"""
function _pauli_basis_label(label::AbstractString)
  normalized = uppercase(strip(label))
  normalized == "I" && return 1
  normalized == "X" && return 2
  normalized == "Y" && return 3
  normalized == "Z" && return 4
  throw(ArgumentError("unsupported Pauli-basis label: $(label)"))
end

"""
    _spinhalf_span(dim)

Infer how many spin-1/2 sites are represented by a Hilbert-space dimension `dim`.
"""
function _spinhalf_span(dim::Integer)
  dim < 1 && throw(ArgumentError("spin-1/2 operator-space helpers require a positive matrix dimension"))
  span = round(Int, log2(dim))
  2^span == dim || throw(ArgumentError("matrix dimension must be compatible with spin-1/2 sites"))
  return span
end

"""
    _pauli_basis_operators(nsites)

Enumerate the normalized Pauli-basis operators for `nsites` spin-1/2 sites.
"""
function _pauli_basis_operators(nsites::Integer)
  local_basis = [matrix / sqrt(2) for matrix in values(pauli_matrices())]
  operators = local_basis
  for _ in 2:nsites
    operators = [kron(left, right) for left in operators for right in local_basis]
  end
  return operators
end

"""
    _lindblad_jump_list(jumps)

Normalize the `jumps` argument into a plain Julia vector of dense matrices.
"""
function _lindblad_jump_list(jumps::AbstractMatrix)
  return AbstractMatrix[jumps]
end

"""
    _lindblad_jump_list(jumps::AbstractVector)

Vector overload of [`_lindblad_jump_list`](@ref).
"""
function _lindblad_jump_list(jumps::AbstractVector)
  return collect(jumps)
end
