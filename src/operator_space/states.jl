"""
    operator_siteinds(nsites; d=2, tagprefix="OperatorSpace")

Construct site indices for a vectorized operator-space `MPS` on `nsites` sites of local
Hilbert space dimension `d`.

# Arguments
- `nsites`: Number of operator-space sites.

# Keyword Arguments
- `d`: Local Hilbert space dimension (`2S + 1` for spin `S`). The generated indices have
  dimension `d^2`.
- `tagprefix`: Prefix used when naming the generated `Index` tags.

# Returns
- A vector of length `nsites` containing dimension-`d^2` `Index` objects.

# Notes
- The basis ordering on each site is [`operator_basis_matrices`](@ref)`(d)`, whose first
  element is the normalized identity. Downstream helpers recover `d` from the site dimension
  via [`local_dimension`](@ref), so it never has to be passed again.
"""
function operator_siteinds(nsites::Integer; d::Integer=2, tagprefix::AbstractString="OperatorSpace")
  nsites >= 1 || throw(ArgumentError("number of operator-space sites must be positive"))
  local_dim = Int(d)
  local_dim >= 2 || throw(ArgumentError("local Hilbert space dimension must be at least 2"))
  return [Index(local_dim^2, "$(tagprefix),n=$(n)") for n in 1:Int(nsites)]
end

"""
    _operator_basis_label(label, d)

Normalize one local operator-basis label into its integer index in the
[`operator_basis_matrices`](@ref) ordering. Pauli letters are accepted only at `d == 2`.
"""
function _operator_basis_label(label::Integer, d::Integer)
  1 <= label <= Int(d)^2 || throw(ArgumentError("operator-basis labels must lie in 1:$(Int(d)^2) for local dimension $(d)"))
  return Int(label)
end

function _operator_basis_label(label::Symbol, d::Integer)
  return _operator_basis_label(String(label), d)
end

function _operator_basis_label(label::AbstractString, d::Integer)
  normalized = uppercase(strip(label))
  normalized == "I" && return 1
  Int(d) == 2 || throw(ArgumentError("letter operator-basis labels are only defined for local dimension 2; use an integer label in 1:$(Int(d)^2)"))
  normalized == "X" && return 2
  normalized == "Y" && return 3
  normalized == "Z" && return 4
  throw(ArgumentError("unsupported operator-basis label: $(label)"))
end

"""
    _operator_coefficients(op, d)

Return the coefficient vector of a dense `d x d` operator in the normalized basis, i.e. the
entries `tr(P_mu' * op)` for `P_mu = operator_basis_matrices(d)[mu]`.
"""
function _operator_coefficients(op::AbstractMatrix, d::Integer)
  local_dim = Int(d)
  size(op) == (local_dim, local_dim) ||
    throw(ArgumentError("local operator must be $(local_dim) x $(local_dim), got $(size(op))"))
  basis = operator_basis_matrices(local_dim)
  dense = Matrix{ComplexF64}(op)
  return ComplexF64[tr(adjoint(matrix) * dense) for matrix in basis]
end

"""
    operator_basis_state(sites, labels; coefficient=1.0)

Build a product `MPS` in operator space selecting one basis element per site.

# Arguments
- `sites`: Operator-space site indices, typically from [`operator_siteinds`](@ref).
- `labels`: One local basis label per site: an integer in `1:d^2`, or `"I"` (any `d`), or
  `"X"`, `"Y"`, `"Z"` when `d == 2`.

# Keyword Arguments
- `coefficient`: Overall scalar prefactor stored on the first tensor.

# Returns
- A product-state `MPS` in operator space.
"""
function operator_basis_state(sites, labels::AbstractVector; coefficient::Number=1.0)
  length(sites) == length(labels) || throw(ArgumentError("operator-basis labels must have one entry per site"))
  tensors = ITensor[]
  for (index, (site, label)) in enumerate(zip(sites, labels))
    d = local_dimension(site)
    tensor = ITensor(site)
    tensor[site => _operator_basis_label(label, d)] = index == 1 ? coefficient : 1.0
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

"""
    operator_product_state(sites, ops; coefficient=1.0)

Build the operator-space `MPS` for the tensor product `⊗_j O_j` of dense local operators.

# Arguments
- `sites`: Operator-space site indices.
- `ops`: One dense `d x d` matrix per site.

# Keyword Arguments
- `coefficient`: Overall scalar prefactor stored on the first tensor.

# Returns
- A bond-dimension-1 `MPS` whose amplitude on the basis product `alpha` is
  `prod_j tr(P_{alpha_j}' * O_j)`.

# Notes
- Use this rather than [`operator_basis_state`](@ref) whenever the operator you want is not a
  single basis element. At `d = 3` the physical `S^z = diag(1, 0, -1)` is **not** proportional
  to any single Gell-Mann matrix, so a spin-1 `S^z` string must be built this way.
"""
function operator_product_state(sites, ops; coefficient::Number=1.0)
  length(sites) == length(ops) || throw(ArgumentError("operator_product_state needs one local operator per site"))
  tensors = ITensor[]
  for (index, (site, op)) in enumerate(zip(sites, ops))
    d = local_dimension(site)
    coefficients = _operator_coefficients(op, d)
    tensor = ITensor(ComplexF64, site)
    scale = index == 1 ? ComplexF64(coefficient) : one(ComplexF64)
    for mu in eachindex(coefficients)
      iszero(coefficients[mu]) && continue
      tensor[site => mu] = scale * coefficients[mu]
    end
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

"""
    operator_local_sum_state(sites, op, coeffs)

Build the operator-space `MPS` representing the local-density sum `sum_j c_j O_j`, where `O_j`
is the dense local operator `op` placed on site `j` with identity elsewhere.

# Arguments
- `sites`: Operator-space site indices.
- `op`: Dense `d x d` local operator.
- `coeffs`: One coefficient `c_j` per site.

# Returns
- An `MPS` of bond dimension 2 (link state 1 = "operator not yet placed", state 2 = "placed,
  identity henceforth").

# Notes
- This is the transport source builder: uniform `coeffs` give a total charge such as
  `sum_j S^z_j`, and a sign-split profile gives the domain-wall operator whose melting sets the
  dynamical exponent. `pauli_total_sz_state` and `pauli_domain_wall_state` are the `d = 2` cases.
- The identity coefficient inserted on unoccupied sites is `1 / sqrt(d)` times `sqrt(d)`, i.e.
  the amplitude `1` on basis element 1, matching the normalized-basis convention used by
  [`operator_trace`](@ref).
"""
function operator_local_sum_state(sites, op::AbstractMatrix, coeffs::AbstractVector)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("operator_local_sum_state requires at least one site"))
  length(coeffs) == nsites || throw(ArgumentError("operator_local_sum_state needs one coefficient per site"))
  d = local_dimension(first(sites))
  all(local_dimension(site) == d for site in sites) ||
    throw(ArgumentError("operator_local_sum_state requires a uniform local dimension"))
  weights = _operator_coefficients(op, d)

  function placed!(tensor, indices, site, scale)
    for mu in eachindex(weights)
      iszero(weights[mu]) && continue
      tensor[(indices..., site => mu)...] = scale * weights[mu]
    end
    return tensor
  end

  if nsites == 1
    tensor = ITensor(ComplexF64, sites[1])
    placed!(tensor, (), sites[1], ComplexF64(coeffs[1]))
    return MPS([tensor])
  end

  left_link = Index(2, "OperatorStateLink,n=1")
  first_tensor = ITensor(ComplexF64, sites[1], left_link)
  first_tensor[sites[1] => 1, left_link => 1] = 1.0
  placed!(first_tensor, (left_link => 2,), sites[1], ComplexF64(coeffs[1]))
  tensors = ITensor[first_tensor]

  for j in 2:(nsites - 1)
    right_link = Index(2, "OperatorStateLink,n=$(j)")
    tensor = ITensor(ComplexF64, dag(left_link), sites[j], right_link)
    tensor[dag(left_link) => 1, sites[j] => 1, right_link => 1] = 1.0
    tensor[dag(left_link) => 2, sites[j] => 1, right_link => 2] = 1.0
    placed!(tensor, (dag(left_link) => 1, right_link => 2), sites[j], ComplexF64(coeffs[j]))
    push!(tensors, tensor)
    left_link = right_link
  end

  last_tensor = ITensor(ComplexF64, dag(left_link), sites[nsites])
  last_tensor[dag(left_link) => 2, sites[nsites] => 1] = 1.0
  placed!(last_tensor, (dag(left_link) => 1,), sites[nsites], ComplexF64(coeffs[nsites]))
  push!(tensors, last_tensor)
  return MPS(tensors)
end

"""
    pauli_siteinds(nsites; tagprefix="PauliSpace")

Spin-1/2 case of [`operator_siteinds`](@ref): dimension-4 sites in the `(I, X, Y, Z)` ordering.
"""
function pauli_siteinds(nsites::Integer; tagprefix::AbstractString="PauliSpace")
  return operator_siteinds(nsites; d=2, tagprefix=tagprefix)
end

"""
    pauli_basis_state(sites, labels; coefficient=1.0)

Spin-1/2 case of [`operator_basis_state`](@ref), with labels in the `(I, X, Y, Z)` ordering.
"""
function pauli_basis_state(sites, labels::AbstractVector; coefficient::Number=1.0)
  return operator_basis_state(sites, labels; coefficient=coefficient)
end

"""
    pauli_total_sz_state(sites; coefficient=nothing)

Spin-1/2 case of [`operator_local_sum_state`](@ref) for the uniform total-`S^z` operator.
The default per-site coefficient is `2^(N / 2 - 1)`, matching the normalized Pauli convention.
"""
function pauli_total_sz_state(sites; coefficient::Union{Nothing,Number}=nothing)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("pauli_total_sz_state requires at least one site"))
  weight = isnothing(coefficient) ? 2.0^(nsites / 2 - 1) : coefficient
  return operator_local_sum_state(sites, pauli_matrices().Z / sqrt(2), fill(weight, nsites))
end

"""
    pauli_domain_wall_state(sites; kink=length(sites) ÷ 2, coefficient=1.0)

Spin-1/2 case of [`operator_local_sum_state`](@ref) for the signed domain-wall operator
`D = sum_j sign(j - kink) * sigma^z_j`, the infinite-temperature transport source.
`D` is traceless, so measure its profile with `normalize=false`.
"""
function pauli_domain_wall_state(sites; kink::Integer=length(sites) ÷ 2, coefficient::Number=1.0)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("pauli_domain_wall_state requires at least one site"))
  0 <= kink <= nsites || throw(ArgumentError("kink must lie in 0:length(sites)"))
  weights = [j <= kink ? -coefficient : coefficient for j in 1:nsites]
  return operator_local_sum_state(sites, pauli_matrices().Z / sqrt(2), weights)
end
