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

"""
    _operator_basis_label(label::Symbol, d)

Symbol overload of [`_operator_basis_label`](@ref).
"""
function _operator_basis_label(label::Symbol, d::Integer)
  return _operator_basis_label(String(label), d)
end

"""
    _operator_basis_label(label::AbstractString, d)

String overload of [`_operator_basis_label`](@ref). Whitespace is stripped and case is ignored.
"""
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

# Notes
- Real for a Hermitian `op`, because the basis is Hermitian and orthonormal: `tr(P_mu' op)` is
  then the expansion coefficient of a Hermitian operator on a Hermitian basis element (measured
  imaginary part exactly `0.0`; a non-Hermitian `op` gives 0.89-0.98 relative and stays complex).
  This function is the pivot for the element type of every state builder below, so downcasting
  here is what lets a Hermitian density matrix be represented, evolved and truncated entirely in
  `Float64`. See [`_real_operator_data`](@ref) for the tolerance and the guarantee.
"""
function _operator_coefficients(op::AbstractMatrix, d::Integer)
  local_dim = Int(d)
  size(op) == (local_dim, local_dim) ||
    throw(ArgumentError("local operator must be $(local_dim) x $(local_dim), got $(size(op))"))
  basis = operator_basis_matrices(local_dim)
  dense = Matrix{ComplexF64}(op)
  return _real_operator_data(ComplexF64[tr(adjoint(matrix) * dense) for matrix in basis])
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
  # One element type for every tensor. Allocating `ITensor(site)` and letting each assignment pick
  # its own type made this a mixed-`MPS` producer: with a complex `coefficient` only tensor 1
  # promoted, leaving `[ComplexF64, Float64, ...]`, which costs a materialized promoted copy at
  # every contraction that spans the boundary (see `_mps_eltype`).
  elt = float(typeof(coefficient))
  tensors = ITensor[]
  for (index, (site, label)) in enumerate(zip(sites, labels))
    d = local_dimension(site)
    tensor = ITensor(elt, site)
    tensor[site => _operator_basis_label(label, d)] = index == 1 ? coefficient : one(elt)
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
  # Expand every site first, then pick one element type for the whole chain. Per-site typing would
  # make a chain of mostly-Hermitian operators mixed the moment one site is non-Hermitian, and a
  # mixed `MPS` pays a materialized promoted copy at the boundary (see `_mps_eltype`).
  expanded = [_operator_coefficients(op, local_dimension(site)) for (site, op) in zip(sites, ops)]
  elt = float(promote_type(typeof(coefficient), eltype.(expanded)...))
  tensors = ITensor[]
  for (index, (site, coefficients)) in enumerate(zip(sites, expanded))
    tensor = ITensor(elt, site)
    scale = index == 1 ? convert(elt, coefficient) : one(elt)
    for mu in eachindex(coefficients)
      iszero(coefficients[mu]) && continue
      tensor[site => mu] = scale * coefficients[mu]
    end
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

"""
    _local_sum_state(sites, weights, coeffs, identity_amplitude)

Shared bond-dimension-2 tensor-network construction behind [`operator_local_sum_state`](@ref)
and the `d = 2` wrappers [`pauli_total_sz_state`](@ref)/[`pauli_domain_wall_state`](@ref): place
a local operator with normalized-basis coefficients `weights` (as returned by
[`_operator_coefficients`](@ref)) on site `j` with prefactor `coeffs[j]`, and insert amplitude
`identity_amplitude` on basis element 1 ("identity direction") on every other site.

# Notes
- `identity_amplitude` is the only free parameter that distinguishes the callers:
  `operator_local_sum_state` passes `sqrt(d)`, which reproduces a literal `d x d` identity on
  every unoccupied site (so the `MPS` equals the literal operator sum exactly); the `pauli_*`
  wrappers pass `1.0`, preserving their pre-existing bond-dimension-2 amplitudes. The two
  differ by exactly `identity_amplitude^(nsites - 1)`, since every unoccupied site contributes
  one factor and the occupied site contributes none.
- For `nsites == 1` there is no unoccupied site, so `identity_amplitude` never enters and the
  single-site branch is identical for every caller.
"""
function _local_sum_state(sites, weights::AbstractVector, coeffs::AbstractVector, identity_amplitude::Number)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("_local_sum_state requires at least one site"))
  length(coeffs) == nsites || throw(ArgumentError("_local_sum_state needs one coefficient per site"))

  function placed!(tensor, indices, site, scale)
    for mu in eachindex(weights)
      iszero(weights[mu]) && continue
      tensor[(indices..., site => mu)...] = scale * weights[mu]
    end
    return tensor
  end

  # One element type across the chain, promoted over everything that gets stored: a Hermitian
  # local operator has real `weights` (see `_operator_coefficients`), so a real-coefficient sum
  # stays entirely in `Float64`.
  elt = float(promote_type(eltype(weights), eltype(coeffs), typeof(identity_amplitude)))

  if nsites == 1
    tensor = ITensor(elt, sites[1])
    placed!(tensor, (), sites[1], convert(elt, coeffs[1]))
    return MPS([tensor])
  end

  amplitude = convert(elt, identity_amplitude)
  left_link = Index(2, "OperatorStateLink,n=1")
  first_tensor = ITensor(elt, sites[1], left_link)
  first_tensor[sites[1] => 1, left_link => 1] = amplitude
  placed!(first_tensor, (left_link => 2,), sites[1], convert(elt, coeffs[1]))
  tensors = ITensor[first_tensor]

  for j in 2:(nsites - 1)
    right_link = Index(2, "OperatorStateLink,n=$(j)")
    tensor = ITensor(elt, dag(left_link), sites[j], right_link)
    tensor[dag(left_link) => 1, sites[j] => 1, right_link => 1] = amplitude
    tensor[dag(left_link) => 2, sites[j] => 1, right_link => 2] = amplitude
    placed!(tensor, (dag(left_link) => 1, right_link => 2), sites[j], convert(elt, coeffs[j]))
    push!(tensors, tensor)
    left_link = right_link
  end

  last_tensor = ITensor(elt, dag(left_link), sites[nsites])
  last_tensor[dag(left_link) => 2, sites[nsites] => 1] = amplitude
  placed!(last_tensor, (dag(left_link) => 1,), sites[nsites], convert(elt, coeffs[nsites]))
  push!(tensors, last_tensor)
  return MPS(tensors)
end

"""
    operator_local_sum_state(sites, op, coeffs)

Build the operator-space `MPS` representing the literal local-density sum `sum_j c_j O_j`,
where `O_j` is the dense local operator `op` placed on site `j` and the literal `d x d`
identity `I` elsewhere.

# Arguments
- `sites`: Operator-space site indices.
- `op`: Dense `d x d` local operator.
- `coeffs`: One coefficient `c_j` per site.

# Returns
- An `MPS` of bond dimension 2 (link state 1 = "operator not yet placed", state 2 = "placed,
  identity henceforth") that equals `sum_j c_j O_j` exactly, including its trace/identity
  component. Because it uses the same literal-operator convention as
  [`operator_product_state`](@ref), the two compose correctly: e.g. `add(operator_product_state(
  sites, fill(Matrix{ComplexF64}(I, d, d), length(sites))), operator_local_sum_state(sites, sz,
  coeffs))` is the literal `I + sum_j c_j sz_j`.

# Notes
- This is the transport source builder: uniform `coeffs` give a total charge such as
  `sum_j S^z_j`, and a sign-split profile gives the domain-wall operator whose melting sets the
  dynamical exponent. `pauli_total_sz_state` and `pauli_domain_wall_state` are the `d = 2`
  cases, but keep their own pre-existing normalization rather than this function's
  literal-sum convention (see their docstrings for the exact factor).
"""
function operator_local_sum_state(sites, op::AbstractMatrix, coeffs::AbstractVector)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("operator_local_sum_state requires at least one site"))
  length(coeffs) == nsites || throw(ArgumentError("operator_local_sum_state needs one coefficient per site"))
  d = local_dimension(first(sites))
  all(local_dimension(site) == d for site in sites) ||
    throw(ArgumentError("operator_local_sum_state requires a uniform local dimension"))
  weights = _operator_coefficients(op, d)
  return _local_sum_state(sites, weights, coeffs, sqrt(d))
end

"""
    pauli_siteinds(nsites; tagprefix="PauliSpace")

Spin-1/2 case of [`operator_siteinds`](@ref): construct a vector of dimension-4 site indices
suitable for vectorized local Pauli bases, in the `(I, X, Y, Z)` ordering.

# Arguments
- `nsites`: Number of operator-space sites.

# Keyword Arguments
- `tagprefix`: Prefix used when naming the generated `Index` tags.

# Returns
- A vector of length `nsites` containing dimension-4 `Index` objects.
"""
function pauli_siteinds(nsites::Integer; tagprefix::AbstractString="PauliSpace")
  return operator_siteinds(nsites; d=2, tagprefix=tagprefix)
end

"""
    pauli_basis_state(sites, labels; coefficient=1.0)

Spin-1/2 case of [`operator_basis_state`](@ref): build a Pauli-basis product `MPS` in the
default local ordering `(I, X, Y, Z)`. Each entry of `labels` may be an integer basis index or
one of `"I"`, `"X"`, `"Y"`, `"Z"` (or the corresponding symbols).

# Arguments
- `sites`: Pauli-space site indices, typically created by [`pauli_siteinds`](@ref).
- `labels`: One local basis label per site.

# Keyword Arguments
- `coefficient`: Overall scalar prefactor stored on the first tensor.

# Returns
- A product-state `MPS` in operator space.

# Examples
```jldoctest
julia> sites = pauli_siteinds(3);

julia> rho = pauli_basis_state(sites, [:I, :Z, :I]);

julia> length(rho)
3
```
"""
function pauli_basis_state(sites, labels::AbstractVector; coefficient::Number=1.0)
  return operator_basis_state(sites, labels; coefficient=coefficient)
end

"""
    pauli_total_sz_state(sites; coefficient=nothing)

Spin-1/2 case of [`operator_local_sum_state`](@ref)'s tensor-network construction, for the
uniform total-`S^z` operator `sum_j S_j^z`.

# Arguments
- `sites`: Pauli-space site indices.

# Keyword Arguments
- `coefficient`: Optional scalar coefficient inserted directly for each single-`Z` string in
  the normalized Pauli-basis `MPS`. Defaults to `2^(N / 2 - 1)` for `N` sites, matching the
  normalized Pauli convention.

# Returns
- An `MPS` representation of the summed operator in Pauli space.

# Notes
- The returned state is not a simple product state; it uses bond dimension `2` to encode the
  operator sum compactly.
- Unlike [`operator_local_sum_state`](@ref) (which builds the literal `sum_j c_j O_j`), this
  wrapper keeps its pre-existing bond-dimension-2 amplitudes: it inserts amplitude `1` (not
  `sqrt(2)`) on unoccupied sites, so the returned operator carries an extra `(1 / sqrt(2))^N`
  normalization relative to the literal `sum_j coefficient * sigma^z_j`. That factor cancels in
  the ratio observables (e.g. [`pauli_expectation_profile`](@ref)) this state is built for; it
  was always present and is documented here for clarity.
- Requires `sites` of local dimension `2` (throws `ArgumentError` otherwise); use
  [`operator_local_sum_state`](@ref) directly for any other local dimension.
"""
function pauli_total_sz_state(sites; coefficient::Union{Nothing,Number}=nothing)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("pauli_total_sz_state requires at least one site"))
  d = local_dimension(first(sites))
  d == 2 || throw(ArgumentError("pauli_total_sz_state is the spin-1/2 case; use operator_local_sum_state for local dimension $(d)"))
  all(local_dimension(site) == d for site in sites) ||
    throw(ArgumentError("pauli_total_sz_state requires a uniform local dimension"))
  weight = isnothing(coefficient) ? 2.0^(nsites / 2 - 1) : coefficient
  weights = _operator_coefficients(pauli_matrices().Z / sqrt(2), d)
  return _local_sum_state(sites, weights, fill(weight, nsites), 1.0)
end

"""
    pauli_domain_wall_state(sites; kink=length(sites) ÷ 2, coefficient=1.0)

Spin-1/2 case of [`operator_local_sum_state`](@ref)'s tensor-network construction, for the
signed domain-wall operator `D = sum_j sign(j - kink) * sigma^z_j`, the infinite-temperature
transport source.

# Arguments
- `sites`: Pauli-space site indices, typically from [`pauli_siteinds`](@ref).

# Keyword Arguments
- `kink`: Site index after which the sign flips from `-` to `+`. The default `length(sites) ÷ 2`
  places the wall on the central bond (perfectly antisymmetric for even `length(sites)`). Must
  lie in `0:length(sites)`.
- `coefficient`: Magnitude of each single-`Z` coefficient. The overall scale is arbitrary for
  ratio observables (it cancels in a normalized profile).

# Returns
- An `MPS` representation of the operator in Pauli space.

# Notes
- Like [`pauli_total_sz_state`](@ref) the returned state is not a product state; it uses bond
  dimension `2` to encode the operator sum compactly (the running-sum link state `1` = "single
  `Z` not yet placed", state `2` = "placed, identity henceforth"). The per-site sign is folded
  into the `Z`-placement amplitude, so the bond dimension is independent of `kink`.
- `D` is traceless, so measure its profile with `normalize=false` (see
  [`pauli_expectation_profile`](@ref)).
- Unlike [`operator_local_sum_state`](@ref) (which builds the literal `sum_j c_j O_j`), this
  wrapper keeps its pre-existing bond-dimension-2 amplitudes: it inserts amplitude `1` (not
  `sqrt(2)`) on unoccupied sites, so the returned operator carries an extra `(1 / sqrt(2))^N`
  normalization relative to the literal `sum_j sign(j - kink) * coefficient * sigma^z_j`. That
  factor cancels in ratio observables; it was always present and is documented here for
  clarity.
- Requires `sites` of local dimension `2` (throws `ArgumentError` otherwise); use
  [`operator_local_sum_state`](@ref) directly for any other local dimension.

# Examples
```jldoctest
julia> sites = pauli_siteinds(4);

julia> dw = pauli_domain_wall_state(sites);

julia> length(dw)
4
```
"""
function pauli_domain_wall_state(sites; kink::Integer=length(sites) ÷ 2, coefficient::Number=1.0)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("pauli_domain_wall_state requires at least one site"))
  0 <= kink <= nsites || throw(ArgumentError("kink must lie in 0:length(sites)"))
  d = local_dimension(first(sites))
  d == 2 || throw(ArgumentError("pauli_domain_wall_state is the spin-1/2 case; use operator_local_sum_state for local dimension $(d)"))
  all(local_dimension(site) == d for site in sites) ||
    throw(ArgumentError("pauli_domain_wall_state requires a uniform local dimension"))
  weights = _operator_coefficients(pauli_matrices().Z / sqrt(2), d)
  coeffs = [j <= kink ? -coefficient : coefficient for j in 1:nsites]
  return _local_sum_state(sites, weights, coeffs, 1.0)
end
