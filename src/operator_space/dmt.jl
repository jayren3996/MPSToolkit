"""
    DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=480, connector_buffer=8)

Options controlling operator-space density matrix truncation (DMT).

!!! warning "Transport-specific algorithm"
    DMT is a specialized truncation scheme designed for **transport** (e.g. spin or energy
    diffusion) in operator space.  It protects local reduced operator data, including the
    identity/trace component and nearby Pauli components, before truncating connected
    long-range correlations. For general operator-space evolution without this transport
    bias, use ordinary TEBD truncation (`LocalGateEvolution`) instead.

# Fields
- `maxdim`: Target bond dimension after DMT truncation.
- `cutoff`: Truncation cutoff used in the repair SVD.
- `gate_maxdim`: Temporary bond dimension budget used during raw gate application.
- `connector_buffer`: Number of connector directions protected before reduced-matrix
  truncation.
"""
struct DMTOptions
  maxdim::Int
  cutoff::Float64
  gate_maxdim::Int
  connector_buffer::Int

  function DMTOptions(maxdim, cutoff, gate_maxdim, connector_buffer)
    maxdim >= 1 || throw(ArgumentError("DMTOptions requires maxdim >= 1"))
    cutoff >= 0 || throw(ArgumentError("DMTOptions requires cutoff >= 0"))
    gate_maxdim >= 1 || throw(ArgumentError("DMTOptions requires gate_maxdim >= 1"))
    connector_buffer >= 0 || throw(ArgumentError("DMTOptions requires connector_buffer >= 0"))
    connector_buffer <= maxdim || throw(ArgumentError("DMTOptions requires connector_buffer <= maxdim"))
    return new(Int(maxdim), Float64(cutoff), Int(gate_maxdim), Int(connector_buffer))
  end

  function DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=480, connector_buffer=8)
    return DMTOptions(maxdim, cutoff, gate_maxdim, connector_buffer)
  end
end

"""
    _pauli_identity_env(site)

Return the local Pauli-basis identity ket used when building DMT environments.
"""
function _pauli_identity_env(site)
  tensor = ITensor(site)
  tensor[site => 1] = 1.0
  return tensor
end

"""
    _left_identity_environment(psi, stop)

Contract the left identity environment used by DMT truncation up to site `stop`.
"""
function _left_identity_environment(psi::MPS, stop::Integer)
  env = ITensor(1.0)
  for site in 1:Int(stop)
    env *= _pauli_identity_env(siteind(psi, site)) * psi[site]
  end
  return env
end

"""
    _right_identity_environment(psi, start)

Contract the right identity environment used by DMT truncation starting at site `start`.
"""
function _right_identity_environment(psi::MPS, start::Integer)
  env = ITensor(1.0)
  for site in length(psi):-1:Int(start)
    env *= _pauli_identity_env(siteind(psi, site)) * psi[site]
  end
  return env
end

"""
    _dmt_truncation_bonds(start, span, direction)

Return the internal bonds at which DMT truncation should be applied for one local update.
"""
function _dmt_truncation_bonds(start::Integer, span::Integer, direction::Symbol)
  bonds = collect(Int(start):(Int(start) + Int(span) - 2))
  direction === :R && return bonds
  direction === :L && return reverse(bonds)
  throw(ArgumentError("DMT direction must be :R or :L"))
end

function _validate_pauli_operator_space(psi::MPS, start::Integer, span::Integer)
  for site in start:(start + span - 1)
    dim(siteind(psi, site)) == 4 || throw(ArgumentError("DMT assumes Pauli operator-space sites ordered as (I, X, Y, Z) with local dimension 4"))
  end
  return nothing
end

function _validate_dmt_step(psi::MPS, gate::AbstractMatrix, start::Integer, span::Integer, direction::Symbol, maxdim::Integer, connector_buffer::Integer)
  direction === :R || direction === :L || throw(ArgumentError("DMT direction must be :R or :L"))
  maxdim >= 0 || throw(ArgumentError("DMT maxdim must be nonnegative"))
  connector_buffer >= 0 || throw(ArgumentError("DMT connector_buffer must be nonnegative"))
  maxdim == 0 || connector_buffer <= maxdim || throw(ArgumentError("DMT connector_buffer must be <= maxdim"))
  start >= 1 || throw(ArgumentError("local gate bond must be at least 1"))
  last_site = start + span - 1
  last_site <= length(psi) || throw(ArgumentError("local gate support exceeds chain length"))
  span == 1 && return nothing
  start == length(psi) && throw(ArgumentError("periodic boundary DMT is not implemented for local gates"))
  for bond in start:(last_site - 1)
    1 <= bond < length(psi) || throw(ArgumentError("DMT target bonds must lie in 1:length(psi)-1"))
  end
  _validate_pauli_operator_space(psi, start, span)
  return nothing
end

"""
    _mat_trunc!(matrix_data, χ; connector_buffer=8)

Apply the reduced-matrix truncation step used by DMT.

# Arguments
- `matrix_data`: Dense reduced matrix to truncate in place.
- `χ`: Number of singular directions to retain after removing the protected connector block.

# Keyword Arguments
- `connector_buffer`: Size of the connector block left untouched at the top-left corner.

# Returns
- `nothing`. The input matrix is modified in place.
"""
function _mat_trunc!(matrix_data::AbstractMatrix, χ::Integer; connector_buffer::Integer=8)
  size(matrix_data, 1) < connector_buffer + 1 && return nothing
  χ + connector_buffer >= size(matrix_data, 1) && return nothing
  matrix_data[1, 1] == 0 && return nothing

  connector = (matrix_data[:, 1:1] * matrix_data[1:1, :]) / matrix_data[1, 1]
  matrix_data .-= connector

  trailing = (connector_buffer + 1):size(matrix_data, 1)
  factorization = svd(matrix_data[trailing, trailing])
  retained = min(χ, length(factorization.S))
  matrix_data[trailing, trailing] .= factorization.U[:, 1:retained] *
                                     Diagonal(factorization.S[1:retained]) *
                                     factorization.Vt[1:retained, :]
  matrix_data .+= connector
  return nothing
end

function _complete_orthonormal_basis(protected::AbstractMatrix, target_dim::Integer=size(protected, 1))
  ambient_dim = size(protected, 1)
  0 <= target_dim <= ambient_dim || throw(ArgumentError("orthonormal basis target dimension must lie in 0:size(protected, 1)"))
  target_dim == 0 && return zeros(eltype(protected), ambient_dim, 0)

  basis = Matrix{eltype(protected)}(undef, ambient_dim, 0)
  if size(protected, 2) > 0
    factorization = svd(Matrix(protected))
    scale = isempty(factorization.S) ? zero(real(float(one(eltype(factorization.S))))) : maximum(factorization.S)
    tolerance = max(ambient_dim, size(protected, 2)) * eps(real(float(scale == 0 ? one(scale) : scale))) * max(scale, one(scale))
    protected_rank = min(count(>(tolerance), factorization.S), target_dim)
    protected_rank > 0 && (basis = factorization.U[:, 1:protected_rank])
  end

  for column in 1:ambient_dim
    size(basis, 2) == target_dim && break
    candidate = zeros(eltype(protected), ambient_dim)
    candidate[column] = one(eltype(protected))
    for basis_column in 1:size(basis, 2)
      candidate .-= basis[:, basis_column] * (basis[:, basis_column]' * candidate)
    end
    candidate_norm = norm(candidate)
    if candidate_norm > ambient_dim * eps(real(float(candidate_norm == 0 ? one(candidate_norm) : candidate_norm)))
      basis = hcat(basis, candidate / candidate_norm)
    end
  end

  size(basis, 2) == target_dim || throw(ArgumentError("could not complete orthonormal DMT basis"))
  return basis
end

"""
    _dmt_bond_truncate!(psi, bond; maxdim, cutoff, direction=:R, connector_buffer=8)

Perform one DMT-preserving bond truncation step.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `bond`: Bond index to truncate.

# Keyword Arguments
- `maxdim`: Target bond dimension.
- `cutoff`: Truncation cutoff used in the final repair SVD.
- `direction`: Sweep direction, either `:R` or `:L`.
- `connector_buffer`: Number of protected connector directions.

# Returns
- The mutated `psi`.
"""
function _dmt_bond_truncate!(
  psi::MPS,
  bond::Integer;
  maxdim::Integer,
  cutoff::Real,
  direction::Symbol=:R,
  connector_buffer::Integer=8,
  left_env=nothing,
  right_env=nothing,
  orthogonalize::Bool=true,
)
  maxdim > 0 || return psi
  connector_buffer >= 0 || throw(ArgumentError("DMT connector_buffer must be nonnegative"))
  connector_buffer <= maxdim || throw(ArgumentError("DMT connector_buffer must be <= maxdim"))
  1 <= bond < length(psi) || throw(ArgumentError("DMT bond must lie in 1:length(psi)-1"))
  current_link = linkind(psi, bond)
  isnothing(current_link) && return psi
  dim(current_link) <= maxdim && return psi

  orthogonalize && orthogonalize!(psi, bond)
  left_site = siteind(psi, bond)
  right_site = siteind(psi, bond + 1)

  isnothing(left_env) && (left_env = _left_identity_environment(psi, bond - 1))
  isnothing(right_env) && (right_env = _right_identity_environment(psi, bond + 2))

  previous_link = linkind(psi, bond - 1)
  left_inds = isnothing(previous_link) ? (left_site,) : (previous_link, left_site)
  u, s, v = svd(psi[bond], left_inds)
  psi[bond] = u
  psi[bond + 1] = v * psi[bond + 1]

  left_link = commonind(u, s)
  right_link = commonind(v, s)
  left_basis = _complete_orthonormal_basis(matrix(left_env * psi[bond], left_link, left_site), dim(left_link))
  right_basis = _complete_orthonormal_basis(matrix(psi[bond + 1] * right_env, right_link, right_site), dim(right_link))
  singular_values = Matrix(matrix(s, left_link, right_link))

  reduced = left_basis' * singular_values * right_basis
  _mat_trunc!(reduced, maxdim - connector_buffer; connector_buffer=connector_buffer)

  repaired = ITensor(left_basis * reduced * right_basis', left_link, right_link)
  new_u, new_s, new_v = svd(repaired, left_link; maxdim=maxdim, cutoff=cutoff)
  if direction === :R
    psi[bond] *= new_u
    psi[bond + 1] = new_s * new_v * psi[bond + 1]
  else
    psi[bond] = psi[bond] * new_u * new_s
    psi[bond + 1] = new_v * psi[bond + 1]
  end
  return psi
end

function _left_identity_prefixes(psi::MPS)
  prefixes = Vector{ITensor}(undef, length(psi) + 1)
  prefixes[1] = ITensor(1.0)
  for site in 1:length(psi)
    prefixes[site + 1] = prefixes[site] * _pauli_identity_env(siteind(psi, site)) * psi[site]
  end
  return prefixes
end

function _right_identity_suffixes(psi::MPS)
  suffixes = Vector{ITensor}(undef, length(psi) + 2)
  suffixes[length(psi) + 1] = ITensor(1.0)
  suffixes[length(psi) + 2] = ITensor(1.0)
  for site in length(psi):-1:1
    suffixes[site] = suffixes[site + 1] * _pauli_identity_env(siteind(psi, site)) * psi[site]
  end
  return suffixes
end

function _dmt_window_truncate!(psi::MPS, start::Integer, span::Integer; maxdim::Integer, cutoff::Real, direction::Symbol, connector_buffer::Integer)
  span <= 1 && return psi
  bonds = _dmt_truncation_bonds(start, span, direction)
  isempty(bonds) && return psi

  orthogonalize!(psi, first(bonds))
  if direction === :R
    right_suffixes = _right_identity_suffixes(psi)
    left_env = _left_identity_environment(psi, first(bonds) - 1)
    for (index, bond) in pairs(bonds)
      _dmt_bond_truncate!(
        psi,
        bond;
        maxdim=maxdim,
        cutoff=cutoff,
        direction=direction,
        connector_buffer=connector_buffer,
        left_env=left_env,
        right_env=right_suffixes[bond + 2],
        orthogonalize=false,
      )
      index < length(bonds) && (left_env = left_env * _pauli_identity_env(siteind(psi, bond)) * psi[bond])
    end
  else
    left_prefixes = _left_identity_prefixes(psi)
    right_env = _right_identity_environment(psi, first(bonds) + 2)
    for (index, bond) in pairs(bonds)
      _dmt_bond_truncate!(
        psi,
        bond;
        maxdim=maxdim,
        cutoff=cutoff,
        direction=direction,
        connector_buffer=connector_buffer,
        left_env=left_prefixes[bond],
        right_env=right_env,
        orthogonalize=false,
      )
      index < length(bonds) && (right_env = right_env * _pauli_identity_env(siteind(psi, bond + 1)) * psi[bond + 1])
    end
  end
  return psi
end

"""
    dmt_step!(psi, gate, bond; maxdim=30, cutoff=1e-12, direction=:R, gate_maxdim=max(maxdim * 16, 64), connector_buffer=8)

Apply one local operator-space gate and then perform DMT-preserving truncation.

This is a **transport-specific** truncation step.  See [`DMTOptions`](@ref) for when DMT
is (and is not) the appropriate choice.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `gate`: Dense local gate in the Pauli basis.
- `bond`: Left-edge location of the local update.

# Keyword Arguments
- `maxdim`: Target post-truncation bond dimension.
- `cutoff`: Truncation cutoff used in the final repair SVD.
- `direction`: Sweep direction, either `:R` or `:L`.
- `gate_maxdim`: Temporary bond dimension budget used for the raw gate application.
- `connector_buffer`: Number of protected connector directions.

# Returns
- The mutated `psi`.
"""
function dmt_step!(
  psi::MPS,
  gate::AbstractMatrix,
  bond;
  maxdim::Integer=30,
  cutoff::Real=1e-12,
  direction::Symbol=:R,
  gate_maxdim::Integer=max(Int(maxdim) * 16, 64),
  connector_buffer::Integer=8,
)
  start = _bond_start(bond)
  span = _operator_span_at(psi, gate, start)
  _validate_dmt_step(psi, gate, start, span, direction, Int(maxdim), Int(connector_buffer))
  tebd_evolve!(psi, gate, start; maxdim=Int(gate_maxdim), cutoff=0.0)
  _dmt_window_truncate!(
    psi,
    start,
    span;
    maxdim=Int(maxdim),
    cutoff=cutoff,
    direction=direction,
    connector_buffer=Int(connector_buffer),
  )
  return psi
end

"""
    _reverse_gate_for_step(gate_spec, schedule, reverse_schedule, bond, index)

Resolve the gate for a reverse DMT schedule entry. Matrix and callable gate providers keep the
same semantics as forward sweeps. Vector gate providers are mapped back to the corresponding
forward schedule entry.
"""
function _is_default_reverse_schedule(schedule, reverse_schedule)
  length(reverse_schedule) == length(schedule) || return false
  for index in eachindex(reverse_schedule)
    reverse_schedule[index] == schedule[length(schedule) - index + 1] || return false
  end
  return true
end

function _reverse_gate_index(schedule, reverse_schedule, bond, index)
  if _is_default_reverse_schedule(schedule, reverse_schedule)
    return length(schedule) - index + 1
  end

  count(==(bond), schedule) == 1 || throw(ArgumentError("custom reverse DMT schedules with repeated bonds require a callable gate provider that does not depend on reverse indices"))
  forward_index = findfirst(==(bond), schedule)
  isnothing(forward_index) && throw(ArgumentError("reverse DMT schedule contains bond $(bond), which is absent from the forward schedule"))
  return forward_index
end

function _reverse_gate_for_step(gate_spec::AbstractVector, schedule, reverse_schedule, bond, index)
  return _gate_for_step(gate_spec, bond, _reverse_gate_index(schedule, reverse_schedule, bond, index))
end

function _reverse_gate_for_step(gate_spec::Function, schedule, reverse_schedule, bond, index)
  return _gate_for_step(gate_spec, bond, _reverse_gate_index(schedule, reverse_schedule, bond, index))
end

function _reverse_gate_for_step(gate_spec, schedule, reverse_schedule, bond, index)
  return _gate_for_step(gate_spec, bond, index)
end

"""
    dmt_evolve!(psi, evo::DMTGateEvolution)

Run scheduled operator-space DMT evolution.

This driver is intended for **transport simulations** (e.g. spin or energy diffusion).  See
[`DMTOptions`](@ref) for the transport-specific assumptions built into DMT truncation.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `evo`: [`DMTGateEvolution`](@ref) describing the gate specification, schedules, and
  truncation budgets.

# Returns
- The mutated and normalized `psi`.

# Notes
- One call runs `evo.nstep` complete forward-and-reverse sweeps.
- `normalize!(psi)` is applied at the end because operator-space trajectories typically
  interpret the state as a normalized vectorized operator.
"""
function dmt_evolve!(psi::MPS, evo::DMTGateEvolution)
  for _ in 1:evo.nstep
    for (index, bond) in pairs(evo.schedule)
      local_gate = _gate_for_step(evo.gate, bond, index)
      dmt_step!(
        psi,
        local_gate,
        bond;
        maxdim=evo.maxdim,
        cutoff=evo.cutoff,
        direction=:R,
        gate_maxdim=evo.gate_maxdim,
        connector_buffer=evo.connector_buffer,
      )
    end
    for (index, bond) in pairs(evo.reverse_schedule)
      local_gate = _reverse_gate_for_step(evo.gate, evo.schedule, evo.reverse_schedule, bond, index)
      dmt_step!(
        psi,
        local_gate,
        bond;
        maxdim=evo.maxdim,
        cutoff=evo.cutoff,
        direction=:L,
        gate_maxdim=evo.gate_maxdim,
        connector_buffer=evo.connector_buffer,
      )
    end
  end
  normalize!(psi)
  return psi
end

"""
    evolve!(psi, evo::DMTGateEvolution)

Dispatch operator-space evolution through the DMT backend.
"""
function evolve!(psi::MPS, evo::DMTGateEvolution)
  return dmt_evolve!(psi, evo)
end
