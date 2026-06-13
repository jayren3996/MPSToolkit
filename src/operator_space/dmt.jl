"""
    DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=max(maxdim * 16, 64), connector_buffer=8)

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

  function DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=max(maxdim * 16, 64), connector_buffer=8)
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
  maxdim >= 1 || throw(ArgumentError("DMT maxdim must be >= 1"))
  connector_buffer >= 0 || throw(ArgumentError("DMT connector_buffer must be nonnegative"))
  connector_buffer <= maxdim || throw(ArgumentError("DMT connector_buffer must be <= maxdim"))
  start >= 1 || throw(ArgumentError("local gate bond must be at least 1"))
  last_site = start + span - 1
  last_site <= length(psi) || throw(ArgumentError("local gate support exceeds chain length"))
  # Validate Pauli operator-space dimensions for every span, including span == 1, so a
  # single-site gate on a non-Pauli (non-dimension-4) site is rejected rather than silently
  # accepted.
  _validate_pauli_operator_space(psi, start, span)
  span == 1 && return nothing
  start == length(psi) && throw(ArgumentError("periodic boundary DMT is not implemented for local gates"))
  for bond in start:(last_site - 1)
    1 <= bond < length(psi) || throw(ArgumentError("DMT target bonds must lie in 1:length(psi)-1"))
  end
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

  # Column/row 1 is the rank-1 identity/trace connector. Subtract it before truncating so the
  # identity direction is preserved exactly -- but only when it is well-conditioned. A
  # traceless operator (e.g. a transport current in operator space) has (near) zero identity
  # overlap, where `matrix_data[1, 1] ≈ 0`; an exact `== 0` guard would then skip truncation
  # entirely (leaving `maxdim` unenforced) and a tiny-but-nonzero value would blow up the
  # connector via the `1 / matrix_data[1, 1]` scaling. Use a relative tolerance and, when the
  # connector is negligible, skip subtracting it but still truncate so `maxdim` is enforced.
  scale = norm(matrix_data)
  tolerance = size(matrix_data, 1) * eps(real(eltype(matrix_data))) * scale
  has_connector = abs(matrix_data[1, 1]) > tolerance
  if has_connector
    connector = (matrix_data[:, 1:1] * matrix_data[1:1, :]) / matrix_data[1, 1]
    matrix_data .-= connector
  end

  trailing = (connector_buffer + 1):size(matrix_data, 1)
  factorization = svd(matrix_data[trailing, trailing])
  retained = min(χ, length(factorization.S))
  matrix_data[trailing, trailing] .= factorization.U[:, 1:retained] *
                                     Diagonal(factorization.S[1:retained]) *
                                     factorization.Vt[1:retained, :]

  if has_connector
    matrix_data .+= connector
  end
  return nothing
end

function _complete_orthonormal_basis(protected::AbstractMatrix, target_dim::Integer=size(protected, 1))
  ambient_dim = size(protected, 1)
  0 <= target_dim <= ambient_dim || throw(ArgumentError("orthonormal basis target dimension must lie in 0:size(protected, 1)"))
  element_type = eltype(protected)
  target_dim == 0 && return zeros(element_type, ambient_dim, 0)

  # Preallocate the whole basis and fill it column by column in place. The leading columns are
  # the trace/connector-aligned singular vectors of `protected` (in singular-value order, so
  # column 1 is the dominant connector direction `_mat_trunc!` protects); the remaining columns
  # complete the space from the standard basis. The column selection and the leading protected
  # block are identical to the previous per-column Gram-Schmidt sweep, but the orthogonalization
  # is a single BLAS matrix-vector product against the already-filled block instead of a scalar
  # loop, and the basis grows in place rather than via repeated `hcat`. That removes the O(d^3)
  # allocation traffic (`hcat` reallocating the basis and `basis[:, j]` copies every step) that
  # dominated DMT evolution runtime.
  basis = Matrix{element_type}(undef, ambient_dim, target_dim)
  filled = 0
  if size(protected, 2) > 0
    factorization = svd(Matrix(protected))
    scale = isempty(factorization.S) ? zero(real(float(one(element_type)))) : maximum(factorization.S)
    tolerance = max(ambient_dim, size(protected, 2)) * eps(real(float(scale == 0 ? one(scale) : scale))) * max(scale, one(scale))
    protected_rank = min(count(>(tolerance), factorization.S), target_dim)
    if protected_rank > 0
      @views basis[:, 1:protected_rank] .= factorization.U[:, 1:protected_rank]
      filled = protected_rank
    end
  end

  candidate = Vector{element_type}(undef, ambient_dim)
  projection = Vector{element_type}(undef, target_dim)
  column = 0
  while filled < target_dim && column < ambient_dim
    column += 1
    fill!(candidate, zero(element_type))
    candidate[column] = one(element_type)
    if filled > 0
      block = @view basis[:, 1:filled]
      coefficients = @view projection[1:filled]
      mul!(coefficients, block', candidate)                                       # coefficients = block' * e_column
      mul!(candidate, block, coefficients, -one(element_type), one(element_type)) # candidate -= block * coefficients
    end
    candidate_norm = norm(candidate)
    # Linear-dependence test relative to the candidate's *original* unit norm, not the
    # post-projection residual. Each candidate starts as a unit standard-basis vector, so a
    # genuinely new direction leaves an O(1) residual while a dependent one leaves only
    # O(ambient * eps) roundoff. Scaling the threshold by the residual's own magnitude (as a
    # naive `eps(candidate_norm)` would) collapses it to ~eps^2 for roundoff residuals, which
    # then pass as spurious near-duplicate columns and break orthonormality.
    if candidate_norm > ambient_dim * eps(one(real(float(candidate_norm))))
      filled += 1
      @views basis[:, filled] .= candidate ./ candidate_norm
    end
  end

  filled == target_dim || throw(ArgumentError("could not complete orthonormal DMT basis"))
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
- `orthogonalize`: If `true` (default), re-gauge the MPS so the orthogonality center is at
  `bond` before truncating, as the connector-preserving construction requires. Set to `false`
  only when the center is already known to be at `bond`; otherwise the bond SVD no longer
  returns the true Schmidt values and the truncation is invalid.
- `left_env`, `right_env`: Optional precomputed identity/trace environments. When `nothing`
  (default) they are built from the current canonical gauge. If supplied, they MUST match the
  gauge at `bond` (e.g. as produced by `orthogonalize!(psi, bond)`); a stale environment
  silently yields an invalid truncation.

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
  # `s` is the diagonal Schmidt matrix from the bond SVD. Fold it into the basis change as a
  # column scaling (Diagonal) rather than materializing a dense matrix and doing two full
  # matmuls -- one O(n^3) product per bond becomes an O(n^2) scaling.
  singular_values = diag(matrix(s, left_link, right_link))

  reduced = (left_basis' * Diagonal(singular_values)) * right_basis
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

function _dmt_window_truncate!(psi::MPS, start::Integer, span::Integer; maxdim::Integer, cutoff::Real, direction::Symbol, connector_buffer::Integer)
  span <= 1 && return psi

  # Truncate every bond inside the gate window as an independent single-bond DMT update, with
  # the orthogonality center restored to that bond before truncating it. The
  # connector-preserving construction in `_dmt_bond_truncate!` assumes a *canonical* gauge: the
  # bond SVD must return the true Schmidt values and the trace environments must be built from
  # orthonormal blocks. A multi-bond window cannot satisfy that with a single shared sweep,
  # because each DMT truncation rewrites the bond tensors and resets the MPS gauge bookkeeping.
  # Sweeping with cached environments therefore truncates later bonds in a non-canonical gauge:
  # mildly wrong for a `:R` window, and badly wrong for `:L`, where the orthogonality center
  # ends up on the wrong site so `svd(psi[bond])` degenerates to `s ≈ I` and the truncation
  # discards information indiscriminately. Re-gauging per bond keeps a span-`S` window exactly
  # equal to the verified single-bond path applied to each of its bonds.
  for bond in _dmt_truncation_bonds(start, span, direction)
    _dmt_bond_truncate!(
      psi,
      bond;
      maxdim=maxdim,
      cutoff=cutoff,
      direction=direction,
      connector_buffer=connector_buffer,
      orthogonalize=true,
    )
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
- `gate_maxdim`: Temporary bond dimension budget used for the raw gate application before DMT
  truncates the bond back to `maxdim`. A large `gate_maxdim` lets the gate inflate the bond
  prior to truncation, which can reduce accuracy at small `connector_buffer`; choose it
  together with `maxdim` and `connector_buffer` rather than independently.
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
    dmt_step!(psi, gate, bond, opts::DMTOptions; direction=:R)

Apply one DMT step using the truncation settings bundled in a [`DMTOptions`](@ref). This is a
convenience overload of [`dmt_step!`](@ref) that forwards `opts` fields to the keyword form.

# Returns
- The mutated `psi`.
"""
function dmt_step!(psi::MPS, gate::AbstractMatrix, bond, opts::DMTOptions; direction::Symbol=:R)
  return dmt_step!(
    psi,
    gate,
    bond;
    maxdim=opts.maxdim,
    cutoff=opts.cutoff,
    direction=direction,
    gate_maxdim=opts.gate_maxdim,
    connector_buffer=opts.connector_buffer,
  )
end

"""
    _is_default_reverse_schedule(schedule, reverse_schedule)

Return `true` when `reverse_schedule` is exactly `reverse(schedule)`.
"""
function _is_default_reverse_schedule(schedule, reverse_schedule)
  length(reverse_schedule) == length(schedule) || return false
  for index in eachindex(reverse_schedule)
    reverse_schedule[index] == schedule[length(schedule) - index + 1] || return false
  end
  return true
end

function _reverse_gate_index(is_default, schedule, bond, index)
  if is_default
    return length(schedule) - index + 1
  end

  count(==(bond), schedule) == 1 || throw(ArgumentError("custom reverse DMT schedules with repeated bonds require a callable gate provider that does not depend on reverse indices"))
  forward_index = findfirst(==(bond), schedule)
  isnothing(forward_index) && throw(ArgumentError("reverse DMT schedule contains bond $(bond), which is absent from the forward schedule"))
  return forward_index
end

"""
    _reverse_gate_for_step(gate_spec, is_default, schedule, bond, index)

Resolve the gate for a reverse DMT schedule entry. Matrix and callable gate providers keep the
same semantics as forward sweeps. Vector gate providers are mapped back to the corresponding
forward schedule entry. `is_default` is the precomputed
`_is_default_reverse_schedule(schedule, reverse_schedule)` flag, hoisted out of the per-step
loop by the caller so the loop-invariant schedule scan runs once per `dmt_evolve!`, not once
per reverse step.
"""
function _reverse_gate_for_step(gate_spec::AbstractVector, is_default, schedule, bond, index)
  return _gate_for_step(gate_spec, bond, _reverse_gate_index(is_default, schedule, bond, index))
end

function _reverse_gate_for_step(gate_spec::Function, is_default, schedule, bond, index)
  return _gate_for_step(gate_spec, bond, _reverse_gate_index(is_default, schedule, bond, index))
end

function _reverse_gate_for_step(gate_spec, is_default, schedule, bond, index)
  return _gate_for_step(gate_spec, bond, index)
end

"""
    dmt_evolve!(psi, evo::DMTGateEvolution; normalize=true)

Run scheduled operator-space DMT evolution.

This driver is intended for **transport simulations** (e.g. spin or energy diffusion).  See
[`DMTOptions`](@ref) for the transport-specific assumptions built into DMT truncation.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `evo`: [`DMTGateEvolution`](@ref) describing the gate specification, schedules, and
  truncation budgets.

# Keyword Arguments
- `normalize`: If `true` (default, taken from `evo.normalize`), `normalize!(psi)` at the end, the usual convention for
  trajectories that interpret the state as a normalized vectorized operator. Set `false`
  when tracking **unnormalized** traces of a traceless operator (e.g. the conserved
  `tr(H O(t))` of a Heisenberg-evolved energy density): truncation sheds Hilbert-Schmidt
  norm faster than it sheds the DMT-protected components, so per-step renormalization
  silently inflates absolute traces over time even though ratio observables are unaffected.

# Returns
- The mutated (and, by default, normalized) `psi`.

# Notes
- One call runs `evo.nstep` complete forward-and-reverse sweeps.
"""
function dmt_evolve!(psi::MPS, evo::DMTGateEvolution; normalize::Bool=evo.normalize)
  reverse_is_default = _is_default_reverse_schedule(evo.schedule, evo.reverse_schedule)
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
      local_gate = _reverse_gate_for_step(evo.gate, reverse_is_default, evo.schedule, bond, index)
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
  normalize && normalize!(psi)
  return psi
end

"""
    evolve!(psi, evo::DMTGateEvolution)

Dispatch operator-space evolution through the DMT backend.

The `normalize` keyword is forwarded to [`dmt_evolve!`](@ref); set it `false` when tracking
unnormalized traces of a traceless operator (e.g. conserved `tr(H O(t))`).
"""
function evolve!(psi::MPS, evo::DMTGateEvolution; normalize::Bool=evo.normalize)
  return dmt_evolve!(psi, evo; normalize=normalize)
end
