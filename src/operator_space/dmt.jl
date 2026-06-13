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
    _DMTEnvCache(psi)

Amortizing cache for the DMT identity/trace environments (perf-1). `_dmt_bond_truncate!` needs,
at each bond `b`, the left environment over sites `1:b-1` and the right environment over
`b+2:N`. Rebuilding both from scratch on every call is the O(N) per-bond cost that makes a
schedule sweep O(N^2) (see the note on `_dmt_window_truncate!`). This cache stores the prefix
(`left`) and suffix (`right`) contractions and rebuilds only the part invalidated since the last
use, so a *local* sweep -- the regime every gate schedule produces (verified: each gate's
footprint is its O(span) window, never the whole chain) -- costs O(1) amortized per bond.

Validity is tracked by two watermarks: `left[0..lvalid]` and `right[rvalid..N+1]` are current. A
mutation over sites `[lo,hi]` lowers `lvalid` to `lo-2` and raises `rvalid` to `hi+2`. The
one-site-beyond ("-2/+2") margin is required because a left env over `1:k` carries the OPEN index
`linkind(psi,k)`: a factorization at site `m` re-creates the adjacent link, so the env whose open
index touches a mutated bond is stale even when its contracted sites are not.

The cache memoizes *only* the environment rebuilds; `_dmt_bond_truncate!` keeps its exact
operation sequence (including `orthogonalize!`), so a cached env that equals the from-scratch
rebuild yields bit-for-bit identical truncation output. `_DMT_VERIFY_ENVS[] = true` asserts that
equality at every bond and is exercised by the ED-oracle test suite.
"""
mutable struct _DMTEnvCache
  left::Vector{ITensor}    # left[k+1]  = _left_identity_environment(psi, k),  k in 0:n
  right::Vector{ITensor}   # right[k]   = _right_identity_environment(psi, k), k in 1:n+1
  lvalid::Int              # left[0..lvalid] valid
  rvalid::Int              # right[rvalid..n+1] valid
  n::Int
end

function _DMTEnvCache(psi::MPS)
  n = length(psi)
  left = Vector{ITensor}(undef, n + 1)
  right = Vector{ITensor}(undef, n + 1)
  left[1] = ITensor(1.0)        # L[0]
  right[n + 1] = ITensor(1.0)   # R[n+1]
  return _DMTEnvCache(left, right, 0, n + 1, n)
end

# When true, every cached environment is checked against the from-scratch rebuild in
# `_dmt_bond_truncate!` (index identity, then norm of the difference). Off in production; the
# ED-oracle tests flip it on to prove the threaded path reproduces the rebuild bit-for-bit.
const _DMT_VERIFY_ENVS = Ref(false)

"""
    _left_env_at!(cache, psi, k)

Return `_left_identity_environment(psi, k)`, extending the cached prefix from the current
watermark by one site at a time. O(1) amortized when `k` advances locally.
"""
function _left_env_at!(cache::_DMTEnvCache, psi::MPS, k::Integer)
  k = Int(k)
  0 <= k <= cache.n || throw(ArgumentError("DMT left env index $k out of range 0:$(cache.n)"))
  while cache.lvalid < k
    j = cache.lvalid + 1
    cache.left[j + 1] = cache.left[j] * (_pauli_identity_env(siteind(psi, j)) * psi[j])
    cache.lvalid = j
  end
  return cache.left[k + 1]
end

"""
    _right_env_at!(cache, psi, k)

Return `_right_identity_environment(psi, k)`, extending the cached suffix from the current
watermark by one site at a time. O(1) amortized when `k` retreats locally.
"""
function _right_env_at!(cache::_DMTEnvCache, psi::MPS, k::Integer)
  k = Int(k)
  1 <= k <= cache.n + 1 || throw(ArgumentError("DMT right env index $k out of range 1:$(cache.n + 1)"))
  while cache.rvalid > k
    j = cache.rvalid - 1
    cache.right[j] = cache.right[j + 1] * (_pauli_identity_env(siteind(psi, j)) * psi[j])
    cache.rvalid = j
  end
  return cache.right[k]
end

"""
    _invalidate_env!(cache, lo, hi)

Mark the environments stale for a mutation that touched sites `lo:hi`. Lowers the left watermark
to `lo-2` and raises the right watermark to `hi+2` (the one-site-beyond margin covers the open
link index; see [`_DMTEnvCache`](@ref)).
"""
function _invalidate_env!(cache::_DMTEnvCache, lo::Integer, hi::Integer)
  cache.lvalid = clamp(min(cache.lvalid, Int(lo) - 2), 0, cache.n)
  cache.rvalid = clamp(max(cache.rvalid, Int(hi) + 2), 1, cache.n + 1)
  return cache
end

"""
    _orthogonalize_env!(cache, psi, bond)

Move the orthogonality center to `bond`, invalidating the re-gauged range. `orthogonalize!`
re-gauges the tensors between the old center and `bond`; the old center lies in
`[leftlim, rightlim]`, so `[min(leftlim,bond), max(rightlim,bond)]` bounds the touched sites.
"""
function _orthogonalize_env!(cache::_DMTEnvCache, psi::MPS, bond::Integer)
  ll = ITensorMPS.leftlim(psi)
  rl = ITensorMPS.rightlim(psi)
  orthogonalize!(psi, bond)
  lo = clamp(min(ll, Int(bond)), 1, cache.n)
  hi = clamp(max(rl, Int(bond)), 1, cache.n)
  _invalidate_env!(cache, lo, hi)
  return psi
end

function _assert_env_matches(label, cached::ITensor, fresh::ITensor)
  ITensors.hassameinds(cached, fresh) ||
    error("DMT env verify ($label): index mismatch\n  cached=$(inds(cached))\n  fresh =$(inds(fresh))")
  difference = norm(cached - fresh)
  scale = max(norm(fresh), one(real(float(one(eltype(fresh))))))
  difference <= sqrt(eps(Float64)) * scale ||
    error("DMT env verify ($label): norm(diff)=$difference exceeds tolerance (scale=$scale)")
  return nothing
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
  orthogonalize::Bool=true,
  cache::Union{Nothing,_DMTEnvCache}=nothing,
)
  maxdim > 0 || return psi
  connector_buffer >= 0 || throw(ArgumentError("DMT connector_buffer must be nonnegative"))
  connector_buffer <= maxdim || throw(ArgumentError("DMT connector_buffer must be <= maxdim"))
  1 <= bond < length(psi) || throw(ArgumentError("DMT bond must lie in 1:length(psi)-1"))
  current_link = linkind(psi, bond)
  isnothing(current_link) && return psi
  dim(current_link) <= maxdim && return psi

  if orthogonalize
    isnothing(cache) ? orthogonalize!(psi, bond) : _orthogonalize_env!(cache, psi, bond)
  end
  left_site = siteind(psi, bond)
  right_site = siteind(psi, bond + 1)

  # The identity/trace environments depend only on the current `psi[1:bond-1]` / `psi[bond+2:N]`
  # tensors. Rebuilding both from scratch every call is the dominant per-sweep cost (the O(N^2)
  # note on `_dmt_window_truncate!`); a supplied `cache` memoizes them and rebuilds only the part
  # invalidated since the last bond, amortizing the cost to O(1) for a local sweep. Because the
  # rest of this function is byte-for-byte unchanged, a cached env equal to the rebuild yields
  # identical truncation output -- asserted under `_DMT_VERIFY_ENVS[]`.
  if isnothing(cache)
    left_env = _left_identity_environment(psi, bond - 1)
    right_env = _right_identity_environment(psi, bond + 2)
  else
    left_env = _left_env_at!(cache, psi, bond - 1)
    right_env = _right_env_at!(cache, psi, bond + 2)
    if _DMT_VERIFY_ENVS[]
      _assert_env_matches("left b=$bond", left_env, _left_identity_environment(psi, bond - 1))
      _assert_env_matches("right b=$bond", right_env, _right_identity_environment(psi, bond + 2))
    end
  end

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
  # The repair SVD rewrote sites `bond` and `bond+1`; invalidate them so the cache stays correct
  # independent of what runs next (each operation owns its own footprint).
  isnothing(cache) || _invalidate_env!(cache, bond, bond + 1)
  return psi
end

function _dmt_window_truncate!(psi::MPS, start::Integer, span::Integer; maxdim::Integer, cutoff::Real, direction::Symbol, connector_buffer::Integer, cache::Union{Nothing,_DMTEnvCache}=nothing)
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
  #
  # PERFORMANCE (perf-1, implemented). Rebuilding both identity/trace environments from scratch
  # in every `_dmt_bond_truncate!` (`_left_identity_environment` over `1:bond-1`,
  # `_right_identity_environment` over `bond+2:N`) costs O(N) per bond, so a full ~N-bond
  # schedule sweep is O(N^2 * chi^2) -- the dominant term for long chains at moderate chi (the
  # regime reached pushing PXP runs to the times needed for the asymptotic z=3/2). When a
  # `cache::_DMTEnvCache` is threaded through the sweep (`dmt_evolve!` builds one and passes it
  # down), each bond instead extends a memoized running environment by only the sites mutated
  # since the previous bond, amortizing to O(1) per bond -> O(N * chi^2) per sweep.
  #
  # The cache does NOT change what this function computes: `orthogonalize=true` and the entire
  # `_dmt_bond_truncate!` operation sequence are unchanged, and a cached environment is by
  # construction the same contraction of the same `psi` tensors the rebuild would produce -- so
  # the truncation output is bit-for-bit identical. This sidesteps the three problems that made
  # the optimization look hard:
  #   1. Gauge consistency: the env value depends only on the current tensors, not on which gauge
  #      produced them, and the per-bond `orthogonalize!` is retained -- so the cache is correct
  #      in any gauge, not only a frozen monotonic prefix/suffix.
  #   2. The gate primitive: `tebd_evolve!` (ITensorMPS `product`) re-gauges opaquely, but it
  #      only ever touches the path between the old orthocenter and the gate window; the caller
  #      invalidates that bounded range (via `leftlim`/`rightlim`) rather than tracking `product`
  #      internally.
  #   3. Schedule shape: invalidation is per-operation and footprint-based, so the non-monotonic,
  #      mixed-span, overlapping PXP schedule needs no special handling -- it stays O(1)/bond
  #      because every gate's footprint is its O(span) window (verified empirically).
  # Correctness is guarded by `_DMT_VERIFY_ENVS[] = true`, which asserts every cached env equals
  # the from-scratch rebuild (index identity + norm of difference) across the ED-oracle suite.
  for bond in _dmt_truncation_bonds(start, span, direction)
    _dmt_bond_truncate!(
      psi,
      bond;
      maxdim=maxdim,
      cutoff=cutoff,
      direction=direction,
      connector_buffer=connector_buffer,
      orthogonalize=true,
      cache=cache,
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
  cache::Union{Nothing,_DMTEnvCache}=nothing,
)
  start = _bond_start(bond)
  span = _operator_span_at(psi, gate, start)
  _validate_dmt_step(psi, gate, start, span, direction, Int(maxdim), Int(connector_buffer))
  # `tebd_evolve!` (ITensorMPS `product`) re-gauges the path between the old orthocenter and the
  # gate window. Capture the limits first so the env cache can invalidate that bounded range.
  if !isnothing(cache)
    ll = ITensorMPS.leftlim(psi)
    rl = ITensorMPS.rightlim(psi)
  end
  tebd_evolve!(psi, gate, start; maxdim=Int(gate_maxdim), cutoff=0.0)
  if !isnothing(cache)
    _invalidate_env!(
      cache,
      clamp(min(ll, start), 1, length(psi)),
      clamp(max(rl, start + span - 1), 1, length(psi)),
    )
  end
  _dmt_window_truncate!(
    psi,
    start,
    span;
    maxdim=Int(maxdim),
    cutoff=cutoff,
    direction=direction,
    connector_buffer=Int(connector_buffer),
    cache=cache,
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
  # One environment cache threaded through the whole call: each `dmt_step!` mutates `psi`
  # locally and invalidates only the touched range, so the cache stays consistent across the
  # forward sweep, the forward->reverse turnaround, and successive `nstep` passes (no mutation
  # happens between steps that the step itself does not record). See `_DMTEnvCache` / the
  # perf-1 note on `_dmt_window_truncate!`.
  cache = _DMTEnvCache(psi)
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
        cache=cache,
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
        cache=cache,
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
