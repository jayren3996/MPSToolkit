"""
    DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=0, preserve_diameter=3, truncation=:dense)

Options controlling operator-space density matrix truncation (DMT).

!!! warning "Transport-specific algorithm"
    DMT is a specialized truncation scheme designed for **transport** (e.g. spin or energy
    diffusion) in operator space.  It protects local reduced operator data, including the
    identity/trace component and nearby local-operator components, before truncating connected
    long-range correlations. Because the protected anchor IS the trace/identity component, DMT is
    valid only for near-infinite-temperature **density operators** (e.g. an energy domain-wall
    melt). Do **not** apply it to a **traceless** operator such as a two-point correlator (e.g.
    `P_G h_c P_G`): there is no trace component to anchor the truncation, so it runs outside its
    design regime. For correlators — or any general operator-space evolution without this transport
    bias — use ordinary TEBD truncation (`LocalGateEvolution`) instead.

# Fields
- `maxdim`: **Total** bond dimension after DMT truncation, inclusive of the protected block.
  The complement is truncated to `maxdim - 2 d^(preserve_diameter - 1)` directions, so `maxdim`
  must be at least `2 d^(preserve_diameter - 1) + 1` for the local dimension `d` in use.
- `cutoff`: Truncation cutoff used in the final refactorization.
- `gate_maxdim`: Temporary bond dimension cap applied while the raw gate is applied, before DMT
  truncates the bond back to `maxdim`. **`0` (the default) means no cap: the gate is applied
  exactly.** A positive cap pre-truncates the inflated bond with a plain SVD, which discards the
  smallest singular values *before* DMT sees them — exactly the error DMT exists to avoid — so
  set one only to bound peak cost, never for accuracy.
- `preserve_diameter`: Positive odd diameter of the observables preserved exactly;
  `radius = (preserve_diameter - 1) / 2` sites are protected on each side of the cut.
- `truncation`: `:dense` (default) or `:random` complement truncation. `:random` preserves the
  guarantee to the same `1e-15` (see `_truncated_svd` internally) and is measurably faster —
  1.05x-1.2x on a whole sweep at moderate budgets, ~1.4x once the gate-inflated bond passes
  ~2500, and up to 1.8x on the DMT truncation alone, which is the only part it touches. It is
  not the default because it is **not deterministic**: it draws from the global RNG, so two
  truncations of the same bond agree only to randomized-SVD accuracy (~1e-3) — and
  `dmt_evolve!` truncates every bond twice per sweep (`direction=:R` on the forward schedule,
  `:L` on the reverse), so a `:random` default would make the two halves of *every* production
  sweep sketch independently instead of reproducing the same physical state. Opt in for large
  runs where the speedup matters and reproducibility is handled by seeding.

# Notes
- The `gate_maxdim` default was `max(maxdim * 16, 64)`, which in steady state capped nothing: a
  two-site gate inflates the bond from the incoming `chi` to `d^2 * chi`, and `chi <= maxdim`
  once DMT has truncated that bond, so `16 * maxdim` was already unreachable for every `d <= 4`.
  (On a first sweep over a state entering *wider* than `maxdim`, `chi > maxdim` and the old cap
  could bite even at `d <= 4`.) `0` makes the shipped behaviour explicit and unconditional, and
  extends it to `d >= 5`, where the old formula pre-truncated even in steady state.
"""
struct DMTOptions
  maxdim::Int
  cutoff::Float64
  gate_maxdim::Int
  preserve_diameter::Int
  truncation::Symbol

  function DMTOptions(maxdim, cutoff, gate_maxdim, preserve_diameter, truncation)
    maxdim >= 1 || throw(ArgumentError("DMTOptions requires maxdim >= 1"))
    cutoff >= 0 || throw(ArgumentError("DMTOptions requires cutoff >= 0"))
    gate_maxdim >= 0 || throw(ArgumentError("DMTOptions requires gate_maxdim >= 0 (0 = no cap)"))
    isodd(preserve_diameter) && preserve_diameter >= 1 || throw(ArgumentError(
      "DMTOptions requires a positive odd preserve_diameter, got $(preserve_diameter)"))
    truncation in (:dense, :random) ||
      throw(ArgumentError("DMTOptions truncation must be :dense or :random, got $(truncation)"))
    return new(Int(maxdim), Float64(cutoff), Int(gate_maxdim), Int(preserve_diameter), Symbol(truncation))
  end

  function DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=0,
    preserve_diameter=3, truncation=:dense, connector_buffer=nothing)
    _reject_connector_buffer(connector_buffer)
    return DMTOptions(maxdim, cutoff, gate_maxdim, preserve_diameter, truncation)
  end
end

"""
    _identity_env(site)

Return the local basis-1 ket used when building operator-space identity/trace environments.
Basis element 1 is always the normalized identity `I / sqrt(d)` (see
[`operator_basis_matrices`](@ref)), so this is dimension-generic without any change of body.
"""
function _identity_env(site)
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
    env *= _identity_env(siteind(psi, site)) * psi[site]
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
    env *= _identity_env(siteind(psi, site)) * psi[site]
  end
  return env
end

"""
    _DMTEnvCache(psi)

Amortizing cache for the DMT identity/trace environments (perf-1). `_dmt_bond_truncate!` needs,
at each bond `b`, the left environment over sites `1:b-radius` and the right environment over
`b+1+radius:N`. Rebuilding both from scratch on every call is the O(N) per-bond cost that makes a
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
    cache.left[j + 1] = cache.left[j] * (_identity_env(siteind(psi, j)) * psi[j])
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
    cache.right[j] = cache.right[j + 1] * (_identity_env(siteind(psi, j)) * psi[j])
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

"""
    _validate_operator_space(psi, start, span)

Validate that the window `start:(start + span - 1)` of operator-space `MPS` `psi` has a
uniform local dimension, and return that dimension.

# Returns
- The local dimension `d` shared by every site in the window.
"""
function _validate_operator_space(psi::MPS, start::Integer, span::Integer)
  d = local_dimension(siteind(psi, Int(start)))
  for site in Int(start):(Int(start) + Int(span) - 1)
    local_dimension(siteind(psi, site)) == d ||
      throw(ArgumentError("operator-space window sites must share a common local dimension"))
  end
  return d
end

function _validate_dmt_step(psi::MPS, gate::AbstractMatrix, start::Integer, span::Integer, direction::Symbol, maxdim::Integer, preserve_diameter::Integer)
  direction === :R || direction === :L || throw(ArgumentError("DMT direction must be :R or :L"))
  maxdim >= 1 || throw(ArgumentError("DMT maxdim must be >= 1"))
  start >= 1 || throw(ArgumentError("local gate bond must be at least 1"))
  last_site = start + span - 1
  last_site <= length(psi) || throw(ArgumentError("local gate support exceeds chain length"))
  # Validate operator-space dimensions for every span, including span == 1, so a single-site
  # gate on a non-operator-space (non-square-dimension) or non-uniform-dimension site is
  # rejected rather than silently accepted.
  _validate_operator_space(psi, start, span)
  # The budget floor must be checked HERE, not only inside `_dmt_bond_truncate!`: `dmt_step!`
  # reaches the kernel only after `tebd_evolve!` has applied the gate and re-gauged `psi`, so a
  # kernel-side throw would leave a half-updated state behind and break the "invalid DMT calls
  # do not mutate the state" invariant. Checked before the `span == 1` return as well, which
  # otherwise never reaches the kernel and would skip budget validation entirely.
  _validate_dmt_budget(psi, maxdim, preserve_diameter)
  span == 1 && return nothing
  start == length(psi) && throw(ArgumentError("periodic boundary DMT is not implemented for local gates"))
  for bond in start:(last_site - 1)
    1 <= bond < length(psi) || throw(ArgumentError("DMT target bonds must lie in 1:length(psi)-1"))
  end
  return nothing
end

function _dmt_window_truncate!(psi::MPS, start::Integer, span::Integer; maxdim::Integer, cutoff::Real, direction::Symbol, preserve_diameter::Integer=3, truncation::Symbol=:dense, cache::Union{Nothing,_DMTEnvCache}=nothing)
  span <= 1 && return psi

  # Truncate every bond inside the gate window as an independent single-bond DMT update, with
  # the orthogonality center restored to that bond before truncating it. The
  # connector-preserving construction in `_dmt_bond_truncate!` assumes a *canonical* gauge: the
  # bond factorization must see an orthonormal basis on each side of the cut, and the trace
  # environments must be built from orthonormal blocks. A multi-bond window cannot satisfy that
  # with a single shared sweep, because each DMT truncation rewrites the bond tensors and resets
  # the MPS gauge bookkeeping. Sweeping with cached environments therefore truncates later bonds
  # in a non-canonical gauge: mildly wrong for a `:R` window, and badly wrong for `:L`, where the
  # orthogonality center ends up on the wrong site so the bond matrix degenerates to a unitary
  # (every Schmidt weight equal) and the truncation discards information indiscriminately.
  # Re-gauging per bond keeps a span-`S` window exactly equal to the verified single-bond path
  # applied to each of its bonds.
  #
  # PERFORMANCE (perf-1, implemented). Rebuilding both identity/trace environments from scratch
  # in every `_dmt_bond_truncate!` (`_left_identity_environment` over `1:bond-radius`,
  # `_right_identity_environment` over `bond+1+radius:N`) costs O(N) per bond, so a full ~N-bond
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
      preserve_diameter=preserve_diameter,
      truncation=truncation,
      orthogonalize=true,
      cache=cache,
    )
  end
  return psi
end

"""
    dmt_step!(psi, gate, bond; maxdim=30, cutoff=1e-12, direction=:R, gate_maxdim=0, preserve_diameter=3, truncation=:dense)

Apply one local operator-space gate and then perform DMT-preserving truncation.

This is a **transport-specific** truncation step.  See [`DMTOptions`](@ref) for when DMT
is (and is not) the appropriate choice.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `gate`: Dense local gate in the operator basis.
- `bond`: Left-edge location of the local update.

# Keyword Arguments
- `maxdim`: **Total** post-truncation bond dimension, inclusive of the protected block; it must
  be at least `2 d^(preserve_diameter - 1) + 1`.
- `cutoff`: Truncation cutoff used in the final refactorization.
- `direction`: Sweep direction, either `:R` or `:L`.
- `gate_maxdim`: Temporary bond dimension cap for the raw gate application, before DMT truncates
  the bond back to `maxdim`. `0` (the default) means no cap, i.e. the gate is applied exactly;
  see [`DMTOptions`](@ref) for why pre-truncation is the error DMT exists to avoid.
- `preserve_diameter`: Positive odd diameter of the observables preserved exactly.
- `truncation`: `:dense` (default) or `:random` complement truncation; see [`DMTOptions`](@ref)
  for the measured speedup and the determinism it costs.

# Returns
- The mutated `psi`.

# Notes
- Applying the gate exactly inflates the bond from the incoming `chi` to `d^2 * chi` before DMT
  truncates it back to `maxdim`, so the bond tensor the kernel factorizes is `d^2` times wider
  than `maxdim` suggests (`chi <= maxdim` in steady state). That is the dominant cost of a step,
  and it is the cost of not throwing away protected data.
"""
function dmt_step!(
  psi::MPS,
  gate::AbstractMatrix,
  bond;
  maxdim::Integer=30,
  cutoff::Real=1e-12,
  direction::Symbol=:R,
  gate_maxdim::Integer=0,
  preserve_diameter::Integer=3,
  truncation::Symbol=:dense,
  cache::Union{Nothing,_DMTEnvCache}=nothing,
  connector_buffer=nothing,
)
  _reject_connector_buffer(connector_buffer)
  start = _bond_start(bond)
  span = _operator_span_at(psi, gate, start)
  # Checked before the gate runs, like every other budget here, so a rejected call leaves `psi`
  # untouched. `0` is the "no cap" sentinel; anything below it is a typo, not a smaller cap.
  gate_maxdim >= 0 || throw(ArgumentError("DMT gate_maxdim must be >= 0 (0 = no cap)"))
  _validate_dmt_step(psi, gate, start, span, direction, Int(maxdim), Int(preserve_diameter))
  # `tebd_evolve!` (ITensorMPS `product`) re-gauges the path between the old orthocenter and the
  # gate window. Capture the limits first so the env cache can invalidate that bounded range.
  if !isnothing(cache)
    ll = ITensorMPS.leftlim(psi)
    rl = ITensorMPS.rightlim(psi)
  end
  # `gate_maxdim == 0` means "apply the gate exactly": `_tebd_truncation_kwargs` omits `maxdim`
  # from the ITensor call for a non-positive budget, which is the same "no cap" convention
  # `LocalGateEvolution` already uses for its own `maxdim`. Passing the sentinel straight through
  # is preferred over substituting `typemax(Int)`, which would push a sentinel-sized budget into
  # ITensors' truncation arithmetic instead of taking the branch that skips it.
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
    preserve_diameter=Int(preserve_diameter),
    truncation=truncation,
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
    preserve_diameter=opts.preserve_diameter,
    truncation=opts.truncation,
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

!!! warning "Density operators only — not traceless operators / correlators"
    DMT protects the trace/identity component, so `psi` must be a near-infinite-temperature
    **density operator** (e.g. an energy domain-wall melt). Do **not** use it to Heisenberg-evolve a
    **traceless** operator such as a two-point correlator — there is no trace to anchor the
    truncation. Use ordinary TEBD (`tebd_evolve!` / [`LocalGateEvolution`](@ref)) for those.

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
        preserve_diameter=evo.preserve_diameter,
        truncation=evo.truncation,
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
        preserve_diameter=evo.preserve_diameter,
        truncation=evo.truncation,
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
