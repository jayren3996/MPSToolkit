"""
    operator_trace(rho)

Return the physical trace of a vectorized operator stored as an operator-space `MPS`, via
`tr(ρ) = (√d)^N c_{I…I}`, where `c_{I…I}` is the amplitude of the all-identity string and `d`
is the local Hilbert space dimension recovered from `rho`'s site indices — one
identity-environment contraction.

# Arguments
- `rho`: Operator-space `MPS` on sites of a uniform local dimension `d`, typically from
  [`operator_siteinds`](@ref).

# Returns
- The trace as a `ComplexF64` (real up to numerical noise for Hermitian operators).

# Notes
- The `(√d)^N` prefactor overflows `Float64` for `N ≳ 2048 / log2(d)`, so absolute traces of
  very large operators are not representable (far beyond any feasible operator-space MPS
  size); normalized ratio observables are unaffected because the factor cancels.
"""
function operator_trace(rho::MPS)
  d = _validate_operator_space(rho, 1, length(rho))
  return Float64(d)^(length(rho) / 2) * scalar(_right_identity_environment(rho, 1))
end

"""
    pauli_trace(rho)

Spin-1/2 case of [`operator_trace`](@ref): return the physical trace of a vectorized operator
stored as a Pauli-basis `MPS`.

# Arguments
- `rho`: Operator-space `MPS` on dimension-4 Pauli sites.

# Returns
- The trace as a `ComplexF64` (real up to numerical noise for Hermitian operators).
"""
pauli_trace(rho::MPS) = operator_trace(rho)

"""
    _operator_window_cap(window_sites, op, d)

Build the window cap tensor whose entry at the multi-site basis label `α` is `tr(P_α O)`, so
that contracting it against an operator-space MPS window (with identity caps elsewhere) yields
`tr(ρ O)` up to the `(√d)` factors handled by the callers.

# Notes
- `tr(P_{α_1} ⊗ … ⊗ P_{α_s} · O)` is *separable* in the site label `α_k`, so the cap is built by
  contracting `O` — reshaped into its `2 s` physical legs — against one `d^2 x d x d` weight
  tensor per site, at `O(s d^(2s + 2))` cost. Enumerating the `d^(2s)` labels and forming
  `tr(kron(...) * O)` for each is `O(d^(5s))` and is what makes higher `d` unaffordable: at
  `d = 4, s = 3` the label loop measures 0.254 s per term against 0.000254 s here, which takes
  the diameter-3 preservation sweep at `d = 4` from 122 s to 18 s per pass. Both forms agree to
  1-2 ulp: the largest relative difference measured is 3.39e-16, over `d` in `2:5` at spans 1-3
  (complex, Hermitian, real and sparse inputs), `d = 3` at span 4, and `d = 2` at spans 4-5.
- `kron` orders the first factor slowest, so the reshaped row legs run `(s, s-1, …, 1)`; hence
  the `reverse` when the legs are attached.
"""
function _operator_window_cap(window_sites, op::AbstractMatrix, d::Integer)
  span = length(window_sites)
  local_dim = Int(d)
  size(op) == (local_dim^span, local_dim^span) ||
    throw(ArgumentError("window cap operator size must match the window span"))
  local_basis = operator_basis_matrices(local_dim)
  row = [Index(local_dim, "OperatorCapRow,n=$(k)") for k in 1:span]
  col = [Index(local_dim, "OperatorCapCol,n=$(k)") for k in 1:span]
  legs = reshape(Matrix{ComplexF64}(op), ntuple(_ -> local_dim, 2 * span))
  cap = ITensor(legs, reverse(row)..., reverse(col)...)
  for k in 1:span
    # weights[α, j, i] = P_α[i, j], the transpose that turns the contraction into a trace.
    weights = ITensor(ComplexF64, window_sites[k], row[k], col[k])
    for label in 1:(local_dim^2)
      element = local_basis[label]
      for i in 1:local_dim, j in 1:local_dim
        iszero(element[i, j]) && continue
        weights[window_sites[k] => label, row[k] => j, col[k] => i] = element[i, j]
      end
    end
    cap *= weights
  end
  return cap
end

"""
    _validated_operator_windows(rho, terms, d)

Validate a collection of `(start, op)` measurement terms against the chain length of `rho`,
inferring each term's site span from `op`'s size at local dimension `d`, and return the
`(start, span)` windows in input order.
"""
function _validated_operator_windows(rho::MPS, terms, d::Integer)
  windows = Vector{Tuple{Int,Int}}(undef, length(terms))
  for (k, (start, op)) in enumerate(terms)
    size(op, 1) == size(op, 2) || throw(ArgumentError("dense local operator must be square"))
    1 <= start <= length(rho) || throw(ArgumentError("term start must lie in 1:$(length(rho))"))
    span = _operator_span(size(op, 1), d)
    start + span - 1 <= length(rho) || throw(ArgumentError("term support exceeds chain length"))
    windows[k] = (Int(start), span)
  end
  return windows
end

"""
    _log_trace_resolution(rho, d, identity_coefficient)

Return `(log|tr(ρ)|, log ||ρ||_HS)` for an operator-space `MPS`, both in natural log and both
computed without ever forming the quantity itself.

# Arguments
- `rho`: Operator-space `MPS` on sites of a uniform local dimension `d`.
- `d`: That local dimension.
- `identity_coefficient`: The all-identity amplitude `c_{I…I}` of `rho`, i.e. the scalar of its
  full identity-environment contraction, so that `tr(ρ) = (√d)^N c_{I…I}`.

# Returns
- A `Tuple{Float64,Float64}` of the two logs. `log|tr(ρ)|` is `-Inf` for an exactly traceless
  operator; `log ||ρ||_HS` is `-Inf` only for the zero operator.

# Notes
- Log space is not cosmetic here. Both quantities overflow `Float64` for states that are
  otherwise perfectly ordinary: `(√d)^N` alone passes `1e308` at `N ≈ 2048 / log2(d)`, and for
  the cold thermal product state `ρ = ⊗_j e^{-Z_j}` at `N = 400` the Hilbert-Schmidt norm is
  `1e175`, so `norm(rho)` — which forms `||ρ||²` — overflows to `NaN`.
- The two logs are returned rather than their difference so that a caller can distinguish
  `NaN`/`Inf` inputs (an unrepresentable operator) from a legitimately tiny ratio.
- Only the norm side is protected. `lognorm` rescales progressively, but `identity_coefficient`
  arrives as a plain linear-space contraction, so it can overflow before this function is
  reached: for the *unnormalized* product operator `ρ = ⊗_j e^{-β Z_j}`, whose identity
  coefficient is exactly `(√d cosh β)^N`, that happens at `β ≥ 2.11` for `N = 400` and `β ≥ 3.90`
  for `N = 200` (`d = 2`). `_reject_unresolvable_trace` reports the resulting `Inf`/`NaN`
  as its own error. A `normalize!`d operator is immune at any feasible size: positivity confines
  its identity coefficient to `[d^(-N/2), 1]`, which stays inside `Float64` for the same
  `N ≲ 2048 / log2(d)` that bounds the `(√d)^N` prefactor everywhere else in this file. This is
  why [`operator_gibbs_state`](@ref), which normalizes, never meets it, and why an operator built
  directly with [`operator_product_state`](@ref) should be normalized before it is measured.
"""
function _log_trace_resolution(rho::MPS, d::Integer, identity_coefficient::Number)
  log_trace = (length(rho) / 2) * log(Float64(d)) + log(abs(identity_coefficient))
  return (log_trace, lognorm(rho))
end

"""
    _reject_unresolvable_trace(rho, d, identity_coefficient)

Throw an `ArgumentError` unless `tr(ρ)` is resolvable against `ρ`'s own Hilbert-Schmidt norm,
i.e. unless `|tr(ρ)| > sqrt(eps) ||ρ||_HS`. Used only on the `normalize=true` path of
[`operator_expectation_profile`](@ref), where every result is divided by `tr(ρ)`.

# Arguments
- `rho`, `d`, `identity_coefficient`: See [`_log_trace_resolution`](@ref).

# Notes
- The comparison is in *physical* units: `tr(ρ) = (√d)^N c_{I…I}` carries a `(√d)^N` that the
  identity coefficient does not. Dropping that factor — comparing `|c_{I…I}|` against
  `sqrt(eps) ||ρ||_HS`, as this check did until 2026-08-10 — puts a `d^(N/2)` on the wrong side
  and rejects legitimate cold thermal states as the chain grows: for the positive,
  bond-dimension-1 product state `ρ = ⊗_j e^{-β Z_j}`, whose trace is exactly `(2 cosh β)^N`,
  the old form threw at `(β, N) = (0.35, 400)`, `(0.5, 240)` and `(1.0, 120)`.
- In physical units the check cannot reject a positive operator *on conditioning grounds*: for
  `λ_i ≥ 0`, `tr(ρ) = Σ_i λ_i ≥ sqrt(Σ_i λ_i²) = ||ρ||_HS`, so every positive `ρ` clears the
  threshold by the full `1 / sqrt(eps) ≈ 6.7e7`, independent of `N`, `d` and temperature. The
  bound is tight only at rank 1 (`β → ∞`), where the ratio is measured at exactly `1.0` for
  `N = 400`. That guarantee is what a magnitude threshold cannot have and is why the check is
  stated this way rather than as a cancellation test on the environment sweep: the identity
  *component* of a cold Gibbs state decays exponentially in `N`, but its *trace* does not.
- The guarantee is about the *comparison*, not about `Float64` range. An **unnormalized** operator
  can still be rejected because its identity coefficient overflows before the comparison happens
  — `β ≥ 2.11` at `N = 400` for `ρ = ⊗_j e^{-β Z_j}` — which is reported as the representability
  error below, not as a traceless operator. See [`_log_trace_resolution`](@ref); `normalize!` the
  operator (or build it with [`operator_gibbs_state`](@ref), which does) and it cannot arise.
- What it still rejects is the case the check exists for: a traceless operator (a transport
  current, an evolved two-point correlator) whose post-truncation trace is an `O(eps)`
  cancellation residue rather than exactly zero, where normalizing amplifies every entry by
  `~1/eps`. A DMT-evolved [`pauli_domain_wall_state`](@ref) measures `|tr(ρ)| / ||ρ||_HS` at
  7e-14 (`N = 12`) to 2e-12 (`N = 20`) after one sweep, still 4e-11 after twelve.
- The residue is not bounded uniformly, and the check is a statement about *half precision*, not
  about tracelessness: the same `N = 40` melt at `maxdim = 24` reaches `3e-8` by sweep twelve and
  is accepted. That is the criterion behaving as written — by then the trace has more than half
  its digits, because the truncation error itself has grown past `sqrt(eps) ||ρ||_HS` — and it
  cannot be tuned away without rejecting positive operators again, whose guaranteed floor is
  `1.0`. Traceless transport operators should carry `normalize=false` by construction rather
  than rely on being caught.
- A non-finite trace or norm is an error, not a pass. The `NaN` that `norm(rho)` used to return
  for the `N = 400` state above made `x <= sqrt(eps) * NaN` evaluate `false`, so the old check
  went silently vacuous exactly where it was firing hardest at smaller `N`.
"""
function _reject_unresolvable_trace(rho::MPS, d::Integer, identity_coefficient::Number)
  log_trace, log_norm = _log_trace_resolution(rho, d, identity_coefficient)
  (isnan(log_trace) || isnan(log_norm) || log_trace == Inf || log_norm == Inf) &&
    throw(ArgumentError("normalized operator-space expectations require a representable trace and norm; got log|tr(rho)| = $(log_trace) and log||rho||_HS = $(log_norm), so the operator carries NaN or overflowed entries"))
  log_trace <= log_norm + log(sqrt(eps(Float64))) &&
    throw(ArgumentError("normalized operator-space expectations require a trace resolvable against the operator's own norm; log10|tr(rho)| = $(log_trace / log(10)) against log10||rho||_HS = $(log_norm / log(10)) leaves |tr(rho)| at or below sqrt(eps) * ||rho||_HS, i.e. a cancellation residue rather than a trace (use normalize=false for a traceless operator)"))
  return nothing
end

"""
    operator_expectation_profile(rho, terms; normalize=true)

Evaluate `tr(ρ O_k)` (optionally over `tr(ρ)`) for a list of dense local operators against a
vectorized operator `rho`, in one O(N) sweep with cumulative identity environments.

# Arguments
- `rho`: Operator-space `MPS` on sites of a uniform local dimension `d`.
- `terms`: Collection of `(start, op)` pairs; each `op` is a dense physical operator on
  `d^span` dimensions acting on the sites `start:(start + span - 1)`. The pairs may be given
  in any order; results follow the input order.

# Keyword Arguments
- `normalize`: If `true` (default), return `tr(ρ O_k) / tr(ρ)`; otherwise return the
  unnormalized `tr(ρ O_k)`. The unnormalized branch carries a `(√d)^N` factor that overflows
  `Float64` for `N ≳ 2048 / log2(d)` (far beyond feasible MPS sizes); the normalized ratio is
  immune. `true` throws an `ArgumentError` when `|tr(ρ)| ≤ sqrt(eps) ||ρ||_HS`, which is the
  traceless case (see `_reject_unresolvable_trace`); measure a traceless operator with
  `normalize=false`. A positive `ρ` — any thermal state, at any temperature and chain length —
  satisfies `tr(ρ) ≥ ||ρ||_HS` and is never rejected for being ill-conditioned. It can still be
  rejected for being unrepresentable: an *unnormalized* `⊗_j e^{-β Z_j}` overflows its identity
  coefficient at `β ≥ 2.11` for `N = 400`, so normalize an operator built directly with
  [`operator_product_state`](@ref) before measuring it ([`operator_gibbs_state`](@ref) already
  does).

# Returns
- A `Vector{ComplexF64}` of expectation values. For Hermitian `ρ` and `O_k` the entries are
  real up to numerical noise.

# Notes
- This is the energy-density-profile measurement for operator-space transport runs: pass the
  Hamiltonian terms (e.g. from [`pxp_term_hamiltonian`](@ref)) and take the real part.
- Every site of `rho` must share the same local dimension `d`; it is recovered once from the
  chain (not re-derived per term) and reused for every window in `terms`.
"""
function operator_expectation_profile(rho::MPS, terms; normalize::Bool=true)
  nsites = length(rho)
  d = _validate_operator_space(rho, 1, nsites)
  # Materialize so a lazy (non-indexable) generator of (start, op) pairs is accepted; the windowed
  # sweep below indexes terms[k] after sorting by window start.
  terms = collect(terms)
  isempty(terms) && return ComplexF64[]
  windows = _validated_operator_windows(rho, terms, d)

  right = Vector{ITensor}(undef, nsites + 1)
  right[nsites + 1] = ITensor(1.0)
  for site in nsites:-1:1
    right[site] = right[site + 1] * (_identity_env(siteind(rho, site)) * rho[site])
  end
  denominator = scalar(right[1])
  # Reject a numerically-negligible (not just exactly-zero) trace relative to the operator
  # scale, which is the O(eps) cancellation residue a traceless operator has after truncation.
  # Evaluated only when `normalize` is true: the `normalize=false` branch never divides by the
  # trace and must not pay for the norm.
  normalize && _reject_unresolvable_trace(rho, d, denominator)

  results = Vector{ComplexF64}(undef, length(terms))
  left = ITensor(1.0)
  absorbed = 0
  # `sortperm(windows)` orders lexicographically, so `start` is still non-decreasing (which the
  # cumulative `left` prefix requires) and terms sharing a window are now adjacent. The window
  # contraction `left * rho[start:start+span-1] * right` does not depend on the operator, so it
  # is built once per *window* rather than once per term: it is the reduced operator of that
  # window, and every term on it costs only the `d^(2 span)` cap contraction afterwards. A
  # diameter sweep measures many operators per window (5115 width-5 probes over 5 windows in the
  # DMT preservation tests), where this is the difference between 205 s and 12 s per pass.
  window_key = (0, 0)
  window_sites = Vector{typeof(siteind(rho, 1))}()
  reduced = ITensor(1.0)
  for k in sortperm(windows)
    start, span = windows[k]
    while absorbed < start - 1
      absorbed += 1
      left = left * (_identity_env(siteind(rho, absorbed)) * rho[absorbed])
    end
    if (start, span) != window_key
      window_block = rho[start]
      for site in (start + 1):(start + span - 1)
        window_block *= rho[site]
      end
      window_sites = [siteind(rho, s) for s in start:(start + span - 1)]
      reduced = left * window_block * right[start + span]
      window_key = (start, span)
    end
    cap = _operator_window_cap(window_sites, last(terms[k]), d)
    raw = scalar(cap * reduced)
    results[k] = normalize ? raw / (Float64(d)^(span / 2) * denominator) : Float64(d)^((nsites - span) / 2) * raw
  end
  return results
end

"""
    pauli_expectation_profile(rho, terms; normalize=true)

Spin-1/2 case of [`operator_expectation_profile`](@ref): evaluate `tr(ρ O_k)` (optionally over
`tr(ρ)`) for a list of dense local operators against a vectorized operator `rho` on Pauli sites.

# Arguments
- `rho`: Operator-space `MPS` on dimension-4 Pauli sites.
- `terms`: Collection of `(start, op)` pairs; each `op` is a dense physical operator on `2^span`
  dimensions.

# Keyword Arguments
- `normalize`: See [`operator_expectation_profile`](@ref).

# Returns
- A `Vector{ComplexF64}` of expectation values.
"""
pauli_expectation_profile(rho::MPS, terms; normalize::Bool=true) =
  operator_expectation_profile(rho, terms; normalize=normalize)

"""
    operator_expectation(rho, op, start; normalize=true)

Evaluate `tr(ρ O) / tr(ρ)` (or the unnormalized `tr(ρ O)` with `normalize=false`) for one
dense local operator `op` applied at `start` against a vectorized operator `rho`.

See [`operator_expectation_profile`](@ref) for the batched O(N) version and the conventions.

# Arguments
- `rho`: Operator-space `MPS` on sites of a uniform local dimension `d`.
- `op`: Dense local operator of size `d^span x d^span`.
- `start`: Left-edge site index of `op`'s support.

# Keyword Arguments
- `normalize`: See [`operator_expectation_profile`](@ref).

# Returns
- A `ComplexF64` expectation value (real up to numerical noise for Hermitian inputs).
"""
function operator_expectation(rho::MPS, op::AbstractMatrix, start::Integer; normalize::Bool=true)
  return only(operator_expectation_profile(rho, [(Int(start), op)]; normalize=normalize))
end

"""
    pauli_expectation(rho, op, start; normalize=true)

Spin-1/2 case of [`operator_expectation`](@ref): evaluate `tr(ρ O) / tr(ρ)` for one dense local
operator `op` applied at `start` against a vectorized operator `rho` on Pauli sites.

# Arguments
- `rho`: Operator-space `MPS` on dimension-4 Pauli sites.
- `op`: Dense local Pauli-space operator.
- `start`: Left-edge site index of `op`'s support.

# Keyword Arguments
- `normalize`: See [`operator_expectation_profile`](@ref).

# Returns
- A `ComplexF64` expectation value (real up to numerical noise for Hermitian inputs).
"""
pauli_expectation(rho::MPS, op::AbstractMatrix, start::Integer; normalize::Bool=true) =
  operator_expectation(rho, op, start; normalize=normalize)
