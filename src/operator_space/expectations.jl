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
  machine precision (measured relative difference at most 1.7e-16 for `d` in `2:5` at spans
  1-3, and for `d = 2` at spans 4-5).
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
  immune.

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
  isempty(terms) && return ComplexF64[]
  windows = _validated_operator_windows(rho, terms, d)

  right = Vector{ITensor}(undef, nsites + 1)
  right[nsites + 1] = ITensor(1.0)
  for site in nsites:-1:1
    right[site] = right[site + 1] * (_identity_env(siteind(rho, site)) * rho[site])
  end
  denominator = scalar(right[1])
  # Reject a numerically-negligible (not just exactly-zero) trace relative to the operator
  # scale. For a traceless operator the post-truncation trace is an O(eps) residue rather than
  # exactly zero, and normalizing by it silently amplifies every entry by ~1/eps. This mirrors
  # the relative tolerance the DMT kernel (`_dmt_connector`) already uses; traceless operators
  # should be measured with `normalize=false`.
  normalize && abs(denominator) <= sqrt(eps(Float64)) * norm(rho) &&
    throw(ArgumentError("normalized operator-space expectations require a nonzero trace; the operator trace is numerically negligible relative to its norm (use normalize=false for a traceless operator)"))

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
