"""
    pauli_trace(rho)

Return the physical trace of a vectorized operator stored as a Pauli-basis `MPS`.

In the normalized Pauli convention `tr(ρ) = (√2)^N c_{I…I}`, where `c_{I…I}` is the
amplitude of the all-identity string, so the trace is one identity-environment contraction.

# Arguments
- `rho`: Operator-space `MPS` on dimension-4 Pauli sites.

# Returns
- The trace as a `ComplexF64` (real up to numerical noise for Hermitian operators).
"""
function pauli_trace(rho::MPS)
  _validate_pauli_operator_space(rho, 1, length(rho))
  return 2.0^(length(rho) / 2) * scalar(_right_identity_environment(rho, 1))
end

"""
    _pauli_window_cap(window_sites, op)

Build the window cap tensor whose entry at the multi-site Pauli label `α` is `tr(P_α O)`,
so that contracting it against an operator-space MPS window (with identity caps elsewhere)
yields `tr(ρ O)` up to the `(√2)` factors handled by the callers.
"""
function _pauli_window_cap(window_sites, op::AbstractMatrix;
                          local_basis=[matrix / sqrt(2) for matrix in values(pauli_matrices())])
  span = length(window_sites)
  size(op) == (2^span, 2^span) || throw(ArgumentError("window cap operator size must match the window span"))
  cap = ITensor(ComplexF64, window_sites...)
  for labels in Iterators.product(ntuple(_ -> 1:4, span)...)
    pauli = foldl(kron, (local_basis[l] for l in labels))
    value = tr(pauli * op)
    iszero(value) && continue
    cap[(site => label for (site, label) in zip(window_sites, labels))...] = value
  end
  return cap
end

function _validated_pauli_windows(rho::MPS, terms)
  windows = Vector{Tuple{Int,Int}}(undef, length(terms))
  for (k, (start, op)) in enumerate(terms)
    size(op, 1) == size(op, 2) || throw(ArgumentError("dense local operator must be square"))
    span = _spinhalf_span(size(op, 1))
    start >= 1 || throw(ArgumentError("term start must be at least 1"))
    start + span - 1 <= length(rho) || throw(ArgumentError("term support exceeds chain length"))
    windows[k] = (Int(start), span)
  end
  return windows
end

"""
    pauli_expectation_profile(rho, terms; normalize=true)

Evaluate `tr(ρ O_k)` (optionally over `tr(ρ)`) for a list of dense local operators against a
vectorized operator `rho`, in one O(N) sweep with cumulative identity environments.

# Arguments
- `rho`: Operator-space `MPS` on dimension-4 Pauli sites.
- `terms`: Collection of `(start, op)` pairs; each `op` is a dense physical operator on
  `2^span` dimensions acting on the sites `start:(start + span - 1)`. The pairs may be given
  in any order; results follow the input order.

# Keyword Arguments
- `normalize`: If `true` (default), return `tr(ρ O_k) / tr(ρ)`; otherwise return the
  unnormalized `tr(ρ O_k)`.

# Returns
- A `Vector{ComplexF64}` of expectation values. For Hermitian `ρ` and `O_k` the entries are
  real up to numerical noise.

# Notes
- This is the energy-density-profile measurement for operator-space transport runs: pass the
  Hamiltonian terms (e.g. from [`pxp_term_hamiltonian`](@ref)) and take the real part.
"""
function pauli_expectation_profile(rho::MPS, terms; normalize::Bool=true)
  nsites = length(rho)
  _validate_pauli_operator_space(rho, 1, nsites)
  isempty(terms) && return ComplexF64[]
  windows = _validated_pauli_windows(rho, terms)
  local_basis = [matrix / sqrt(2) for matrix in values(pauli_matrices())]

  right = Vector{ITensor}(undef, nsites + 1)
  right[nsites + 1] = ITensor(1.0)
  for site in nsites:-1:1
    right[site] = right[site + 1] * (_pauli_identity_env(siteind(rho, site)) * rho[site])
  end
  denominator = scalar(right[1])
  normalize && iszero(denominator) && throw(ArgumentError("normalized Pauli expectations require a nonzero trace"))

  results = Vector{ComplexF64}(undef, length(terms))
  left = ITensor(1.0)
  absorbed = 0
  for k in sortperm(windows; by=first)
    start, span = windows[k]
    while absorbed < start - 1
      absorbed += 1
      left = left * (_pauli_identity_env(siteind(rho, absorbed)) * rho[absorbed])
    end
    window_block = rho[start]
    for site in (start + 1):(start + span - 1)
      window_block *= rho[site]
    end
    cap = _pauli_window_cap([siteind(rho, s) for s in start:(start + span - 1)], last(terms[k]); local_basis=local_basis)
    raw = scalar(left * (cap * window_block) * right[start + span])
    results[k] = normalize ? raw / (2.0^(span / 2) * denominator) : 2.0^((nsites - span) / 2) * raw
  end
  return results
end

"""
    pauli_expectation(rho, op, start; normalize=true)

Evaluate `tr(ρ O) / tr(ρ)` (or the unnormalized `tr(ρ O)` with `normalize=false`) for one
dense local operator `op` applied at `start` against a vectorized operator `rho`.

See [`pauli_expectation_profile`](@ref) for the batched O(N) version and the conventions.

# Returns
- A `ComplexF64` expectation value (real up to numerical noise for Hermitian inputs).
"""
function pauli_expectation(rho::MPS, op::AbstractMatrix, start::Integer; normalize::Bool=true)
  return only(pauli_expectation_profile(rho, [(Int(start), op)]; normalize=normalize))
end
