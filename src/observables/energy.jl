"""
    _finite_local_wavefunction(psi, start, span)

Assemble a local wavefunction tensor for a finite `MPS`.

# Arguments
- `psi`: Finite matrix-product state.
- `start`: First site in the local window.
- `span`: Number of consecutive sites to include.

# Returns
- A local `ITensor` representing the contracted wavefunction on the requested window.

# Notes
- The input `psi` is orthogonalized around `start` before the local block is contracted.
"""
function _finite_local_wavefunction(psi::MPS, start::Int, span::Int)
  centered = orthogonalize(psi, start)
  theta = centered[start]
  for n in (start + 1):(start + span - 1)
    theta *= centered[n]
  end
  return theta
end

"""
    energy_density(psi, op; span=_operator_span(psi, op))

Estimate the finite-chain energy per site of an `MPS` for a dense local operator.

# Arguments
- `psi`: Finite matrix-product state.
- `op`: Dense local operator acting on `span` consecutive sites.

# Keyword Arguments
- `span`: Explicit support size of `op`. If omitted, the support is inferred from the local
  Hilbert-space dimension of `psi`.

# Returns
- The sum of local expectation values over valid windows divided by `length(psi)`.

# Notes
- No translation invariance is assumed; the routine explicitly sums over all valid open-chain
  positions and normalizes by the number of sites to match the `MPO` overload.
- The returned value is norm-invariant: expectations are divided by `<psi|psi>`.
"""
function energy_density(psi::MPS, op::AbstractMatrix; span::Int=_operator_span(psi, op))
  span > 0 || throw(ArgumentError("operator span must be positive"))
  span <= length(psi) || throw(ArgumentError("operator span exceeds chain length"))
  norm2 = real(inner(psi, psi))
  norm2 > 0 || throw(ArgumentError("energy_density requires a nonzero-norm MPS"))
  last_start = length(psi) - span + 1
  # Sweep the orthogonality center along a single working copy rather than re-orthogonalizing
  # (and copying) the whole MPS at every window. Visiting starts in increasing order moves the
  # center one site at a time, so the cost is O(N) center moves instead of O(N^2). With the
  # center at `start` the block theta = work[start..start+span-1] carries the full norm (the rest
  # is left/right orthonormal), so <theta|op|theta> is the windowed expectation -- identical to
  # the per-window orthogonalize the previous implementation performed.
  work = orthogonalize(psi, 1)
  total = 0.0
  for start in 1:last_start
    orthogonalize!(work, start)
    sites = [siteind(work, n) for n in start:(start + span - 1)]
    op_tensor = _dense_local_operator(sites, op)
    theta = work[start]
    for n in (start + 1):(start + span - 1)
      theta *= work[n]
    end
    total += real(inner(theta, apply(op_tensor, theta)))
  end
  return total / (norm2 * length(psi))
end

"""
    energy_density(psi, op::MPO)

Return the finite-chain average energy density of an `MPS` for an `MPO`.

# Arguments
- `psi`: Finite matrix-product state.
- `op`: Hamiltonian or observable represented as an `MPO`.

# Returns
- `real(<psi|op|psi>) / (<psi|psi> * length(psi))`.
"""
function energy_density(psi::MPS, op::MPO)
  norm2 = real(inner(psi, psi))
  norm2 > 0 || throw(ArgumentError("energy_density requires a nonzero-norm MPS"))
  return real(inner(psi', op, psi)) / (norm2 * length(psi))
end
