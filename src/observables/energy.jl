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

Estimate the finite-chain average energy density of an `MPS` for a dense local operator.

# Arguments
- `psi`: Finite matrix-product state.
- `op`: Dense local operator acting on `span` consecutive sites.

# Keyword Arguments
- `span`: Explicit support size of `op`. If omitted, the support is inferred from the local
  Hilbert-space dimension of `psi`.

# Returns
- The arithmetic mean of the local expectation value of `op` over all valid windows.

# Notes
- No translation invariance is assumed; the routine explicitly averages over all positions.
"""
function energy_density(psi::MPS, op::AbstractMatrix; span::Int=_operator_span(psi, op))
  last_start = length(psi) - span + 1
  values = Float64[]
  for start in 1:last_start
    sites = [siteind(psi, n) for n in start:(start + span - 1)]
    op_tensor = _dense_local_operator(sites, op)
    theta = _finite_local_wavefunction(psi, start, span)
    push!(values, real(inner(theta, apply(op_tensor, theta))))
  end
  return sum(values) / length(values)
end

"""
    energy_density(psi, op::MPO)

Return the finite-chain average energy density of an `MPS` for an `MPO`.

# Arguments
- `psi`: Finite matrix-product state.
- `op`: Hamiltonian or observable represented as an `MPO`.

# Returns
- `real(<psi|op|psi>) / length(psi)`.
"""
function energy_density(psi::MPS, op::MPO)
  return real(inner(psi', op, psi)) / length(psi)
end
