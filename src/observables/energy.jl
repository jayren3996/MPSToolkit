"""
    _finite_local_wavefunction(psi, start, span)

Assemble a local wavefunction tensor for a finite MPS over `span` sites.
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

Estimate the finite-chain average energy density of an MPS for a dense local operator.
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

Return the finite-chain average energy density of an MPS for an MPO Hamiltonian.
"""
function energy_density(psi::MPS, op::MPO)
  return real(inner(psi', op, psi)) / length(psi)
end
