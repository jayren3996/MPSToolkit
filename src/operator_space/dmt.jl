"""
    DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=480, connector_buffer=8)

Options controlling operator-space density matrix truncation.
"""
Base.@kwdef struct DMTOptions
  maxdim::Int = 30
  cutoff::Float64 = 1e-12
  gate_maxdim::Int = 480
  connector_buffer::Int = 8
end

function _pauli_identity_env(site)
  tensor = ITensor(site)
  tensor[site => 1] = 1.0
  return tensor
end

function _left_identity_environment(psi::MPS, stop::Integer)
  env = ITensor(1.0)
  for site in 1:Int(stop)
    env *= _pauli_identity_env(siteind(psi, site)) * psi[site]
  end
  return env
end

function _right_identity_environment(psi::MPS, start::Integer)
  env = ITensor(1.0)
  for site in length(psi):-1:Int(start)
    env *= _pauli_identity_env(siteind(psi, site)) * psi[site]
  end
  return env
end

function _dmt_truncation_bond(start::Integer, span::Integer, direction::Symbol)
  direction === :R && return Int(start)
  direction === :L && return Int(start + span - 2)
  throw(ArgumentError("DMT direction must be :R or :L"))
end

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

function _dmt_bond_truncate!(psi::MPS, bond::Integer; maxdim::Integer, cutoff::Real, direction::Symbol=:R, connector_buffer::Integer=8)
  maxdim > 0 || return psi
  1 <= bond < length(psi) || throw(ArgumentError("DMT bond must lie in 1:length(psi)-1"))
  current_link = linkind(psi, bond)
  isnothing(current_link) && return psi
  dim(current_link) <= maxdim && return psi

  orthogonalize!(psi, bond)
  left_site = siteind(psi, bond)
  right_site = siteind(psi, bond + 1)

  left_env = _left_identity_environment(psi, bond - 1)
  right_env = _right_identity_environment(psi, bond + 2)

  u, s, v = svd(psi[bond], (linkind(psi, bond - 1), left_site))
  psi[bond] = u
  psi[bond + 1] = v * psi[bond + 1]

  left_link = commonind(u, s)
  right_link = commonind(v, s)
  left_basis = Matrix(qr(matrix(left_env * psi[bond], left_link, left_site)).Q)
  right_basis = Matrix(qr(matrix(psi[bond + 1] * right_env, right_link, right_site)).Q)
  singular_values = Matrix(matrix(s, left_link, right_link))

  reduced = transpose(left_basis) * singular_values * right_basis
  _mat_trunc!(reduced, maxdim - connector_buffer; connector_buffer=connector_buffer)

  repaired = ITensor(left_basis * reduced * transpose(right_basis), left_link, right_link)
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

"""
    dmt_step!(psi, gate, bond; maxdim=30, cutoff=1e-12, direction=:R, gate_maxdim=max(maxdim * 16, 64), connector_buffer=8)

Apply one local operator-space gate and then perform DMT-preserving truncation on the associated bond.
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
  span = _operator_span(psi, gate)
  tebd_evolve!(psi, gate, start; maxdim=Int(gate_maxdim), cutoff=0.0)
  truncation_bond = _dmt_truncation_bond(start, span, direction)
  _dmt_bond_truncate!(
    psi,
    truncation_bond;
    maxdim=Int(maxdim),
    cutoff=cutoff,
    direction=direction,
    connector_buffer=Int(connector_buffer),
  )
  return psi
end

"""
    dmt_evolve!(psi, evo::DMTGateEvolution)

Run scheduled operator-space DMT evolution using the same gate/schedule structure as local-gate TEBD.
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
      local_gate = _gate_for_step(evo.gate, bond, index)
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

function evolve!(psi::MPS, evo::DMTGateEvolution)
  return dmt_evolve!(psi, evo)
end
