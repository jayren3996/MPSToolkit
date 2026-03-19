"""
    DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=480, connector_buffer=8)

Options controlling operator-space density matrix truncation.

# Fields
- `maxdim`: Target bond dimension after DMT truncation.
- `cutoff`: Truncation cutoff used in the repair SVD.
- `gate_maxdim`: Temporary bond dimension budget used during raw gate application.
- `connector_buffer`: Number of connector directions protected before reduced-matrix
  truncation.
"""
Base.@kwdef struct DMTOptions
  maxdim::Int = 30
  cutoff::Float64 = 1e-12
  gate_maxdim::Int = 480
  connector_buffer::Int = 8
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
    _dmt_truncation_bond(start, span, direction)

Return the bond index at which DMT truncation should be applied for one local update.
"""
function _dmt_truncation_bond(start::Integer, span::Integer, direction::Symbol)
  direction === :R && return Int(start)
  direction === :L && return Int(start + span - 2)
  throw(ArgumentError("DMT direction must be :R or :L"))
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

# Returns
- The mutated `psi`.
"""
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

Apply one local operator-space gate and then perform DMT-preserving truncation.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `gate`: Dense local gate in the Pauli basis.
- `bond`: Left-edge location of the local update.

# Keyword Arguments
- `maxdim`: Target post-truncation bond dimension.
- `cutoff`: Truncation cutoff used in the final repair SVD.
- `direction`: Sweep direction, either `:R` or `:L`.
- `gate_maxdim`: Temporary bond dimension budget used for the raw gate application.
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

Run scheduled operator-space DMT evolution.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `evo`: [`DMTGateEvolution`](@ref) describing the gate specification, schedules, and
  truncation budgets.

# Returns
- The mutated and normalized `psi`.

# Notes
- One call runs `evo.nstep` complete forward-and-reverse sweeps.
- `normalize!(psi)` is applied at the end because operator-space trajectories typically
  interpret the state as a normalized vectorized operator.
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

"""
    evolve!(psi, evo::DMTGateEvolution)

Dispatch operator-space evolution through the DMT backend.
"""
function evolve!(psi::MPS, evo::DMTGateEvolution)
  return dmt_evolve!(psi, evo)
end
