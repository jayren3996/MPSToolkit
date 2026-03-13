"""
    _bond_start(bond)

Normalize a scheduled finite-MPS bond specifier to its left-site index.
"""
function _bond_start(bond::Integer)
  return Int(bond)
end

"""
    _default_gate_from_hamiltonian(h, step_dt)

Convert a dense local Hamiltonian into a dense TEBD gate over step `step_dt`.
"""
function _default_gate_from_hamiltonian(h::AbstractMatrix, step_dt)
  return exp(-im * step_dt * h)
end

"""
    _dense_local_operator(sites, op)

Convert a dense local operator matrix into an ITensor acting on `sites`.
"""
function _dense_local_operator(sites, op::AbstractMatrix)
  dims = dim.(sites)
  prod(dims) == size(op, 1) || throw(ArgumentError("operator dimension does not match site dimensions"))
  size(op, 1) == size(op, 2) || throw(ArgumentError("dense local operator must be square"))
  tensor_data = reshape(op, Tuple(vcat(dims, dims)))
  return itensor(tensor_data, prime.(sites)..., dag.(sites)...)
end

"""
    _operator_span(psi, op)

Infer the support size of a dense local operator from the site dimension of `psi`.
"""
function _operator_span(psi::MPS, op::AbstractMatrix)
  size(op, 1) == size(op, 2) || throw(ArgumentError("dense local operator must be square"))
  d = dim(siteind(psi, 1))
  d > 1 || throw(ArgumentError("site dimension must be larger than one"))
  span = round(Int, log(size(op, 1)) / log(d))
  d^span == size(op, 1) || throw(ArgumentError("operator size is incompatible with the site dimension"))
  return span
end

"""
    _gate_for_step(gate_spec, bond, index)

Resolve the concrete dense gate to apply for one finite-TEBD schedule entry.
"""
function _gate_for_step(gate_spec::AbstractMatrix, bond, index)
  return gate_spec
end

"""
    _gate_for_step(gate_spec, bond, index)

Resolve a per-entry gate from an indexable gate collection.
"""
function _gate_for_step(gate_spec::AbstractVector, bond, index)
  return gate_spec[index]
end

"""
    _gate_for_step(gate_spec, bond, index)

Resolve a gate from a callable gate provider.
"""
function _gate_for_step(gate_spec::Function, bond, index)
  return gate_spec(bond, index)
end

"""
    tebd_strang_schedule(nsites)

Return the nearest-neighbor odd-even-odd Strang schedule and per-entry weights for a finite OBC chain.
"""
function tebd_strang_schedule(nsites::Integer)
  nsites < 2 && throw(ArgumentError("Strang TEBD requires at least two sites"))
  odd_bonds = collect(1:2:(nsites - 1))
  even_bonds = collect(2:2:(nsites - 1))
  schedule = vcat(odd_bonds, even_bonds, odd_bonds)
  weights = vcat(fill(0.5, length(odd_bonds)), fill(1.0, length(even_bonds)), fill(0.5, length(odd_bonds)))
  return schedule, weights
end

"""
    local_gates_from_hamiltonians(hamiltonians, dt; map_hamiltonian=_default_gate_from_hamiltonian)

Convert dense local Hamiltonian data into dense TEBD gates.
`hamiltonians` may be one dense matrix, a per-step vector of matrices, or a callable provider.
"""
function local_gates_from_hamiltonians(hamiltonians::AbstractMatrix, dt; map_hamiltonian::Function=_default_gate_from_hamiltonian)
  return map_hamiltonian(hamiltonians, dt)
end

function local_gates_from_hamiltonians(hamiltonians::AbstractVector, dt; map_hamiltonian::Function=_default_gate_from_hamiltonian)
  return [map_hamiltonian(h, dt) for h in hamiltonians]
end

function local_gates_from_hamiltonians(hamiltonians::Function, dt; map_hamiltonian::Function=_default_gate_from_hamiltonian)
  return (bond, index) -> map_hamiltonian(hamiltonians(bond, index), dt)
end

"""
    tebd_evolution_from_hamiltonians(hamiltonians, dt; map_hamiltonian=_default_gate_from_hamiltonian, schedule=nothing, nstep=1, maxdim=0, cutoff=0.0)

Construct a `LocalGateEvolution` by converting dense local Hamiltonians into TEBD gates.
"""
function tebd_evolution_from_hamiltonians(
  hamiltonians,
  dt;
  map_hamiltonian::Function=_default_gate_from_hamiltonian,
  schedule=nothing,
  nstep=1,
  maxdim=0,
  cutoff=0.0,
)
  gates = local_gates_from_hamiltonians(hamiltonians, dt; map_hamiltonian=map_hamiltonian)
  return LocalGateEvolution(gates, dt; schedule=schedule, nstep=nstep, maxdim=maxdim, cutoff=cutoff)
end

"""
    tebd_strang_evolution(nsites, dt; local_hamiltonian, map_hamiltonian=_default_gate_from_hamiltonian, nstep=1, maxdim=0, cutoff=0.0)

Construct a nearest-neighbor odd-even-odd Strang `LocalGateEvolution` from a local Hamiltonian builder.
`local_hamiltonian(bond, weight)` must return the dense local Hamiltonian for one scheduled update.
"""
function tebd_strang_evolution(
  nsites::Integer,
  dt;
  local_hamiltonian::Function,
  map_hamiltonian::Function=_default_gate_from_hamiltonian,
  nstep=1,
  maxdim=0,
  cutoff=0.0,
)
  schedule, weights = tebd_strang_schedule(nsites)
  hamiltonians = [local_hamiltonian(bond, weight) for (bond, weight) in zip(schedule, weights)]
  return tebd_evolution_from_hamiltonians(
    hamiltonians,
    dt;
    map_hamiltonian=map_hamiltonian,
    schedule=schedule,
    nstep=nstep,
    maxdim=maxdim,
    cutoff=cutoff,
  )
end

"""
    tebd_evolve!(psi, gate, bond; maxdim, cutoff)

Apply one dense local TEBD gate to a finite MPS block starting at `bond`.
If `bond == length(psi)` and the gate is two-site, the update is applied across the periodic
boundary between the last and first sites.
"""
function tebd_evolve!(psi::MPS, gate::AbstractMatrix, bond; maxdim::Int, cutoff::Real)
  n = _bond_start(bond)
  span = _operator_span(psi, gate)
  if span == 1
    n <= length(psi) || throw(ArgumentError("local gate support exceeds chain length"))
    sites = [siteind(psi, n)]
    gate_tensor = _dense_local_operator(sites, gate)
    updated = product(gate_tensor, psi, [n]; maxdim=maxdim, cutoff=cutoff)
    psi[:] = updated
    return psi
  end
  if n == length(psi)
    span == 2 || throw(ArgumentError("periodic boundary TEBD currently supports only two-site gates"))
    reordered = movesite(psi, 1 => length(psi); orthocenter=length(psi), maxdim=maxdim, cutoff=cutoff)
    sites = [siteind(reordered, length(psi) - 1), siteind(reordered, length(psi))]
    gate_tensor = _dense_local_operator(sites, gate)
    updated = product(gate_tensor, reordered, [length(psi) - 1, length(psi)]; maxdim=maxdim, cutoff=cutoff)
    restored = movesite(updated, length(psi) => 1; orthocenter=1, maxdim=maxdim, cutoff=cutoff)
    psi[:] = restored
    return psi
  end
  last_site = n + span - 1
  last_site <= length(psi) || throw(ArgumentError("local gate support exceeds chain length"))
  sites = [siteind(psi, j) for j in n:last_site]
  gate_tensor = _dense_local_operator(sites, gate)
  updated = product(gate_tensor, psi, collect(n:last_site); maxdim=maxdim, cutoff=cutoff)
  psi[:] = updated
  return psi
end

"""
    evolve!(psi, evo::LocalGateEvolution)

Run scheduled local-gate TEBD evolution on a finite OBC MPS.
"""
function evolve!(psi::MPS, evo::LocalGateEvolution)
  for _ in 1:evo.nstep
    for (index, bond) in pairs(evo.schedule)
      local_gate = _gate_for_step(evo.gate, bond, index)
      tebd_evolve!(psi, local_gate, bond; maxdim=evo.maxdim, cutoff=evo.cutoff)
    end
  end
  return psi
end
