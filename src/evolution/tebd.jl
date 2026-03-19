"""
    _bond_start(bond)

Normalize a scheduled finite-`MPS` bond specifier to its left-site index.

# Arguments
- `bond`: Bond identifier stored in a TEBD schedule.

# Returns
- The left-site index at which the local gate should be applied.

# Notes
- The current implementation accepts integer bond labels unchanged, but the helper exists so
  richer schedule encodings can be introduced without rewriting callers.
"""
function _bond_start(bond::Integer)
  return Int(bond)
end

"""
    _default_gate_from_hamiltonian(h, step_dt)

Convert a dense local Hamiltonian into a dense real-time TEBD gate.

# Arguments
- `h`: Dense local Hamiltonian matrix.
- `step_dt`: Time increment attached to that local update.

# Returns
- `exp(-im * step_dt * h)`.
"""
function _default_gate_from_hamiltonian(h::AbstractMatrix, step_dt)
  return exp(-im * step_dt * h)
end

"""
    _dense_local_operator(sites, op)

Convert a dense local operator matrix into an `ITensor` acting on `sites`.

# Arguments
- `sites`: Site indices defining the local Hilbert space ordering.
- `op`: Dense square matrix whose dimension matches `prod(dim.(sites))`.

# Returns
- An `ITensor` with primed output legs and dagged input legs.

# Notes
- This helper is shared by gate application and dense local observable evaluation.
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

# Arguments
- `psi`: Finite `MPS` whose local site dimension defines the base Hilbert-space dimension.
- `op`: Dense square local operator matrix.

# Returns
- The number of sites acted on by `op`.

# Notes
- The calculation assumes a uniform local dimension across the chain.
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

Resolve the concrete dense gate to apply for one TEBD schedule entry.

# Arguments
- `gate_spec`: Gate specification stored inside an evolution object.
- `bond`: Bond label for the current schedule entry.
- `index`: Position of the current entry inside the schedule.

# Returns
- The dense matrix that should be applied for this schedule entry.
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

Return the nearest-neighbor odd-even-odd Strang schedule for a finite OBC chain.

# Arguments
- `nsites`: Number of sites in the chain.

# Returns
- `(schedule, weights)` where:
  - `schedule` is the ordered bond list
  - `weights` stores the Strang half-step / full-step prefactors for each entry

# Notes
- The weights are intended to be passed into a local Hamiltonian builder so that half steps
  on odd bonds and full steps on even bonds are handled explicitly.
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

# Arguments
- `hamiltonians`: Hamiltonian specification. Supported forms match the `gate` conventions of
  [`LocalGateEvolution`](@ref).
- `dt`: Time increment passed to `map_hamiltonian`.

# Keyword Arguments
- `map_hamiltonian`: Function `(h, dt) -> gate`, defaulting to `exp(-im * dt * h)`.

# Returns
- One gate, a vector of gates, or a callable gate provider with the same structural shape
  as `hamiltonians`.
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

Construct a [`LocalGateEvolution`](@ref) from dense local Hamiltonians.

# Arguments
- `hamiltonians`: Hamiltonian specification accepted by [`local_gates_from_hamiltonians`](@ref).
- `dt`: Logical TEBD time step.

# Keyword Arguments
- `map_hamiltonian`: Function used to convert each Hamiltonian into a gate.
- `schedule`: Explicit bond schedule.
- `nstep`: Number of complete schedule passes per `evolve!` call.
- `maxdim`: Bond dimension cap for gate application.
- `cutoff`: Truncation cutoff for gate application.

# Returns
- A `LocalGateEvolution` whose `gate` field contains the converted dense TEBD gates.
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

# Arguments
- `nsites`: Number of sites in the open chain.
- `dt`: Logical time step associated with one Strang sweep.

# Keyword Arguments
- `local_hamiltonian`: Function `(bond, weight) -> h_local` used to build each scheduled
  Hamiltonian term.
- `map_hamiltonian`: Function used to exponentiate each local Hamiltonian.
- `nstep`: Number of full Strang sweeps per `evolve!` call.
- `maxdim`: Bond dimension cap for TEBD gate application.
- `cutoff`: Truncation cutoff for TEBD gate application.

# Returns
- A `LocalGateEvolution` with explicit Strang schedule metadata.
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

# Arguments
- `psi`: Finite matrix-product state to mutate in place.
- `gate`: Dense one-site, two-site, or few-site gate matrix.
- `bond`: Left-edge location of the update window.

# Keyword Arguments
- `maxdim`: Maximum bond dimension allowed during gate application.
- `cutoff`: Truncation cutoff used by ITensor operations.

# Returns
- The mutated `psi`.

# Notes
- Periodic wraparound is currently supported only for two-site gates and only when the
  schedule explicitly uses `bond == length(psi)`.
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

Run scheduled local-gate TEBD evolution on a finite `MPS`.

# Arguments
- `psi`: State to mutate in place.
- `evo`: [`LocalGateEvolution`](@ref) describing gates, schedule, and truncation settings.

# Returns
- The mutated `psi`.

# Notes
- One call runs `evo.nstep` complete traversals of `evo.schedule`.
- This function is the implementation behind the generic [`evolve!`](@ref) dispatch for
  `LocalGateEvolution`.
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
