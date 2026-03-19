"""
    LocalGateEvolution

Configuration for dense local-gate evolution on a finite `MPS`.

# Fields
- `gate`: One dense gate, a per-step gate vector, or a callable gate provider.
- `dt`: Logical time increment associated with one call to [`evolve!`](@ref).
- `schedule`: Bond schedule that determines where each gate application acts.
- `nstep`: Number of complete passes through `schedule` per `evolve!` call.
- `maxdim`: Maximum bond dimension passed to the underlying ITensor gate application.
- `cutoff`: Singular-value cutoff passed to the underlying ITensor gate application.

# Notes
- `LocalGateEvolution` does not hide the schedule. A caller can inspect or modify it
  directly, which is useful for explicit TEBD-style workflows.
"""
struct LocalGateEvolution{TG,TS}
  gate::TG
  dt::Float64
  schedule::TS
  nstep::Int
  maxdim::Int
  cutoff::Float64
end

"""
    LocalGateEvolution(gate, dt; schedule=nothing, nstep=1, maxdim=0, cutoff=0.0)

Construct a [`LocalGateEvolution`](@ref).

# Arguments
- `gate`: Dense gate specification. This may be:
  - one matrix reused at every schedule entry
  - a vector of matrices indexed by schedule position
  - a callable `(bond, index) -> gate`
- `dt`: Logical time increment associated with one `evolve!` call.

# Keyword Arguments
- `schedule`: Ordered collection of bond specifiers. The TEBD helpers fill this in
  automatically.
- `nstep`: Number of complete schedule traversals per `evolve!` call.
- `maxdim`: Maximum bond dimension for gate application.
- `cutoff`: Truncation cutoff for gate application.

# Returns
- A concrete `LocalGateEvolution` object with normalized numeric field types.
"""
function LocalGateEvolution(gate, dt; schedule=nothing, nstep=1, maxdim=0, cutoff=0.0)
  return LocalGateEvolution(gate, Float64(dt), schedule, Int(nstep), Int(maxdim), Float64(cutoff))
end

"""
    DMTGateEvolution

Configuration for scheduled operator-space DMT evolution.

# Fields
- `gate`: Dense local gate specification in the Pauli basis.
- `dt`: Logical time increment associated with one full DMT evolution call.
- `schedule`: Forward update schedule.
- `reverse_schedule`: Reverse update schedule used for the backward sweep.
- `nstep`: Number of complete forward-plus-reverse sweeps per `evolve!` call.
- `maxdim`: Target post-DMT bond dimension.
- `cutoff`: Truncation cutoff used in the final repair SVD.
- `gate_maxdim`: Temporary bond dimension budget used for raw gate application.
- `connector_buffer`: Number of connector modes protected during reduced-matrix truncation.
"""
struct DMTGateEvolution{TG,TS,TR}
  gate::TG
  dt::Float64
  schedule::TS
  reverse_schedule::TR
  nstep::Int
  maxdim::Int
  cutoff::Float64
  gate_maxdim::Int
  connector_buffer::Int
end

"""
    DMTGateEvolution(gate, dt; schedule, reverse_schedule=reverse(schedule), nstep=1, maxdim=30, cutoff=1e-12, gate_maxdim=max(maxdim * 16, 64), connector_buffer=8)

Construct a [`DMTGateEvolution`](@ref).

# Arguments
- `gate`: Dense local gate specification in the operator-space Pauli basis.
- `dt`: Logical time increment associated with one `dmt_evolve!` call.

# Keyword Arguments
- `schedule`: Forward update schedule.
- `reverse_schedule`: Reverse update schedule. By default the forward schedule is reversed.
- `nstep`: Number of complete forward-plus-reverse sweeps per evolution call.
- `maxdim`: Target bond dimension after DMT truncation.
- `cutoff`: Truncation cutoff used when repairing the compressed bond.
- `gate_maxdim`: Temporary gate-application bond dimension budget.
- `connector_buffer`: Number of connector directions preserved before reduced-matrix
  truncation is applied.

# Returns
- A concrete `DMTGateEvolution` object with normalized numeric field types.
"""
function DMTGateEvolution(
  gate,
  dt;
  schedule,
  reverse_schedule=reverse(schedule),
  nstep=1,
  maxdim=30,
  cutoff=1e-12,
  gate_maxdim=max(Int(maxdim) * 16, 64),
  connector_buffer=8,
)
  return DMTGateEvolution(
    gate,
    Float64(dt),
    schedule,
    reverse_schedule,
    Int(nstep),
    Int(maxdim),
    Float64(cutoff),
    Int(gate_maxdim),
    Int(connector_buffer),
  )
end

"""
    TDVPEvolution

Configuration for finite-`MPS` TDVP evolution driven by an MPO generator.

# Fields
- `generator`: MPO or compatible object passed to `tdvp`.
- `t`: Total evolution interval passed to `tdvp`.
- `time_step`: Internal TDVP time step.
- `nsteps`: Preferred number of TDVP steps.
- `nsweeps`: Legacy fallback step-count field used when `nsteps` is `nothing`.
- `reverse_step`: Whether TDVP should alternate sweep direction.
- `updater_backend`: Backend name forwarded to `tdvp`.
- `updater`: Optional custom TDVP updater.
- `normalize`: Whether TDVP should normalize the state.
- `solver_kwargs`: Additional keyword arguments forwarded to `tdvp`.
- `schedule`: Optional metadata field kept for interface symmetry with TEBD-based code.

# Notes
- `tdvp_evolve!` interprets the effective step count as `nsteps` if present, otherwise
  `nsweeps`.
"""
struct TDVPEvolution{TH,TT,TS,TB,TU,TK}
  generator::TH
  t::TT
  time_step
  nsteps::Union{Nothing,Int}
  nsweeps::Union{Nothing,Int}
  reverse_step::Bool
  updater_backend::TB
  updater::TU
  normalize::Bool
  solver_kwargs::TK
  schedule::TS
end

"""
    TDVPEvolution(generator, t; kwargs...)

Construct a [`TDVPEvolution`](@ref) for finite OBC `MPS` states.

# Arguments
- `generator`: MPO-style generator passed to `tdvp`.
- `t`: Total evolution interval passed to `tdvp`.

# Keyword Arguments
- `time_step`: Internal TDVP time step.
- `nsteps`: Preferred number of internal TDVP steps.
- `nsweeps`: Fallback internal step count used if `nsteps` is `nothing`.
- `reverse_step`: Whether to alternate left-to-right and right-to-left sweeps.
- `updater_backend`: Backend string forwarded to `tdvp`.
- `updater`: Optional custom updater callback.
- `normalize`: Whether to normalize the state after evolution.
- `solver_kwargs`: Additional keyword arguments forwarded to `tdvp`.
- `schedule`: Optional metadata used by higher-level workflows.

# Returns
- A concrete `TDVPEvolution` object.
"""
function TDVPEvolution(
  generator,
  t;
  time_step=nothing,
  nsteps=nothing,
  nsweeps=nothing,
  reverse_step=true,
  updater_backend="exponentiate",
  updater=nothing,
  normalize=false,
  solver_kwargs=(;),
  schedule=nothing,
)
  return TDVPEvolution(
    generator,
    t,
    time_step,
    isnothing(nsteps) ? nothing : Int(nsteps),
    isnothing(nsweeps) ? nothing : Int(nsweeps),
    Bool(reverse_step),
    updater_backend,
    updater,
    Bool(normalize),
    solver_kwargs,
    schedule,
  )
end
