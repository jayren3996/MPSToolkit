"""
    LocalGateEvolution

Configuration for local-gate evolution driven by a dense gate matrix over an explicit bond schedule.
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

Construct a local-gate evolution configuration with explicit schedule and truncation controls.
"""
function LocalGateEvolution(gate, dt; schedule=nothing, nstep=1, maxdim=0, cutoff=0.0)
  return LocalGateEvolution(gate, Float64(dt), schedule, Int(nstep), Int(maxdim), Float64(cutoff))
end

"""
    DMTGateEvolution

Configuration for scheduled operator-space DMT driven by a dense gate matrix over explicit forward and reverse schedules.
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

Construct a scheduled operator-space DMT configuration sharing the same gate/schedule structure as local-gate TEBD.
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

Configuration for finite-MPS TDVP evolution driven by an MPO generator and explicit solver options.
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

Construct a TDVP evolution configuration for finite OBC `MPS` states.
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
