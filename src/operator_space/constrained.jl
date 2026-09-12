"""
    ConstrainedDMTRunState(; completed_steps=0, physical_time=0.0,
                          steps_since_projection=0, projection_count=0)

Persistent cadence and clock state for [`constrained_dmt_evolve!`](@ref). Pass the same object
to successive calls, and store it in a DMT checkpoint, so splitting a run or restarting it does
not move constraint projections to different physical steps.
"""
mutable struct ConstrainedDMTRunState
  completed_steps::Int
  physical_time::Float64
  steps_since_projection::Int
  projection_count::Int

  function ConstrainedDMTRunState(completed_steps, physical_time, steps_since_projection,
                                  projection_count)
    completed_steps >= 0 || throw(ArgumentError("completed_steps must be nonnegative"))
    physical_time >= 0 || throw(ArgumentError("physical_time must be nonnegative"))
    steps_since_projection >= 0 ||
      throw(ArgumentError("steps_since_projection must be nonnegative"))
    projection_count >= 0 || throw(ArgumentError("projection_count must be nonnegative"))
    return new(Int(completed_steps), Float64(physical_time), Int(steps_since_projection),
               Int(projection_count))
  end
end

ConstrainedDMTRunState(; completed_steps=0, physical_time=0.0,
                       steps_since_projection=0, projection_count=0) =
  ConstrainedDMTRunState(completed_steps, physical_time, steps_since_projection,
                         projection_count)

function _copy_dmt_evolution(evo::DMTGateEvolution; nstep::Integer=evo.nstep)
  return DMTGateEvolution(
    evo.gate,
    evo.dt;
    schedule=evo.schedule,
    reverse_schedule=evo.reverse_schedule,
    nstep=nstep,
    maxdim=evo.maxdim,
    cutoff=evo.cutoff,
    gate_maxdim=evo.gate_maxdim,
    preserve_diameter=evo.preserve_diameter,
    preserve_operators=evo.preserve_operators,
    truncation=evo.truncation,
    gate_backend=evo.gate_backend,
    normalize=evo.normalize,
  )
end

function _constraint_project!(rho::MPS, projector::MPO, projector_maxdim::Integer,
                              projector_cutoff::Real, normalize::Bool)
  rho[:] = apply(projector, rho; maxdim=Int(projector_maxdim), cutoff=projector_cutoff)
  normalize && normalize!(rho)
  return rho
end

"""
    constraint_leakage_squared(rho, projector)

Return `1 - real(<rho, P rho>)/<rho,rho>` by direct MPS-MPO-MPS contraction. If `P` is an
orthogonal projector, this is the squared relative Hilbert--Schmidt leakage. No compressed
projected state is formed, so the diagnostic has no projection-application truncation error.
"""
function constraint_leakage_squared(rho::MPS, projector::MPO)
  length(projector) == length(rho) ||
    throw(ArgumentError("projector and state must have matching lengths"))
  denominator = real(inner(rho, rho))
  denominator > 0 || throw(ArgumentError("constraint leakage requires a nonzero state"))
  return 1 - real(inner(prime(rho), projector, rho)) / denominator
end

"""
    constrained_dmt_evolve!(rho, evo, projector; project_every=1,
                            projector_maxdim=2*evo.maxdim,
                            projector_cutoff=evo.cutoff, normalize=evo.normalize,
                            run_state=nothing, step_time=2*evo.dt, final_project=false)

Run scheduled operator-space DMT evolution with periodic constraint-projection checkpoints.

!!! warning "Density operators only"
    Like [`dmt_evolve!`](@ref), this is a DMT driver: `rho` must be a near-infinite-temperature
    **density operator** (e.g. a constrained energy domain-wall melt). DMT protects the trace
    component, so it must **not** be used to Heisenberg-evolve a **traceless** operator / two-point
    correlator — use ordinary TEBD for those.

The driver executes the `evo.nstep` forward-plus-reverse DMT sweeps of
[`dmt_evolve!`](@ref). With a `run_state`, projection cadence persists across calls and a short
call does not force a projection. The legacy stateless call keeps its established behavior: it
runs chunks of `project_every` sweeps and projects the final short chunk. For a constrained model
such as PXP the checkpoints remove weight that truncation leaks out of the constrained sector.

# Arguments
- `rho`: Operator-space `MPS` to mutate in place.
- `evo`: [`DMTGateEvolution`](@ref) describing gates, schedules, and truncation budgets.
- `projector`: Operator-space MPO enforcing the constraint, e.g.
  [`pauli_pxp_constraint_projector`](@ref).

# Keyword Arguments
- `project_every`: Number of complete DMT sweeps between checkpoints.
- `projector_maxdim`: Bond-dimension cap for the projector application, default
  `2 * evo.maxdim`. The checkpoint state is typically an `O(leakage)` perturbation of the
  pre-projection state (which DMT just compressed to `evo.maxdim`), and the modest `2x` buffer
  gives the cutoff room to act. This projector compression is a plain SVD and has no DMT
  preservation guarantee; validate it for weak signals or set `projector_cutoff=0`.
  Raising this (in the limit, dropping the cap the way `evo.gate_maxdim = 0` does for the gate)
  hands the compression decision to the next DMT sweep instead, but A/B benchmarks (PXP
  correlator, N=64) show no measurable accuracy gain while the inflated bonds add substantial
  cost per sweep (1.3-2x at moderate `maxdim`, growing with `maxdim` as the cubic SVD scaling
  takes over).
- `projector_cutoff`: Truncation cutoff for the projector application. This is **ITensors'**
  cutoff — discarded sum of squares relative to the total weight — not the DMT one, which bounds
  the discarded complement relative to the leading complement singular value. Defaulting it to
  `evo.cutoff` therefore reuses one number under two meanings that can sit decades apart. The
  projector step is a plain SVD and carries no preservation guarantee, so on a state whose signal
  is `eps` above an infinite-temperature background it can discard weight `eps^2` that the DMT
  sweeps around it are protecting exactly; pass a smaller `projector_cutoff` (or `0.0`) when
  `eps^2` approaches `evo.cutoff`.
- `normalize`: Whether to renormalize at each projection boundary. Stateless calls also let each
  legacy chunk normalize. Defaults to `evo.normalize`, so a
  `DMTGateEvolution(...; normalize=false)` is honored without re-passing the keyword here —
  matching [`dmt_evolve!`](@ref) and [`evolve!`](@ref). Set `false` to preserve absolute
  operator scales, which is required when tracking unnormalized traces of a **traceless**
  operator (e.g. conserved `tr(H O(t))` in correlator-protocol transport runs).
- `run_state`: Optional persistent [`ConstrainedDMTRunState`](@ref). New long-running drivers
  should pass one and save it with the checkpoint.
- `step_time`: Physical time advanced by one complete forward-plus-reverse step. The legacy PXP
  schedule uses `2 * evo.dt`; the clock is advanced from this explicit value rather than from a
  sum of signed gate durations.
- `final_project`: Apply and record one explicit projection after this call, even if the regular
  cadence just projected. This also resets `steps_since_projection`.

# Returns
- The mutated `rho`. Persistent mode normalizes only at recorded projection boundaries, so a
  call ending between boundaries can return a state whose norm is not one.
"""
function constrained_dmt_evolve!(
  rho::MPS,
  evo::DMTGateEvolution,
  projector::MPO;
  project_every::Integer=1,
  projector_maxdim::Integer=2 * evo.maxdim,
  projector_cutoff::Real=evo.cutoff,
  normalize::Bool=evo.normalize,
  run_state::Union{Nothing,ConstrainedDMTRunState}=nothing,
  step_time::Real=2 * evo.dt,
  final_project::Bool=false,
)
  project_every >= 1 || throw(ArgumentError("constrained_dmt_evolve! requires project_every >= 1"))
  projector_maxdim >= 1 || throw(ArgumentError("constrained_dmt_evolve! requires projector_maxdim >= 1"))
  projector_cutoff >= 0 || throw(ArgumentError("constrained_dmt_evolve! requires projector_cutoff >= 0"))
  length(projector) == length(rho) || throw(ArgumentError("projector and state must have matching lengths"))
  step_time >= 0 || throw(ArgumentError("constrained_dmt_evolve! requires step_time >= 0"))
  final_project && isnothing(run_state) && throw(ArgumentError(
    "constrained_dmt_evolve! requires run_state to record final_project"))

  if isnothing(run_state)
    remaining = evo.nstep
    while remaining > 0
      chunk = min(Int(project_every), remaining)
      dmt_evolve!(rho, _copy_dmt_evolution(evo; nstep=chunk); normalize=normalize)
      _constraint_project!(rho, projector, projector_maxdim, projector_cutoff, normalize)
      remaining -= chunk
    end
    return rho
  end

  run_state.steps_since_projection < project_every || throw(ArgumentError(
    "run_state.steps_since_projection must be less than project_every"))
  remaining = evo.nstep
  while remaining > 0
    chunk = min(Int(project_every) - run_state.steps_since_projection, remaining)
    dmt_evolve!(rho, _copy_dmt_evolution(evo; nstep=chunk); normalize=false)
    run_state.completed_steps += chunk
    run_state.physical_time += chunk * Float64(step_time)
    run_state.steps_since_projection += chunk
    remaining -= chunk
    if run_state.steps_since_projection == project_every
      _constraint_project!(rho, projector, projector_maxdim, projector_cutoff, normalize)
      run_state.steps_since_projection = 0
      run_state.projection_count += 1
    end
  end
  if final_project
    _constraint_project!(rho, projector, projector_maxdim, projector_cutoff, normalize)
    run_state.steps_since_projection = 0
    run_state.projection_count += 1
  end
  return rho
end

function constrained_dmt_evolve!(rho::MPS, plan::PXPDMTPlan, projector::MPO;
                                 step_time::Real=plan.tau, kwargs...)
  return constrained_dmt_evolve!(rho, plan.evolution, projector; step_time=step_time,
    kwargs...)
end
