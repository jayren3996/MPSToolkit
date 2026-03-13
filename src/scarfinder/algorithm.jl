"""
    _assign_state!(psi, updated)

Overwrite the storage of `psi` with the contents of `updated`.
"""
function _assign_state!(psi::MPS, updated::MPS)
  psi[:] = updated
  return psi
end

"""
    _clone_state(psi)

Create a working copy of a state for refinement or rollback operations.
"""
function _clone_state(psi)
  return deepcopy(psi)
end

"""
    _match_energy_dense!(psi, evolution, truncation, target)

Apply a dense-local-operator post-evolution correction loop that nudges the state toward a target energy.
"""
function _match_energy_dense!(psi, evolution, truncation, target)
  energy_error = energy_density(psi, target.operator) - target.target
  abs(energy_error) < target.tol && return psi

  correction_dt = target.alpha * energy_error
  correction_gate = exp(-correction_dt * target.operator)
  correction_evolution = LocalGateEvolution(
    correction_gate,
    correction_dt;
    schedule=evolution.schedule,
    nstep=1,
    maxdim=truncation.maxdim,
    cutoff=truncation.cutoff,
  )

  for _ in 1:target.maxstep
    evolve!(psi, correction_evolution)
    project!(psi, truncation)

    next_error = energy_density(psi, target.operator) - target.target
    if energy_error * next_error < 0
      mix = abs(next_error) / (abs(energy_error) + abs(next_error))
      rollback_gate = exp(mix * correction_dt * target.operator)
      rollback_evolution = LocalGateEvolution(
        rollback_gate,
        correction_dt;
        schedule=evolution.schedule,
        nstep=1,
        maxdim=truncation.maxdim,
        cutoff=truncation.cutoff,
      )
      evolve!(psi, rollback_evolution)
      project!(psi, truncation)
      return psi
    elseif abs(next_error) > abs(energy_error)
      rollback_gate = exp(correction_dt * target.operator)
      rollback_evolution = LocalGateEvolution(
        rollback_gate,
        correction_dt;
        schedule=evolution.schedule,
        nstep=1,
        maxdim=truncation.maxdim,
        cutoff=truncation.cutoff,
      )
      evolve!(psi, rollback_evolution)
      project!(psi, truncation)
      return psi
    end
    energy_error = next_error
  end
  return psi
end

"""
    _correction_time(target, energy_error)

Map an energy error to the signed imaginary-time interval used by an energy-correction step.
"""
function _correction_time(target, energy_error)
  return -1im * target.alpha * energy_error
end

"""
    _solver_kwargs(evolution, truncation)

Build TDVP solver keyword settings for an MPO-based energy-correction step.
"""
function _solver_kwargs(evolution::TDVPEvolution, truncation)
  return evolution.solver_kwargs
end

"""
    _solver_kwargs(evolution::TDVPEvolution, truncation)

Return TDVP solver keyword settings for an MPO-based energy-correction step.
"""
function _solver_kwargs(evolution::TDVPEvolution, truncation::BondDimTruncation)
  return evolution.solver_kwargs
end

"""
    _solver_kwargs(evolution, truncation)

Build TDVP solver keyword settings for an MPO-based energy-correction step from a non-TDVP evolution object.
"""
function _solver_kwargs(evolution, truncation::BondDimTruncation)
  maxdim = hasproperty(evolution, :maxdim) ? getproperty(evolution, :maxdim) : truncation.maxdim
  cutoff = hasproperty(evolution, :cutoff) ? getproperty(evolution, :cutoff) : truncation.cutoff
  return (; maxdim=maxdim, cutoff=cutoff)
end

"""
    _tdvp_correction_evolution(operator, correction_time, evolution, truncation)

Construct a finite-MPS TDVP evolution object used for MPO-based energy correction.
"""
function _tdvp_correction_evolution(operator::MPO, correction_time, evolution, truncation)
  return TDVPEvolution(
    operator,
    correction_time;
    time_step=correction_time,
    nsteps=1,
    reverse_step=hasproperty(evolution, :reverse_step) ? getproperty(evolution, :reverse_step) : true,
    updater_backend=hasproperty(evolution, :updater_backend) ? getproperty(evolution, :updater_backend) : "exponentiate",
    updater=hasproperty(evolution, :updater) ? getproperty(evolution, :updater) : nothing,
    normalize=hasproperty(evolution, :normalize) ? getproperty(evolution, :normalize) : false,
    solver_kwargs=_solver_kwargs(evolution, truncation),
  )
end

"""
    _match_energy_mpo!(psi, evolution, truncation, target)

Apply an MPO-based TDVP correction loop that nudges a finite MPS toward a target energy.
"""
function _match_energy_mpo!(psi::MPS, evolution, truncation, target)
  energy_error = energy_density(psi, target.operator) - target.target
  abs(energy_error) < target.tol && return psi

  correction_time = _correction_time(target, energy_error)
  correction_evolution = _tdvp_correction_evolution(target.operator, correction_time, evolution, truncation)

  for _ in 1:target.maxstep
    evolve!(psi, correction_evolution)
    project!(psi, truncation)

    next_error = energy_density(psi, target.operator) - target.target
    abs(next_error) < target.tol && return psi

    if energy_error * next_error < 0 || abs(next_error) > abs(energy_error)
      rollback_evolution = _tdvp_correction_evolution(target.operator, -correction_time, evolution, truncation)
      evolve!(psi, rollback_evolution)
      project!(psi, truncation)
      return psi
    end

    energy_error = next_error
    correction_time = _correction_time(target, energy_error)
    correction_evolution = _tdvp_correction_evolution(target.operator, correction_time, evolution, truncation)
  end
  return psi
end

"""
    match_energy!(psi, evolution, truncation, target)

Apply a post-evolution correction loop that nudges the state toward a target energy.
"""
function match_energy!(psi, evolution, truncation, target)
  isnothing(target.operator) && return psi
  target.operator isa MPO && return _match_energy_mpo!(psi, evolution, truncation, target)
  return _match_energy_dense!(psi, evolution, truncation, target)
end

"""
    trajectory_refine!(psi, evolution, truncation, selector; kwargs...)

Scan a short projected trajectory and keep the state with the best selector score.
"""
function trajectory_refine!(psi, evolution, truncation, selector; kwargs...)
  isnothing(selector) && return psi

  refine_steps = get(kwargs, :refine_steps, 10)
  target_energy = get(kwargs, :target_energy, nothing)
  selector_context = get(kwargs, :selector_context, SelectionContext())
  candidate = _clone_state(psi)
  best_state = _clone_state(psi)
  best_score = score(selector, best_state, selector_context)

  for _ in 1:refine_steps
    scarfinder_step!(candidate, evolution, truncation; target_energy=target_energy, selector=selector)
    candidate_score = score(selector, candidate, selector_context)
    if candidate_score < best_score
      best_state = _clone_state(candidate)
      best_score = candidate_score
    end
  end

  _assign_state!(psi, best_state)
  return psi
end

"""
    scarfinder_step!(psi, evolution, truncation; target_energy=nothing, selector=nothing, kwargs...)

Run one ScarFinder step consisting of evolution, projection, and optional energy correction.
"""
function scarfinder_step!(psi, evolution, truncation; target_energy=nothing, selector=nothing, kwargs...)
  evolve!(psi, evolution)
  project!(psi, truncation)
  isnothing(target_energy) || match_energy!(psi, evolution, truncation, target_energy)
  return psi
end

"""
    scarfinder!(psi, evolution, truncation, niter; refine=false, selector=nothing, kwargs...)

Run multiple ScarFinder steps and optionally refine the final state by selector-based trajectory search.
"""
function scarfinder!(psi, evolution, truncation, niter; refine=false, selector=nothing, kwargs...)
  target_energy = get(kwargs, :target_energy, nothing)
  for _ in 1:niter
    scarfinder_step!(psi, evolution, truncation; target_energy=target_energy, selector=selector, kwargs...)
  end
  if refine
    trajectory_refine!(psi, evolution, truncation, selector; kwargs...)
  end
  return psi
end
