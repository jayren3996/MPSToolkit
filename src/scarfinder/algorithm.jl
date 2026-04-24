"""
    _assign_state!(psi, updated)

Overwrite the storage of `psi` with the contents of `updated`.

# Arguments
- `psi`: Mutable destination state.
- `updated`: Source state whose tensor data should replace `psi`.

# Returns
- The mutated `psi`.

# Notes
- Keeping this helper separate makes it easy for tests and custom state wrappers to
  override how assignment happens.
"""
function _assign_state!(psi::MPS, updated::MPS)
  psi[:] = updated
  return psi
end

"""
    _scarfinder_default_steps()

Return the effective evolution step count that ScarFinder should use when a caller
provides a single-step evolution object.

# Returns
- The integer step count that ScarFinder substitutes for a main-loop evolution configured
  with exactly one step.
"""
function _scarfinder_default_steps()
  return 10
end

"""
    _scarfinder_effective_steps(evolution)

Return the effective per-call evolution step count that ScarFinder will observe for
`evolution`, or `nothing` when the object does not expose a known step-count field.

# Arguments
- `evolution`: Evolution object that may or may not expose `nstep`-style metadata.

# Returns
- An integer step count, or `nothing` if no supported field is present.
"""
function _scarfinder_effective_steps(evolution)
  return nothing
end

function _scarfinder_effective_steps(evolution::LocalGateEvolution)
  return evolution.nstep
end

function _scarfinder_effective_steps(evolution::DMTGateEvolution)
  return evolution.nstep
end

function _scarfinder_effective_steps(evolution::TDVPEvolution)
  return isnothing(evolution.nsteps) ? evolution.nsweeps : evolution.nsteps
end

"""
    _scarfinder_rebuild_evolution(evolution, steps)

Return an evolution object equivalent to `evolution` except for the effective step count,
which is replaced with `steps`.

# Arguments
- `evolution`: Existing evolution configuration.
- `steps`: Replacement effective step count.

# Returns
- A rebuilt evolution object, or `evolution` unchanged for unsupported types.
"""
function _scarfinder_rebuild_evolution(evolution, steps)
  return evolution
end

function _scarfinder_rebuild_evolution(evolution::LocalGateEvolution, steps)
  return LocalGateEvolution(
    evolution.gate,
    evolution.dt;
    schedule=evolution.schedule,
    nstep=steps,
    maxdim=evolution.maxdim,
    cutoff=evolution.cutoff,
  )
end

function _scarfinder_rebuild_evolution(evolution::DMTGateEvolution, steps)
  return DMTGateEvolution(
    evolution.gate,
    evolution.dt;
    schedule=evolution.schedule,
    reverse_schedule=evolution.reverse_schedule,
    nstep=steps,
    maxdim=evolution.maxdim,
    cutoff=evolution.cutoff,
    gate_maxdim=evolution.gate_maxdim,
    connector_buffer=evolution.connector_buffer,
  )
end

function _scarfinder_rebuild_evolution(evolution::TDVPEvolution, steps)
  return TDVPEvolution(
    evolution.generator,
    evolution.t;
    time_step=evolution.time_step,
    nsteps=steps,
    nsweeps=evolution.nsweeps,
    reverse_step=evolution.reverse_step,
    updater_backend=evolution.updater_backend,
    updater=evolution.updater,
    normalize=evolution.normalize,
    solver_kwargs=evolution.solver_kwargs,
    schedule=evolution.schedule,
  )
end

"""
    _scarfinder_evolution(evolution; warn=true)

Normalize an evolution object for ScarFinder use. If the effective step count is `1`,
ScarFinder treats that as an unsupported main-loop setting, emits a warning, and uses
`10` steps instead for this ScarFinder call only.

# Arguments
- `evolution`: Evolution configuration passed into a ScarFinder entry point.

# Keyword Arguments
- `warn`: Whether to emit the ScarFinder step-normalization warning when a replacement is
  made.

# Returns
- Either the original evolution object or a rebuilt one with step count `10`.

# Notes
- This helper intentionally does not touch the internal one-step correction evolutions used
  by `match_energy!`. Those are narrow post-step correction loops rather than the main
  ScarFinder trajectory evolution.
"""
function _scarfinder_evolution(evolution; warn::Bool=true)
  effective_steps = _scarfinder_effective_steps(evolution)
  effective_steps == 1 || return evolution

  replacement_steps = _scarfinder_default_steps()
  if warn
    @warn "ScarFinder received an evolution with effective step count 1 and is using $(replacement_steps) instead. Pass nstep=$(replacement_steps) or nsteps=$(replacement_steps) explicitly to silence this warning."
  end
  return _scarfinder_rebuild_evolution(evolution, replacement_steps)
end

"""
    _clone_state(psi)

Create a working copy of a state for refinement or rollback operations.

# Arguments
- `psi`: State to duplicate.

# Returns
- A deep copy of `psi`.
"""
function _clone_state(psi)
  return deepcopy(psi)
end

"""
    _match_energy_dense!(psi, evolution, truncation, target)

Apply a dense-local-operator post-evolution correction loop.

# Arguments
- `psi`: State to mutate in place.
- `evolution`: Main ScarFinder evolution object whose schedule is reused for correction.
- `truncation`: Projection settings applied after each correction substep.
- `target`: [`EnergyTarget`](@ref) defining the desired energy and correction hyperparameters.

# Returns
- The mutated `psi`.

# Notes
- The correction gate is always applied with `nstep=1` because it represents a single local
  rollback-or-improvement move, not the main ScarFinder trajectory evolution.
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

Map an energy error to the signed imaginary-time interval used by an MPO correction step.

# Arguments
- `target`: [`EnergyTarget`](@ref) object.
- `energy_error`: Current signed deviation from the target energy.

# Returns
- A purely imaginary correction interval proportional to the error.
"""
function _correction_time(target, energy_error)
  return -1im * target.alpha * energy_error
end

"""
    _solver_kwargs(evolution, truncation)

Build TDVP solver keyword settings for an MPO-based energy-correction step.

# Arguments
- `evolution`: Reference evolution object.
- `truncation`: Projection settings that may provide fallback truncation metadata.

# Returns
- A named tuple of keyword arguments suitable for `TDVPEvolution(...; solver_kwargs=...)`.
"""
function _solver_kwargs(evolution::TDVPEvolution, truncation)
  return evolution.solver_kwargs
end

function _solver_kwargs(evolution::TDVPEvolution, truncation::BondDimTruncation)
  return evolution.solver_kwargs
end

"""
    _solver_kwargs(evolution, truncation)

Build TDVP solver keyword settings from a non-TDVP evolution object.

# Notes
- This fallback extracts `maxdim` and `cutoff` from either the evolution object or the
  truncation settings so that an MPO-based correction step still has sensible solver
  parameters.
"""
function _solver_kwargs(evolution, truncation::BondDimTruncation)
  maxdim = hasproperty(evolution, :maxdim) ? getproperty(evolution, :maxdim) : truncation.maxdim
  cutoff = hasproperty(evolution, :cutoff) ? getproperty(evolution, :cutoff) : truncation.cutoff
  return (; maxdim=maxdim, cutoff=cutoff)
end

"""
    _tdvp_correction_evolution(operator, correction_time, evolution, truncation)

Construct the TDVP evolution object used for one MPO-based energy-correction move.

# Arguments
- `operator`: Correction `MPO`.
- `correction_time`: Imaginary-time interval produced by [`_correction_time`](@ref).
- `evolution`: Reference evolution object from the main ScarFinder loop.
- `truncation`: Projection settings used to infer solver defaults when needed.

# Returns
- A `TDVPEvolution` with `nsteps=1`.
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

Apply an MPO-based TDVP correction loop that nudges a finite `MPS` toward a target energy.

# Arguments
- `psi`: State to mutate in place.
- `evolution`: Reference evolution object from the main ScarFinder loop.
- `truncation`: Projection settings applied after each correction substep.
- `target`: [`EnergyTarget`](@ref) controlling the correction loop.

# Returns
- The mutated `psi`.
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

# Arguments
- `psi`: State to mutate in place.
- `evolution`: Main ScarFinder evolution object.
- `truncation`: Projection settings reused after each correction substep.
- `target`: [`EnergyTarget`](@ref) configuration.

# Returns
- The mutated `psi`.

# Notes
- Dense operators use the local-gate correction path; `MPO`s use a one-step TDVP correction
  path.
"""
function match_energy!(psi, evolution, truncation, target)
  isnothing(target.operator) && return psi
  target.operator isa MPO && return _match_energy_mpo!(psi, evolution, truncation, target)
  return _match_energy_dense!(psi, evolution, truncation, target)
end

"""
    trajectory_refine!(psi, evolution, truncation, selector; kwargs...)

Scan a short projected trajectory and keep the state with the best selector score.

# Arguments
- `psi`: Initial state and eventual output state.
- `evolution`: Evolution object used to advance the candidate trajectory.
- `truncation`: Projection settings applied after each evolution step.
- `selector`: Selector configuration used to rank candidate states.

# Keyword Arguments
- `refine_steps`: Number of extra projected steps to scan.
- `target_energy`: Optional [`EnergyTarget`](@ref) reused during refinement.
- `selector_context`: Optional [`SelectionContext`](@ref) shared across scoring calls.

# Returns
- The mutated `psi`, overwritten with the best state found during refinement.

# Notes
- Public `trajectory_refine!` first normalizes the evolution object through the internal
  `_scarfinder_evolution` helper, so the ScarFinder single-step warning rule also applies
  here when the function is called directly.
"""
function _trajectory_refine!(psi, evolution, truncation, selector; kwargs...)
  isnothing(selector) && return psi

  refine_steps = get(kwargs, :refine_steps, 10)
  target_energy = get(kwargs, :target_energy, nothing)
  selector_context = get(kwargs, :selector_context, SelectionContext())
  candidate = _clone_state(psi)
  best_state = _clone_state(psi)
  best_score = score(selector, best_state, selector_context)

  for _ in 1:refine_steps
    _scarfinder_step!(candidate, evolution, truncation; target_energy=target_energy)
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

Run one bare ScarFinder step consisting of evolution, projection, and optional energy correction.

# Arguments
- `psi`: State to mutate in place.
- `evolution`: Normalized evolution object to execute.
- `truncation`: Projection settings applied immediately after evolution.

# Keyword Arguments
- `target_energy`: Optional [`EnergyTarget`](@ref) applied after projection.

# Returns
- The mutated `psi`.
"""
function _scarfinder_step!(psi, evolution, truncation; target_energy=nothing)
  evolve!(psi, evolution)
  project!(psi, truncation)
  isnothing(target_energy) || match_energy!(psi, evolution, truncation, target_energy)
  return psi
end

"""
    trajectory_refine!(psi, evolution, truncation, selector; kwargs...)

Public ScarFinder refinement entry point.

# Notes
- This wrapper applies ScarFinder step-count normalization before dispatching to the
  internal refinement routine.
"""
function trajectory_refine!(psi, evolution, truncation, selector; kwargs...)
  return _trajectory_refine!(psi, _scarfinder_evolution(evolution), truncation, selector; kwargs...)
end

"""
    scarfinder_step!(psi, evolution, truncation; target_energy=nothing, selector=nothing, kwargs...)

Run one public ScarFinder step.

# Arguments
- `psi`: State to mutate in place.
- `evolution`: Evolution object for the main ScarFinder update.
- `truncation`: Projection settings applied after evolution.

# Keyword Arguments
- `target_energy`: Optional [`EnergyTarget`](@ref) used for post-step correction.
- `selector`: Accepted for interface symmetry with [`scarfinder!`](@ref). It is not used by
  a single-step update.
- `kwargs...`: Reserved for forward compatibility with higher-level workflow wrappers.

# Returns
- The mutated `psi`.

# Notes
- If the effective step count of `evolution` is `1`, ScarFinder emits a warning and uses
  `10` internally for this call.
"""
function scarfinder_step!(psi, evolution, truncation; target_energy=nothing, selector=nothing, kwargs...)
  return _scarfinder_step!(psi, _scarfinder_evolution(evolution), truncation; target_energy=target_energy)
end

"""
    scarfinder!(psi, evolution, truncation, niter; refine=false, selector=nothing, kwargs...)

Run a multi-step ScarFinder loop.

# Arguments
- `psi`: State to mutate in place.
- `evolution`: Evolution object for the main ScarFinder updates.
- `truncation`: Projection settings applied after each step.
- `niter`: Number of projected evolution steps to perform.

# Keyword Arguments
- `refine`: If `true`, run selector-based trajectory refinement after the main loop.
- `selector`: Selector configuration used during refinement.
- `target_energy`: Optional [`EnergyTarget`](@ref) reused on every step.
- `kwargs...`: Additional refinement keywords forwarded to [`trajectory_refine!`](@ref),
  such as `refine_steps` or `selector_context`.

# Returns
- The mutated `psi`.

# Notes
- The ScarFinder single-step normalization rule is applied once up front, so repeated loop
  iterations do not emit repeated warnings.
"""
function scarfinder!(psi, evolution, truncation, niter; refine=false, selector=nothing, kwargs...)
  scar_evolution = _scarfinder_evolution(evolution)
  target_energy = get(kwargs, :target_energy, nothing)
  for _ in 1:niter
    _scarfinder_step!(psi, scar_evolution, truncation; target_energy=target_energy)
  end
  if refine
    trajectory_refine!(psi, scar_evolution, truncation, selector; kwargs...)
  end
  return psi
end
