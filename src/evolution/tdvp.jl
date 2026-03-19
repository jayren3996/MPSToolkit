"""
    tdvp_evolve!(psi, evo)

Run one TDVP evolution call on a finite OBC `MPS`.

# Arguments
- `psi`: Finite matrix-product state to evolve in place.
- `evo`: [`TDVPEvolution`](@ref) carrying the MPO generator and TDVP solver settings.

# Returns
- The mutated `psi`.

# Notes
- The effective TDVP step count is `evo.nsteps` if present, otherwise `evo.nsweeps`.
- The underlying `tdvp` call returns a new state, so this wrapper copies the result back
  into the original storage to preserve the package-wide in-place API convention.
"""
function tdvp_evolve!(psi::MPS, evo::TDVPEvolution)
  nsteps = isnothing(evo.nsteps) ? evo.nsweeps : evo.nsteps

  kwargs = (
    time_step=evo.time_step,
    nsteps=nsteps,
    reverse_step=evo.reverse_step,
    updater_backend=evo.updater_backend,
    normalize=evo.normalize,
    evo.solver_kwargs...,
  )

  evolved = if isnothing(evo.updater)
    tdvp(evo.generator, evo.t, psi; kwargs...)
  else
    tdvp(evo.generator, evo.t, psi; updater=evo.updater, kwargs...)
  end

  psi[:] = evolved
  return psi
end

"""
    evolve!(psi, evo::TDVPEvolution)

Dispatch finite-`MPS` evolution through the TDVP backend.

# Arguments
- `psi`: State to evolve.
- `evo`: [`TDVPEvolution`](@ref) object.

# Returns
- The mutated `psi`.
"""
function evolve!(psi::MPS, evo::TDVPEvolution)
  return tdvp_evolve!(psi, evo)
end
