"""
    DMTCheckpoint

One snapshot of a long operator-space run: the state, where it had reached, the observables
accumulated so far, and the parameters identifying the run that produced it.

# Notes
- `parameters` is what makes a resume safe rather than merely convenient. Restarting from a
  snapshot of a *different* run silently splices two time series into one plausible-looking
  curve, which no later analysis can detect; [`dmt_checkpoint_resume`](@ref) refuses unless the
  parameters match exactly.
- Written with the `Serialization` stdlib rather than HDF5. Checkpoints are short-lived, are
  read back by the same Julia on the same machine, and exist to survive a wall-clock limit —
  none of which needs a portable container, and this package deliberately carries four
  dependencies.
- Consequently a checkpoint is **not** a durable archive: `Serialization`'s format is not
  guaranteed across Julia versions. Treat it as a resume token, and write anything you want to
  keep to the run's own output file.
"""
struct DMTCheckpoint
  version::Int
  parameters::Dict{Symbol,Any}
  time::Float64
  sweep::Int
  state::MPS
  observables::Dict{Symbol,Any}
end

const _DMT_CHECKPOINT_VERSION = 1

"""
    dmt_checkpoint_path(name, time; dir=".")

Return the conventional checkpoint filename for run `name` at simulation time `time`.

# Notes
- The time is formatted to three decimals with a fixed width so the names sort chronologically
  as plain strings, which is what lets [`dmt_checkpoint_latest`](@ref) pick the newest without
  parsing every file it finds.
"""
function dmt_checkpoint_path(name::AbstractString, time::Real; dir::AbstractString=".")
  return joinpath(dir, @sprintf("%s_t_%09.3f.dmt", name, Float64(time)))
end

"""
    dmt_checkpoint_save(path, rho, parameters; time, sweep, observables=Dict{Symbol,Any}())

Write a [`DMTCheckpoint`](@ref) to `path`, atomically.

# Arguments
- `path`: Destination file, normally from [`dmt_checkpoint_path`](@ref).
- `rho`: Operator-space `MPS` to snapshot.
- `parameters`: Everything that identifies this run (chain length, local dimension, model
  couplings, `dt`, `maxdim`, ...). A resume compares these for equality.

# Keyword Arguments
- `time`: Simulation time reached.
- `sweep`: Sweep index reached.
- `observables`: Accumulated time series to carry across the restart.

# Returns
- `path`.

# Notes
- Written to a temporary file and then `mv`d into place. A checkpoint exists precisely because
  the process may be killed at an arbitrary moment, and the moment most likely to catch it is
  the one where it is writing megabytes of tensors; a half-written file that *looks* like a
  checkpoint is worse than no checkpoint.
"""
function dmt_checkpoint_save(path::AbstractString, rho::MPS, parameters::AbstractDict;
                             time::Real, sweep::Integer,
                             observables::AbstractDict=Dict{Symbol,Any}())
  checkpoint = DMTCheckpoint(_DMT_CHECKPOINT_VERSION, Dict{Symbol,Any}(parameters),
                             Float64(time), Int(sweep), rho, Dict{Symbol,Any}(observables))
  directory = dirname(path)
  isempty(directory) || mkpath(directory)
  scratch = string(path, ".partial")
  open(scratch, "w") do io
    serialize(io, checkpoint)
  end
  mv(scratch, path; force=true)
  return path
end

"""
    dmt_checkpoint_load(path)

Read back a [`DMTCheckpoint`](@ref).

# Notes
- Rejects a checkpoint written by a newer format version rather than deserializing it into the
  current field layout, which would misinterpret the fields silently.
- The site indices come back inside the state's own tensors. A resumed run must take them from
  the loaded state (`siteind(checkpoint.state, n)`) rather than calling `operator_siteinds`
  again: fresh indices carry fresh ids and will not contract with the restored tensors.
"""
function dmt_checkpoint_load(path::AbstractString)
  checkpoint = open(deserialize, path)
  checkpoint isa DMTCheckpoint ||
    throw(ArgumentError("$(path) does not contain a DMT checkpoint"))
  checkpoint.version <= _DMT_CHECKPOINT_VERSION || throw(ArgumentError(
    "checkpoint $(path) has format version $(checkpoint.version) but this MPSToolkit " *
    "understands at most $(_DMT_CHECKPOINT_VERSION)"))
  return checkpoint
end

"""
    dmt_checkpoint_latest(name; dir=".")

Return the path of the newest checkpoint for run `name`, or `nothing` when there is none.
"""
function dmt_checkpoint_latest(name::AbstractString; dir::AbstractString=".")
  isdir(dir) || return nothing
  prefix = string(name, "_t_")
  candidates = filter(f -> startswith(f, prefix) && endswith(f, ".dmt"), readdir(dir))
  isempty(candidates) && return nothing
  return joinpath(dir, maximum(candidates))          # fixed-width times sort chronologically
end

"""
    dmt_checkpoint_resume(name, parameters; dir=".")

Return the newest checkpoint for run `name`, or `nothing` to start from scratch.

Throws if a checkpoint exists but was produced by a run with different `parameters`.

# Notes
- Refusing a mismatch is the point of this function. Silently resuming someone else's run — or
  your own run at a different `maxdim` — produces a continuous-looking time series with a
  discontinuity in the physics at the join, and nothing downstream can tell. Failing loudly
  costs a rerun; succeeding quietly costs a wrong result.
- The comparison is over the whole parameter set, so adding a parameter invalidates old
  checkpoints. That is deliberate: a parameter worth recording is a parameter worth matching.
"""
function dmt_checkpoint_resume(name::AbstractString, parameters::AbstractDict;
                               dir::AbstractString=".")
  path = dmt_checkpoint_latest(name; dir=dir)
  isnothing(path) && return nothing
  checkpoint = dmt_checkpoint_load(path)
  wanted = Dict{Symbol,Any}(parameters)
  if checkpoint.parameters != wanted
    differing = [key for key in union(keys(wanted), keys(checkpoint.parameters))
                 if get(wanted, key, nothing) != get(checkpoint.parameters, key, nothing)]
    throw(ArgumentError(
      "checkpoint $(path) was written by a different run: " *
      join(["$(key) = $(get(checkpoint.parameters, key, "<absent>")) " *
            "(now $(get(wanted, key, "<absent>")))" for key in sort(differing; by=String)], ", ") *
      ". Delete it or choose a different run name; resuming across a parameter change would " *
      "splice two different physics into one time series."))
  end
  return checkpoint
end
