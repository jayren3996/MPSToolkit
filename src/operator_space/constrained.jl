"""
    constrained_dmt_evolve!(rho, evo, projector; project_every=1, projector_maxdim=evo.gate_maxdim, projector_cutoff=evo.cutoff)

Run scheduled operator-space DMT evolution with periodic constraint-projection checkpoints.

The driver executes the `evo.nstep` forward-plus-reverse DMT sweeps of
[`dmt_evolve!`](@ref) in chunks of `project_every` sweeps; after each chunk (including the
final, possibly shorter one) it applies the operator-space `projector` MPO with
`ITensorMPS.apply` and renormalizes. For a constrained model such as PXP the exact gates
commute with the constraint projector, so projection would be a no-op for exact evolution —
the checkpoints exist to remove the weight that **truncation** leaks out of the constrained
sector before it contaminates transport observables.

# Arguments
- `rho`: Operator-space `MPS` to mutate in place.
- `evo`: [`DMTGateEvolution`](@ref) describing gates, schedules, and truncation budgets.
- `projector`: Operator-space MPO enforcing the constraint, e.g.
  [`pauli_pxp_constraint_projector`](@ref).

# Keyword Arguments
- `project_every`: Number of complete DMT sweeps between checkpoints.
- `projector_maxdim`: Bond-dimension cap for the projector application. The default reuses
  `evo.gate_maxdim` deliberately: like a raw gate application, the projector may inflate the
  bond, and the **next DMT sweep** then performs the transport-aware compression instead of
  an ordinary SVD making the truncation decision at the checkpoint.
- `projector_cutoff`: Truncation cutoff for the projector application.

# Returns
- The mutated, projected, and normalized `rho`.
"""
function constrained_dmt_evolve!(
  rho::MPS,
  evo::DMTGateEvolution,
  projector::MPO;
  project_every::Integer=1,
  projector_maxdim::Integer=evo.gate_maxdim,
  projector_cutoff::Real=evo.cutoff,
)
  project_every >= 1 || throw(ArgumentError("constrained_dmt_evolve! requires project_every >= 1"))
  projector_maxdim >= 1 || throw(ArgumentError("constrained_dmt_evolve! requires projector_maxdim >= 1"))
  length(projector) == length(rho) || throw(ArgumentError("projector and state must have matching lengths"))

  remaining = evo.nstep
  while remaining > 0
    chunk = min(Int(project_every), remaining)
    chunk_evo = DMTGateEvolution(
      evo.gate,
      evo.dt;
      schedule=evo.schedule,
      reverse_schedule=evo.reverse_schedule,
      nstep=chunk,
      maxdim=evo.maxdim,
      cutoff=evo.cutoff,
      gate_maxdim=evo.gate_maxdim,
      connector_buffer=evo.connector_buffer,
    )
    dmt_evolve!(rho, chunk_evo)
    rho[:] = apply(projector, rho; maxdim=Int(projector_maxdim), cutoff=projector_cutoff)
    normalize!(rho)
    remaining -= chunk
  end
  return rho
end
