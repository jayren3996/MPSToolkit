"""
    constrained_dmt_evolve!(rho, evo, projector; project_every=1, projector_maxdim=2*evo.maxdim, projector_cutoff=evo.cutoff, normalize=evo.normalize)

Run scheduled operator-space DMT evolution with periodic constraint-projection checkpoints.

!!! warning "Density operators only"
    Like [`dmt_evolve!`](@ref), this is a DMT driver: `rho` must be a near-infinite-temperature
    **density operator** (e.g. a constrained energy domain-wall melt). DMT protects the trace
    component, so it must **not** be used to Heisenberg-evolve a **traceless** operator / two-point
    correlator — use ordinary TEBD for those.

The driver executes the `evo.nstep` forward-plus-reverse DMT sweeps of
[`dmt_evolve!`](@ref) in chunks of `project_every` sweeps; after each chunk (including the
final, possibly shorter one) it applies the operator-space `projector` MPO with
`ITensorMPS.apply`. For a constrained model such as PXP the exact gates commute with the
constraint projector, so projection would be a no-op for exact evolution — the checkpoints
exist to remove the weight that **truncation** leaks out of the constrained sector before it
contaminates transport observables.

# Arguments
- `rho`: Operator-space `MPS` to mutate in place.
- `evo`: [`DMTGateEvolution`](@ref) describing gates, schedules, and truncation budgets.
- `projector`: Operator-space MPO enforcing the constraint, e.g.
  [`pauli_pxp_constraint_projector`](@ref).

# Keyword Arguments
- `project_every`: Number of complete DMT sweeps between checkpoints.
- `projector_maxdim`: Bond-dimension cap for the projector application, default
  `2 * evo.maxdim`. The checkpoint state is only an `O(leakage)` perturbation of the
  pre-projection state (which DMT just compressed to `evo.maxdim`), so compressing it back
  immediately is benign — the kept subspace rotates by `O(leakage)`, it does not lose the
  DMT-protected components wholesale. The modest `2x` buffer gives the cutoff room to act.
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
- `normalize`: Whether to renormalize after each checkpoint (and let the chunked
  `dmt_evolve!` sweeps renormalize). Defaults to `evo.normalize`, so a
  `DMTGateEvolution(...; normalize=false)` is honored without re-passing the keyword here —
  matching [`dmt_evolve!`](@ref) and [`evolve!`](@ref). Set `false` to preserve absolute
  operator scales, which is required when tracking unnormalized traces of a **traceless**
  operator (e.g. conserved `tr(H O(t))` in correlator-protocol transport runs).

# Returns
- The mutated and projected (and, by default, normalized) `rho`.
"""
function constrained_dmt_evolve!(
  rho::MPS,
  evo::DMTGateEvolution,
  projector::MPO;
  project_every::Integer=1,
  projector_maxdim::Integer=2 * evo.maxdim,
  projector_cutoff::Real=evo.cutoff,
  normalize::Bool=evo.normalize,
)
  project_every >= 1 || throw(ArgumentError("constrained_dmt_evolve! requires project_every >= 1"))
  projector_maxdim >= 1 || throw(ArgumentError("constrained_dmt_evolve! requires projector_maxdim >= 1"))
  projector_cutoff >= 0 || throw(ArgumentError("constrained_dmt_evolve! requires projector_cutoff >= 0"))
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
      preserve_diameter=evo.preserve_diameter,
      truncation=evo.truncation,
    )
    dmt_evolve!(rho, chunk_evo; normalize=normalize)
    rho[:] = apply(projector, rho; maxdim=Int(projector_maxdim), cutoff=projector_cutoff)
    normalize && normalize!(rho)
    remaining -= chunk
  end
  return rho
end
