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
  nstep >= 1 || throw(ArgumentError("LocalGateEvolution requires nstep >= 1"))
  maxdim >= 0 || throw(ArgumentError("LocalGateEvolution requires maxdim >= 0"))
  cutoff >= 0 || throw(ArgumentError("LocalGateEvolution requires cutoff >= 0"))
  return LocalGateEvolution(gate, Float64(dt), schedule, Int(nstep), Int(maxdim), Float64(cutoff))
end

"""
    DMTGateEvolution

Configuration for scheduled operator-space DMT evolution.

!!! warning "Transport-specific algorithm"
    DMT is designed for **transport problems** (e.g. spin or energy diffusion) in operator
    space.  Its truncation rule protects local reduced operator data, including the
    identity/trace component and nearby local-operator components, before truncating connected
    long-range correlations. For general-purpose operator-space TEBD, use
    [`LocalGateEvolution`](@ref).

# Fields
- `gate`: Dense local gate specification in the operator basis.
- `dt`: Logical time increment associated with one full DMT evolution call.
- `schedule`: Forward update schedule.
- `reverse_schedule`: Reverse update schedule used for the backward sweep.
- `nstep`: Number of complete forward-plus-reverse sweeps per `evolve!` call.
- `maxdim`: **Total** post-DMT bond dimension, inclusive of the protected block.
- `cutoff`: Truncation cutoff used in the final refactorization.
- `gate_maxdim`: Temporary bond dimension cap for raw gate application; `0` means no cap, i.e.
  the gate is applied exactly.
- `preserve_diameter`: Positive odd diameter of the observables DMT preserves exactly.
- `truncation`: `:dense` (default) or `:random` complement truncation; `:random` is faster at
  large bond dimension but not deterministic. See [`DMTOptions`](@ref).
- `normalize`: Whether `evolve!` / `dmt_evolve!` renormalize the state after evolution.
  Default `true`; set `false` to track unnormalized traces of a traceless operator.
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
  preserve_diameter::Int
  truncation::Symbol
  normalize::Bool
end

"""
    _reject_connector_buffer(connector_buffer)

Throw an `ArgumentError` describing the DMT budget migration unless `connector_buffer` is
`nothing`.

# Notes
- Shared by `DMTOptions`, [`DMTGateEvolution`](@ref) and `dmt_step!`, which all used to accept
  the keyword. A bare `MethodError` would not tell a caller that `maxdim` also changed meaning.
"""
function _reject_connector_buffer(connector_buffer)
  isnothing(connector_buffer) || throw(ArgumentError(
    "connector_buffer was removed: DMT now protects the d^(preserve_diameter - 1) local-operator " *
    "subspace on each side structurally. Use preserve_diameter (odd, default 3) instead, and note " *
    "that maxdim is now the total bond dimension including the protected block " *
    "(floor 2 d^(preserve_diameter - 1) + 1)."))
  return nothing
end

"""
    DMTGateEvolution(gate, dt; schedule, reverse_schedule=reverse(schedule), nstep=1, maxdim=30, cutoff=1e-12, gate_maxdim=0, preserve_diameter=3, truncation=:dense, normalize=true)

Construct a [`DMTGateEvolution`](@ref) for **transport** simulations.

# Arguments
- `gate`: Dense local gate specification in the operator basis.
- `dt`: Logical time increment associated with one `dmt_evolve!` call.

# Keyword Arguments
- `schedule`: Forward update schedule.
- `reverse_schedule`: Reverse update schedule. By default the forward schedule is reversed.
- `nstep`: Number of complete forward-plus-reverse sweeps per evolution call.
- `maxdim`: **Total** bond dimension after DMT truncation, inclusive of the protected block; it
  must be at least `2 d^(preserve_diameter - 1) + 1` for the local dimension `d` in use.
- `cutoff`: Truncation cutoff used when refactorizing the compressed bond.
- `gate_maxdim`: Temporary gate-application bond dimension cap. **`0` (the default) means no
  cap: the gate is applied exactly.** A positive cap pre-truncates the inflated bond with a
  plain SVD, discarding the smallest singular values *before* DMT can protect the local-operator
  content they carry, which is precisely the error DMT exists to avoid. The previous default,
  `max(maxdim * 16, 64)`, capped nothing for any `d <= 4` in steady state (a two-site gate
  inflates the bond from the incoming `chi` to `d^2 * chi`, and `chi <= maxdim` once DMT has
  truncated that bond), so this is a no-op there and a correctness fix at `d >= 5`.
- `preserve_diameter`: Positive odd diameter of the observables preserved exactly;
  `radius = (preserve_diameter - 1) / 2` sites are protected on each side of the cut.
- `truncation`: `:dense` (default) or `:random` complement truncation. `:random` measures
  1.05x-1.2x faster on a whole sweep at moderate budgets and ~1.4x once the gate-inflated bond
  passes ~2500, and preserves the guarantee to the same `1e-15` — but it draws from the global
  RNG, so it is **not deterministic**: two truncations of the same bond agree only to
  randomized-SVD accuracy. See [`DMTOptions`](@ref).
- `normalize`: Default normalization choice carried by the object; `evolve!` / `dmt_evolve!`
  use it unless overridden by their own `normalize` keyword. Set `false` for traceless
  operators (see [`dmt_evolve!`](@ref)).

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
  gate_maxdim=0,
  preserve_diameter=3,
  truncation=:dense,
  normalize=true,
  connector_buffer=nothing,
)
  _reject_connector_buffer(connector_buffer)
  nstep >= 1 || throw(ArgumentError("DMTGateEvolution requires nstep >= 1"))
  maxdim >= 1 || throw(ArgumentError("DMTGateEvolution requires maxdim >= 1"))
  cutoff >= 0 || throw(ArgumentError("DMTGateEvolution requires cutoff >= 0"))
  gate_maxdim >= 0 ||
    throw(ArgumentError("DMTGateEvolution requires gate_maxdim >= 0 (0 = no cap)"))
  isodd(preserve_diameter) && preserve_diameter >= 1 || throw(ArgumentError(
    "DMTGateEvolution requires a positive odd preserve_diameter, got $(preserve_diameter)"))
  truncation in (:dense, :random) || throw(ArgumentError(
    "DMTGateEvolution truncation must be :dense or :random, got $(truncation)"))
  return DMTGateEvolution(
    gate,
    Float64(dt),
    schedule,
    reverse_schedule,
    Int(nstep),
    Int(maxdim),
    Float64(cutoff),
    Int(gate_maxdim),
    Int(preserve_diameter),
    Symbol(truncation),
    Bool(normalize),
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
struct TDVPEvolution{TH,TT,TTS,TS,TB,TU,TK}
  generator::TH
  t::TT
  time_step::TTS
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
- `solver_kwargs`: Additional keyword arguments forwarded to `tdvp` (e.g. `maxdim`, `cutoff`).
  Must not contain a key that already has a dedicated argument above (`time_step`, `nsteps`,
  `reverse_step`, `updater_backend`, `updater`, `normalize`); passing one throws, because it
  would otherwise be forwarded last and silently override the dedicated argument.
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
  isnothing(nsteps) || nsteps >= 1 || throw(ArgumentError("TDVPEvolution requires nsteps >= 1 when provided"))
  isnothing(nsweeps) || nsweeps >= 1 || throw(ArgumentError("TDVPEvolution requires nsweeps >= 1 when provided"))
  reserved_solver_keys = (:time_step, :nsteps, :reverse_step, :updater_backend, :updater, :normalize)
  conflicting_keys = intersect(keys(solver_kwargs), reserved_solver_keys)
  isempty(conflicting_keys) ||
    throw(ArgumentError("TDVPEvolution solver_kwargs must not contain reserved keys $(collect(conflicting_keys)); these are forwarded to `tdvp` last and would silently override the dedicated TDVPEvolution keyword arguments. Set them through the dedicated keyword arguments instead."))
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
