"""
    BondDimTruncation

Truncation settings used by ScarFinder projection routines.

# Fields
- `maxdim`: Maximum allowed bond dimension after projection.
- `cutoff`: Singular-value cutoff used during truncation.
"""
struct BondDimTruncation
  maxdim::Int
  cutoff::Float64
end

"""
    BondDimTruncation(maxdim; cutoff=0.0)

Construct [`BondDimTruncation`](@ref) settings.

# Arguments
- `maxdim`: Maximum bond dimension allowed after projection.

# Keyword Arguments
- `cutoff`: Singular-value cutoff used by the projection backend.

# Returns
- A `BondDimTruncation` object with normalized numeric field types.
"""
function BondDimTruncation(maxdim; cutoff=0.0)
  return BondDimTruncation(Int(maxdim), Float64(cutoff))
end

"""
    EnergyTarget

Parameters controlling post-evolution energy correction toward a target expectation value.

# Fields
- `target`: Desired energy density or expectation value.
- `operator`: Dense operator or `MPO` used to measure and correct the energy.
- `tol`: Early-stop tolerance on the residual energy mismatch.
- `alpha`: Proportional step-size parameter for the correction update.
- `maxstep`: Maximum number of correction iterations per ScarFinder step.
"""
struct EnergyTarget
  target::Float64
  operator
  tol::Float64
  alpha::Float64
  maxstep::Int
end

"""
    EnergyTarget(target; operator=nothing, tol=1e-6, alpha=0.1, maxstep=50)

Construct an [`EnergyTarget`](@ref) configuration for [`match_energy!`](@ref).

# Arguments
- `target`: Desired post-step energy value.

# Keyword Arguments
- `operator`: Dense operator or `MPO` used to measure the current energy.
- `tol`: Early-stop tolerance.
- `alpha`: Proportional update strength.
- `maxstep`: Maximum number of correction iterations.

# Returns
- An `EnergyTarget` object.
"""
function EnergyTarget(target; operator=nothing, tol=1e-6, alpha=0.1, maxstep=50)
  return EnergyTarget(Float64(target), operator, Float64(tol), Float64(alpha), Int(maxstep))
end

"""
    SelectionContext

Additional reference data passed to selector scoring during ScarFinder refinement.

# Fields
- `reference_state`: Optional state used by selectors such as [`FidelitySelector`](@ref).
"""
struct SelectionContext
  reference_state
end

"""
    SelectionContext(; reference_state=nothing)

Construct a [`SelectionContext`](@ref).

# Keyword Arguments
- `reference_state`: Optional reference state shared across selector calls.

# Returns
- A `SelectionContext` object.
"""
function SelectionContext(; reference_state=nothing)
  return SelectionContext(reference_state)
end

"""
    EntropySelector

Selector configuration that scores states using bond entanglement entropy.

# Fields
- `bond`: Bond whose entropy should be minimized. `nothing` means "use the backend default".
"""
struct EntropySelector
  bond::Union{Nothing,Int}
end

"""
    EntropySelector(; bond=nothing)

Construct an [`EntropySelector`](@ref).

# Keyword Arguments
- `bond`: Bond index to score. If `nothing`, the backend default bond is used.

# Returns
- An `EntropySelector` object.
"""
function EntropySelector(; bond=nothing)
  return EntropySelector(bond)
end

"""
    FidelitySelector

Selector configuration that scores states using fidelity against a reference state.

# Notes
- The required reference state is supplied separately through [`SelectionContext`](@ref).
"""
struct FidelitySelector
end
