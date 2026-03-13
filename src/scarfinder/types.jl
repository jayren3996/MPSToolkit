"""
    BondDimTruncation

Truncation settings used by ScarFinder projection routines.
"""
struct BondDimTruncation
  maxdim::Int
  cutoff::Float64
end

"""
    BondDimTruncation(maxdim; cutoff=0.0)

Construct bond-dimension truncation settings.
"""
function BondDimTruncation(maxdim; cutoff=0.0)
  return BondDimTruncation(Int(maxdim), Float64(cutoff))
end

"""
    EnergyTarget

Parameters controlling post-evolution energy correction toward a target expectation value.
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

Construct an energy-targeting configuration for `match_energy!`.
"""
function EnergyTarget(target; operator=nothing, tol=1e-6, alpha=0.1, maxstep=50)
  return EnergyTarget(Float64(target), operator, Float64(tol), Float64(alpha), Int(maxstep))
end

"""
    SelectionContext

Additional reference data passed to selector scoring during ScarFinder refinement.
"""
struct SelectionContext
  reference_state
end

"""
    SelectionContext(; reference_state=nothing)

Construct selector context carrying an optional reference state or other external data.
"""
function SelectionContext(; reference_state=nothing)
  return SelectionContext(reference_state)
end

"""
    EntropySelector

Selector configuration that scores states using bond entanglement entropy.
"""
struct EntropySelector
  bond::Union{Nothing,Int}
end

"""
    EntropySelector(; bond=nothing)

Construct an entropy selector for a specific bond or a backend-default bond if `bond` is `nothing`.
"""
function EntropySelector(; bond=nothing)
  return EntropySelector(bond)
end

"""
    FidelitySelector

Selector configuration that scores states using fidelity against a reference state.
"""
struct FidelitySelector
end
