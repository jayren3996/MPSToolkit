"""
    evolve!(psi, evolution)

Backend-dispatched in-place evolution entry point used throughout `MPSToolkit`.

# Arguments
- `psi`: Mutable state to evolve. In this package that is typically an `MPS`, but
  custom backends used in tests or downstream extensions may define their own state types.
- `evolution`: Evolution configuration object describing how one logical evolution call
  should be carried out.

# Returns
- The same `psi` object after in-place mutation.

# Notes
- This fallback method exists only to provide a common dispatch point.
- Concrete methods are defined for `LocalGateEvolution`, `DMTGateEvolution`, and
  `TDVPEvolution`.
"""
function evolve!(psi, evolution)
  throw(MethodError(evolve!, (psi, evolution)))
end
