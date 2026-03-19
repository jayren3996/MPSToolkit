"""
    score(selector, psi, context)

Score a state using the rule encoded by `selector`.

# Arguments
- `selector`: Selector configuration, such as [`EntropySelector`](@ref) or
  [`FidelitySelector`](@ref).
- `psi`: Candidate state to score.
- `context`: Optional auxiliary information shared across selector calls.

# Returns
- A scalar score, where smaller values are considered better by ScarFinder refinement.
"""
function score(selector::EntropySelector, psi, context::SelectionContext=SelectionContext())
  return bond_entropy(psi, selector.bond)
end

"""
    score(selector::FidelitySelector, psi, context)

Score a state by its fidelity distance to `context.reference_state`.

# Arguments
- `selector`: Fidelity-based selector configuration.
- `psi`: Candidate state to compare.
- `context`: [`SelectionContext`](@ref) whose `reference_state` field must be populated.

# Returns
- `1 - |⟨reference_state|psi⟩|`, so lower values correspond to larger fidelity.
"""
function score(selector::FidelitySelector, psi, context::SelectionContext)
  return fidelity_distance(psi, context.reference_state)
end
