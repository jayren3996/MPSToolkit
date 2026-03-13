"""
    score(selector, psi, context)

Score a state using the rule encoded by `selector`.
"""
function score(selector::EntropySelector, psi, context::SelectionContext=SelectionContext())
  return bond_entropy(psi, selector.bond)
end

"""
    score(selector::FidelitySelector, psi, context)

Score a state by its fidelity distance to `context.reference_state`.
"""
function score(selector::FidelitySelector, psi, context::SelectionContext)
  return fidelity_distance(psi, context.reference_state)
end
