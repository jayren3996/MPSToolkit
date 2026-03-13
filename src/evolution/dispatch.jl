"""
    evolve!(psi, evolution)

Backend-dispatched in-place evolution entry point used by MPSToolkit workflows.
"""
function evolve!(psi, evolution)
  throw(MethodError(evolve!, (psi, evolution)))
end
