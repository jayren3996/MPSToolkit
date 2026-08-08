"""
    _lindblad_jump_list(jumps)

Normalize the `jumps` argument into a plain Julia vector of dense matrices.
"""
function _lindblad_jump_list(jumps::AbstractMatrix)
  return AbstractMatrix[jumps]
end

"""
    _lindblad_jump_list(jumps::AbstractVector)

Vector overload of [`_lindblad_jump_list`](@ref).
"""
function _lindblad_jump_list(jumps::AbstractVector)
  return collect(jumps)
end
