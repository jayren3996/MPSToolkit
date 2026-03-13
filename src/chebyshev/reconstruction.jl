"""
    jackson_damping(n, order)

Return the Jackson damping factor for Chebyshev index `n = 0, 1, ..., order - 1`.
"""
function jackson_damping(n::Integer, order::Integer)
  0 <= n < order || throw(ArgumentError("Jackson index must satisfy 0 <= n < order"))
  θ = π / (order + 1)
  numerator = (order - n + 1) * cos(n * θ) + sin(n * θ) * cot(θ)
  return numerator / (order + 1)
end

"""
    jackson_kernel(order)

Return the full Jackson damping kernel for a Chebyshev series of length `order`.
"""
function jackson_kernel(order::Integer)
  order > 0 || throw(ArgumentError("Chebyshev order must be positive"))
  return [jackson_damping(n, order) for n in 0:(order - 1)]
end

"""
    reconstruct_chebyshev(x, moments; kernel=nothing)

Evaluate the Chebyshev reconstruction at rescaled frequency `x ∈ (-1, 1)`.
If `kernel` is omitted, the raw truncated series is used.
"""
function reconstruct_chebyshev(x::Real, moments::AbstractVector{<:Real}; kernel::Union{Nothing, AbstractVector{<:Real}}=nothing)
  abs(x) < 1 || throw(ArgumentError("Chebyshev reconstruction requires x in (-1, 1)"))
  length(moments) > 0 || throw(ArgumentError("Chebyshev reconstruction requires at least one moment"))
  isnothing(kernel) || length(kernel) == length(moments) || throw(ArgumentError("kernel length must match the number of moments"))

  value = _kernel_weight(kernel, 1) * moments[1]
  t_prev = one(float(x))
  if length(moments) == 1
    return value / (π * sqrt(1 - x^2))
  end

  t_curr = float(x)
  value += 2 * _kernel_weight(kernel, 2) * moments[2] * t_curr
  for n in 3:length(moments)
    t_next = 2 * x * t_curr - t_prev
    value += 2 * _kernel_weight(kernel, n) * moments[n] * t_next
    t_prev = t_curr
    t_curr = t_next
  end
  return value / (π * sqrt(1 - x^2))
end

function _resolve_kernel(kernel, order::Integer)
  if kernel === nothing
    return nothing
  elseif kernel === :jackson
    return jackson_kernel(order)
  elseif kernel isa AbstractVector
    length(kernel) == order || throw(ArgumentError("kernel length must match the number of moments"))
    return collect(float.(kernel))
  end
  throw(ArgumentError("unsupported Chebyshev kernel specification: $(kernel)"))
end

function _kernel_weight(kernel, index::Integer)
  return isnothing(kernel) ? 1.0 : kernel[index]
end
