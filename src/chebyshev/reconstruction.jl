"""
    jackson_damping(n, order)

Return the Jackson damping factor for Chebyshev index `n = 0, 1, ..., order - 1`.

# Arguments
- `n`: Zero-based Chebyshev index.
- `order`: Total number of moments in the truncated expansion.

# Returns
- The scalar Jackson damping factor for index `n`.
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

# Arguments
- `order`: Number of retained Chebyshev moments.

# Returns
- A dense vector whose `n + 1` entry stores the Jackson weight for moment `μ_n`.
"""
function jackson_kernel(order::Integer)
  order > 0 || throw(ArgumentError("Chebyshev order must be positive"))
  return [jackson_damping(n, order) for n in 0:(order - 1)]
end

"""
    reconstruct_chebyshev(x, moments; kernel=nothing)

Evaluate the Chebyshev reconstruction at rescaled frequency `x ∈ (-1, 1)`.
If `kernel` is omitted, the raw truncated series is used.

# Arguments
- `x`: Rescaled Chebyshev coordinate in `(-1, 1)`.
- `moments`: Moment vector in the convention `moments[n + 1] = μ_n`.

# Keyword Arguments
- `kernel`: Optional damping weights with the same length as `moments`.

# Returns
- The reconstructed spectral density at `x`.

# Notes
- The returned value already includes the Chebyshev weight factor `1 / (π * sqrt(1 - x^2))`.
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

"""
    _resolve_kernel(kernel, order)

Normalize a user-facing kernel specification into either `nothing` or a dense weight vector.
"""
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

"""
    _kernel_weight(kernel, index)

Return the weight associated with one moment index, defaulting to `1.0` when no kernel is present.
"""
function _kernel_weight(kernel, index::Integer)
  return isnothing(kernel) ? 1.0 : kernel[index]
end
