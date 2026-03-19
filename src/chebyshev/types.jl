"""
    ChebyshevRescaling(center, halfwidth)

Store the affine map between physical frequencies `ω` and the rescaled Chebyshev
coordinate `x = (ω - center) / halfwidth`.

# Fields
- `center`: Physical frequency mapped to `x = 0`.
- `halfwidth`: Positive scale factor mapping the target energy window to `[-1, 1]`.
"""
struct ChebyshevRescaling
  center::Float64
  halfwidth::Float64

  function ChebyshevRescaling(center::Real, halfwidth::Real)
    halfwidth > 0 || throw(ArgumentError("Chebyshev halfwidth must be positive"))
    return new(Float64(center), Float64(halfwidth))
  end
end

"""
    SpectralFunction(moments, kernel, rescaling)

Represent a reconstructed spectral function from Chebyshev moments.
`moments[n + 1]` stores the moment `μ_n`.

# Fields
- `moments`: Stored Chebyshev moments in Julia's `1`-based indexing convention.
- `kernel`: Optional damping kernel already resolved to a numeric vector.
- `rescaling`: [`ChebyshevRescaling`](@ref) converting between physical and rescaled
  frequencies.
"""
struct SpectralFunction{TM<:AbstractVector{<:Real}, TK<:Union{Nothing, AbstractVector{<:Real}}}
  moments::TM
  kernel::TK
  rescaling::ChebyshevRescaling
end

"""
    spectral_function(moments; center=0.0, halfwidth=1.0, kernel=:jackson)

Build a callable [`SpectralFunction`](@ref) from a Chebyshev moment sequence.

# Arguments
- `moments`: Vector storing `μ_0, μ_1, ..., μ_{N-1}`.

# Keyword Arguments
- `center`: Physical frequency mapped to `x = 0`.
- `halfwidth`: Positive rescaling width.
- `kernel`: Kernel specification passed to the internal `_resolve_kernel` helper, such as
  `:jackson`, `nothing`, or an explicit weight vector.

# Returns
- A `SpectralFunction` object that can be called on scalars or vectors of frequencies.
"""
function spectral_function(moments::AbstractVector{<:Real}; center::Real=0.0, halfwidth::Real=1.0, kernel=:jackson)
  applied_kernel = _resolve_kernel(kernel, length(moments))
  return SpectralFunction(collect(float.(moments)), applied_kernel, ChebyshevRescaling(center, halfwidth))
end

"""
    (spectrum::SpectralFunction)(ω)

Evaluate a reconstructed spectral function at one physical frequency `ω`.
"""
function (spectrum::SpectralFunction)(ω::Real)
  x = (ω - spectrum.rescaling.center) / spectrum.rescaling.halfwidth
  abs(x) < 1 || return 0.0
  return reconstruct_chebyshev(x, spectrum.moments; kernel=spectrum.kernel) / spectrum.rescaling.halfwidth
end

"""
    (spectrum::SpectralFunction)(ωs)

Vectorized overload of [`SpectralFunction`](@ref) evaluation.
"""
function (spectrum::SpectralFunction)(ω::AbstractVector{<:Real})
  return [spectrum(ωi) for ωi in ω]
end
