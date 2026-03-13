"""
    ChebyshevRescaling(center, halfwidth)

Store the affine map between physical frequencies `ω` and the rescaled Chebyshev
coordinate `x = (ω - center) / halfwidth`.
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
"""
struct SpectralFunction{TM<:AbstractVector{<:Real}, TK<:Union{Nothing, AbstractVector{<:Real}}}
  moments::TM
  kernel::TK
  rescaling::ChebyshevRescaling
end

"""
    spectral_function(moments; center=0.0, halfwidth=1.0, kernel=:jackson)

Build a callable `SpectralFunction` from a Chebyshev moment sequence.
"""
function spectral_function(moments::AbstractVector{<:Real}; center::Real=0.0, halfwidth::Real=1.0, kernel=:jackson)
  applied_kernel = _resolve_kernel(kernel, length(moments))
  return SpectralFunction(collect(float.(moments)), applied_kernel, ChebyshevRescaling(center, halfwidth))
end

function (spectrum::SpectralFunction)(ω::Real)
  x = (ω - spectrum.rescaling.center) / spectrum.rescaling.halfwidth
  abs(x) < 1 || return 0.0
  return reconstruct_chebyshev(x, spectrum.moments; kernel=spectrum.kernel) / spectrum.rescaling.halfwidth
end

function (spectrum::SpectralFunction)(ω::AbstractVector{<:Real})
  return [spectrum(ωi) for ωi in ω]
end
