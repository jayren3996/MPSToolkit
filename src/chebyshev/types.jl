"""
    ChebyshevRescaling(center, halfwidth)

Store the affine map between physical frequencies `ω` and the rescaled Chebyshev
coordinate `x = (ω - center) / halfwidth`.

# Fields
- `center`: Physical frequency mapped to `x = 0`.
- `halfwidth`: Positive scale factor mapping the target energy window to `[-1, 1]`.

# Notes
- Choose `halfwidth` with a small safety margin (e.g. `(E_max - E_min) / (2 - ε)`) so that
  every eigenvalue maps strictly inside `(-1, 1)`. Spectral weight that maps to exactly `±1`
  lands on the band edge, where the Chebyshev measure `1 / (π √(1 - x²))` diverges and the
  reconstructed [`SpectralFunction`](@ref) returns `0.0`, silently dropping that weight.
"""
struct ChebyshevRescaling
  center::Float64
  halfwidth::Float64

  function ChebyshevRescaling(center::Real, halfwidth::Real)
    isfinite(center) || throw(ArgumentError("Chebyshev center must be finite"))
    isfinite(halfwidth) || throw(ArgumentError("Chebyshev halfwidth must be finite"))
    halfwidth > 0 || throw(ArgumentError("Chebyshev halfwidth must be positive"))
    return new(Float64(center), Float64(halfwidth))
  end
end

"""
    chebyshev_rescaling(emin, emax; margin=0.05)

Build the [`ChebyshevRescaling`](@ref) that maps the spectral window `[emin, emax]` into the
Chebyshev coordinate so the extreme eigenvalues land at `±(1 - margin)`, leaving a safety
buffer so no spectral weight reaches the band edge `±1`.

# Arguments
- `emin`, `emax`: Lower and upper bounds (or estimates) of the Hamiltonian spectrum.

# Keyword Arguments
- `margin`: Fractional safety buffer in `[0, 1)`. The extreme eigenvalues map to `±(1 - margin)`.

# Returns
- A [`ChebyshevRescaling`](@ref) with `center = (emax + emin) / 2` and
  `halfwidth = (emax - emin) / (2 * (1 - margin))`.
"""
function chebyshev_rescaling(emin::Real, emax::Real; margin::Real=0.05)
  emax > emin || throw(ArgumentError("chebyshev_rescaling requires emax > emin"))
  0 <= margin < 1 || throw(ArgumentError("chebyshev_rescaling requires 0 <= margin < 1"))
  center = (emax + emin) / 2
  halfwidth = (emax - emin) / (2 * (1 - margin))
  return ChebyshevRescaling(center, halfwidth)
end

"""
    rescale_hamiltonian(H::MPO, emin, emax; margin=0.05, cutoff=1e-14)

Rescale an MPO Hamiltonian into the Chebyshev band.

# Arguments
- `H`: Hamiltonian `MPO`.
- `emin`, `emax`: Spectral bounds (or estimates) of `H`.

# Keyword Arguments
- `margin`: Safety buffer forwarded to [`chebyshev_rescaling`](@ref).
- `cutoff`: Truncation cutoff used when subtracting the identity shift.

# Returns
- `(H_rescaled, rescaling)` where `H_rescaled = (H - center) / halfwidth` has its spectrum
  strictly inside `(-1, 1)`, and `rescaling::`[`ChebyshevRescaling`](@ref) inverts the map for
  [`spectral_function`](@ref). Feed `H_rescaled` to [`chebyshev_moments`](@ref) and reuse
  `rescaling.center` / `rescaling.halfwidth` when reconstructing the physical spectral function.
"""
function rescale_hamiltonian(H::MPO, emin::Real, emax::Real; margin::Real=0.05, cutoff::Real=1e-14)
  rescaling = chebyshev_rescaling(emin, emax; margin=margin)
  identity_mpo = MPO(firstsiteinds(H), "Id")
  shifted = add(H, -rescaling.center * identity_mpo; cutoff=cutoff)
  return shifted / rescaling.halfwidth, rescaling
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
  length(moments) > 0 || throw(ArgumentError("spectral_function requires at least one moment"))
  applied_kernel = _resolve_kernel(kernel, length(moments))
  return SpectralFunction(collect(float.(moments)), applied_kernel, ChebyshevRescaling(center, halfwidth))
end

"""
    (spectrum::SpectralFunction)(ω)

Evaluate a reconstructed spectral function at one physical frequency `ω`.

# Notes
- Frequencies that map outside the open interval `(-1, 1)` (i.e. `|x| ≥ 1`, including the band
  edges `x = ±1`) return `0.0`. Keep the spectrum strictly inside the band through the
  rescaling (see [`ChebyshevRescaling`](@ref)).
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
