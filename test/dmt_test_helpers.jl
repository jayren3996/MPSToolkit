using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random

"""
Dense random probes, `nrandom` per `(width, start)` window, as `(start, matrix)` pairs.

# Notes
- A random dense operator has generic components on **every** one of the `d^(2 width)` basis
  products of its window, so these reach the part of the local operator space a capped product
  sweep cannot see. That matters at `d > 2`, where enumerating the complete product basis on a
  wide window is unaffordable and [`diameter_probes`](@ref) has to cap it.
- Scaled to the Frobenius norm of a basis product (`sqrt(d^width)`) so a random probe does not
  dominate the shared denominator in [`preservation_error`](@ref).
- `rng` is seeded independently of the caller's stream, so a probe set is reproducible no matter
  where in a test file it is built.
"""
function random_window_probes(nsites::Int, d::Int, widths; nrandom::Int=3,
  rng::AbstractRNG=MersenneTwister(20260808))
  probes = Tuple{Int,Matrix{ComplexF64}}[]
  for width in widths, start in 1:(nsites - width + 1), _ in 1:nrandom
    op = randn(rng, ComplexF64, d^width, d^width)
    push!(probes, (start, ComplexF64.(op .* (sqrt(Float64(d)^width) / norm(op)))))
  end
  return probes
end

"""
All dense probes of diameter <= `diameter` on an `nsites` chain of local dimension `d`,
as `(start, matrix)` pairs suitable for `operator_expectation_profile`.

The structured product probes come first, then the `nrandom` dense random probes per window, so
a caller that wants to score the two families separately can split on the count returned by
`diameter_probes(...; nrandom = 0)`.

# Keyword Arguments
- `widths`: which support widths to emit; defaults to every width up to `diameter`.
- `full_basis_width`: widths up to this value enumerate the **complete** `d^2`-element onsite
  basis on every site of the window. A capped basis would only prove the guarantee on the
  `span{P_1..P_cap}` corner of the local operator space, which is strictly weaker than what DMT
  claims -- at `d > 2` a kernel protecting only that corner would pass a capped sweep.
- `cap`: onsite basis elements used for wider windows, where the complete product basis is
  unaffordable. `cap = 1` emits no product probes at all, leaving only the random ones.
- `nrandom`: dense random probes per window; see [`random_window_probes`](@ref).
- `rng`: forwarded to [`random_window_probes`](@ref).
"""
function diameter_probes(nsites::Int, d::Int, diameter::Int; widths=1:diameter,
  full_basis_width::Int=2, cap::Int=4, nrandom::Int=3,
  rng::AbstractRNG=MersenneTwister(20260808))
  basis = operator_basis_matrices(d)
  probes = Tuple{Int,Matrix{ComplexF64}}[]
  for width in widths, start in 1:(nsites - width + 1)
    # locals[1] is the identity; the rest are non-identity basis elements. Every combination
    # over the window is probed, so a middle site is not silently left as identity -- the
    # guarantee covers ALL operators of this diameter, so the test must too.
    count = width <= full_basis_width ? length(basis) : min(length(basis), cap)
    locals = [basis[k] * sqrt(d) for k in 1:count]
    for labels in Iterators.product(ntuple(_ -> 1:count, width)...)
      all(isequal(1), labels) && continue          # skip the pure-identity probe (that is the trace)
      op = locals[labels[1]]
      for offset in 2:width
        op = kron(op, locals[labels[offset]])
      end
      push!(probes, (start, ComplexF64.(op)))
    end
  end
  append!(probes, random_window_probes(nsites, d, widths; nrandom=nrandom, rng=rng))
  return probes
end

"""
Whether a probe of support width `span` starting at site `start` falls inside the DMT guarantee
for a truncation at `bond` protecting `radius` sites per side.

# Notes
- The truncation reproduces the bond matrix exactly on the whole protected left row space and on
  the whole protected right column space, so a probe is preserved as soon as **either** side of
  its support fits inside the `radius` protected sites on that side of the cut. That is why a
  probe lying entirely on one side of the cut is preserved at any width (its other side is then
  the trace direction), and why the guarantee runs out exactly at width `2 radius + 2` -- the
  first width able to put more than `radius` sites on both sides at once.
"""
function guarantee_covers(start::Int, span::Int, bond::Int, radius::Int)
  left = clamp(bond - start + 1, 0, span)
  return left <= radius || span - left <= radius
end

"""
Support width in sites of a probe matrix at local dimension `d`.
"""
probe_span(op::AbstractMatrix, d::Int) = round(Int, log(size(op, 1)) / log(d))

"""
Relative sup-norm change of a profile of expectation values under truncation.
"""
function preservation_error(before, after)
  scale = maximum(abs, before)
  scale == 0 && return maximum(abs, after)
  return maximum(abs, after .- before) / scale
end
