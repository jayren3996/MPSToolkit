using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit

"""
All dense probes of diameter <= `diameter` on an `nsites` chain of local dimension `d`,
as `(start, matrix)` pairs suitable for `operator_expectation_profile`.
"""
function diameter_probes(nsites::Int, d::Int, diameter::Int)
  basis = operator_basis_matrices(d)
  # locals[1] is the identity; the rest are a few non-identity basis elements. Every
  # combination over the window is probed, so a middle site is not silently left as identity
  # -- the guarantee covers ALL operators of this diameter, so the test must too.
  locals = [basis[k] * sqrt(d) for k in 1:min(length(basis), 4)]
  probes = Tuple{Int,Matrix{ComplexF64}}[]
  for width in 1:diameter, start in 1:(nsites - width + 1)
    for labels in Iterators.product(ntuple(_ -> eachindex(locals), width)...)
      all(isequal(1), labels) && continue          # skip the pure-identity probe (that is the trace)
      op = locals[labels[1]]
      for offset in 2:width
        op = kron(op, locals[labels[offset]])
      end
      push!(probes, (start, ComplexF64.(op)))
    end
  end
  return probes
end

"""
Relative sup-norm change of a profile of expectation values under truncation.
"""
function preservation_error(before, after)
  scale = maximum(abs, before)
  scale == 0 && return maximum(abs, after)
  return maximum(abs, after .- before) / scale
end
