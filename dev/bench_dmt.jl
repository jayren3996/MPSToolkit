# Per-bond DMT cost across local dimension, bond dimension, and kernel options.
#
# `factorize` (Stage A: how the bond tensor is split, `:qr` or `:svd`) and `truncation` (how the
# doubly-projected complement is truncated, `:dense` or `:random`) are independent knobs, so the
# full cross product is measured and every cell reports the `:qr` speedup at fixed `truncation`.
# The bond factorization is also timed on its own, because that is the part `:qr` replaces: the
# difference between its column and the step columns is what the rest of the kernel costs, and
# under `:qr` the rest includes the `O(chi^2)` triangular products the SVD path gets for free.
#
# Two things the numbers depend on, both deliberate:
#   * The state is gauged to the bond ONCE, outside the timed region, and the kernel runs with
#     `orthogonalize = false`. Re-gauging costs the same whatever the options are, and including
#     it would dilute the comparison.
#   * `(nsites, bond)` is chosen per local dimension so the measured bond tensor really is
#     `(chi d^2) x chi`. On too short a chain the link dimension is clipped by `d^(2 bond)` and
#     the large-`chi` cells silently measure a much smaller problem.
#
# Run: julia --project=. dev/bench_dmt.jl

using ITensors, ITensorMPS, LinearAlgebra, MPSToolkit, Printf, Random

Random.seed!(1)
BLAS.set_num_threads(max(Sys.CPU_THREADS ÷ 2, 1))
println("BLAS threads = ", BLAS.get_num_threads())

# d => (nsites, bond), with d^(2 (bond - 1)) >= the largest chi below on both sides of the cut.
const GEOMETRY = Dict(2 => (14, 7), 3 => (10, 5), 4 => (8, 4))
const CHIS = (400, 800, 1600)

"""
A random operator-space MPS already gauged to its benchmark bond, with the shape of the bond
tensor the kernel will factorize. Returns `(psi, bond, rows, cols)`.
"""
function bench_state(d, chi)
  nsites, bond = GEOMETRY[d]
  sites = operator_siteinds(nsites; d=d)
  psi = random_mps(ComplexF64, sites; linkdims=chi)
  normalize!(psi)
  orthogonalize!(psi, bond)
  return psi, bond, dim(linkind(psi, bond - 1)) * d^2, dim(linkind(psi, bond))
end

"""
Best of `reps` wall-clock times for the bond factorization alone -- Stage A in isolation.
"""
function bench_factorize(psi, bond; factorize, reps=3)
  previous_link = linkind(psi, bond - 1)
  left_inds = isnothing(previous_link) ? (siteind(psi, bond),) : (previous_link, siteind(psi, bond))
  best = Inf
  for _ in 1:reps
    best = min(best, @elapsed MPSToolkit._dmt_bond_factorize(psi, bond, left_inds;
      factorize=factorize))
  end
  return best
end

"""
Best of `reps` wall-clock times for one DMT bond truncation on `psi`.
"""
function bench_bond(psi, bond, maxdim; factorize, truncation, reps=3)
  best = Inf
  for _ in 1:reps
    # The kernel replaces `psi[bond]` and `psi[bond + 1]` rather than mutating them in place, so
    # a shallow copy is enough to give every repetition the same untouched input.
    work = copy(psi)
    elapsed = @elapsed MPSToolkit._dmt_bond_truncate!(work, bond; maxdim=maxdim, cutoff=1e-14,
      factorize=factorize, truncation=truncation, orthogonalize=false)
    best = min(best, elapsed)
  end
  return best
end

# Warm up every branch on a tiny state, so the first real cell does not pay for compilation.
let warmup = bench_state(2, 20)
  for factorize in (:svd, :qr)
    bench_factorize(warmup[1], warmup[2]; factorize=factorize, reps=1)
    for truncation in (:dense, :random)
      bench_bond(warmup[1], warmup[2], 18; factorize=factorize, truncation=truncation, reps=1)
    end
  end
end

@printf("%-3s %-6s %-12s %-7s %25s %25s %25s\n",
  "", "", "", "", "factorization (s)", "step, :dense (s)", "step, :random (s)")
@printf("%-3s %-6s %-12s %-7s %8s %8s %7s %8s %8s %7s %8s %8s %7s\n",
  "d", "chi", "bond tensor", "maxdim",
  "svd", "qr", "gain", "svd", "qr", "gain", "svd", "qr", "gain")
for d in (2, 3, 4), chi in CHIS
  maxdim = 2 * d^2 + 180
  psi, bond, rows, cols = bench_state(d, chi)
  factor_svd = bench_factorize(psi, bond; factorize=:svd)
  factor_qr = bench_factorize(psi, bond; factorize=:qr)
  steps = map((:dense, :random)) do truncation
    return (bench_bond(psi, bond, maxdim; factorize=:svd, truncation=truncation),
      bench_bond(psi, bond, maxdim; factorize=:qr, truncation=truncation))
  end
  @printf("%-3d %-6d %-12s %-7d %8.3f %8.3f %6.2fx %8.3f %8.3f %6.2fx %8.3f %8.3f %6.2fx\n",
    d, chi, "$(rows)x$(cols)", maxdim, factor_svd, factor_qr, factor_svd / factor_qr,
    steps[1][1], steps[1][2], steps[1][1] / steps[1][2],
    steps[2][1], steps[2][2], steps[2][1] / steps[2][2])
  flush(stdout)
end
