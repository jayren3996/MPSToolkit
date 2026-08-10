# DMT cost across local dimension, bond dimension, and kernel options. Three tables:
#
#   1. PER-BOND KERNEL COST -- the `factorize` knob (Stage A: how the bond tensor is split, `:qr`
#      or `:svd`) crossed with the `truncation` knob (Stage B: how the doubly-projected complement
#      is truncated, `:dense` or `:random`). They are independent, so the full cross product is
#      measured. The bond factorization is also timed on its own, because that is the part `:qr`
#      replaces: the difference between its column and the step columns is what the rest of the
#      kernel costs, and under `:qr` the rest includes the `O(chi^2)` triangular products the SVD
#      path gets for free.
#   2. GATE VS TRUNCATION SPLIT -- one gate application followed by one DMT bond truncation,
#      timed separately, across `gate_maxdim`. This is the table that decides the `gate_maxdim`
#      default: it shows what the pre-truncation cap does and does not buy.
#   3. FULL SWEEP -- end-to-end `dmt_evolve!` wall clock across `gate_maxdim` and `truncation`,
#      the number a user actually experiences.
#
# `factorize` appears only in table 1: it is a kernel-internal knob with no observable effect and
# is deliberately not plumbed through `DMTOptions` / `DMTGateEvolution`, so `dmt_evolve!` always
# runs the `:qr` default.
#
# Three things the numbers depend on, all deliberate:
#   * In table 1 the state is gauged to the bond ONCE, outside the timed region, and the kernel
#     runs with `orthogonalize = false`. Re-gauging costs the same whatever the options are, and
#     including it would dilute the comparison. Tables 2 and 3 DO include re-gauging, because
#     there it is part of the cost the caller pays.
#   * `(nsites, bond)` is chosen per local dimension so the measured bond tensor really has the
#     intended shape. On too short a chain the link dimension is clipped by `d^(2 bond)` and the
#     large-`chi` cells silently measure a much smaller problem.
#   * Every table also sweeps `elt in (Float64, ComplexF64)`. A real-arithmetic DMT kernel path is
#     coming in a later change; this file must be able to measure it, and this `elt = ComplexF64`
#     baseline -- taken BEFORE that change, on the shipped kernel -- is what the real path will be
#     compared against. `elt = Float64` today still runs the shipped (complex-only-fast) kernel,
#     so it is not expected to be faster yet; it is here to be re-run once it is.
#
# Table 1 also carries a peak-memory column for `:dense` (bytes allocated by one bond truncation at
# the shipped `:qr` factorize default, via `@allocated` -- the idiom `test/test_dmt_kernel.jl`
# already uses for this same kernel call). `:dense` materializes and factorizes the full
# `chi' x chi'` complement, and halving that transient is the point of the real-arithmetic path:
# see `dmt.jl:39-43`'s ~7 GB `d = 4, preserve_diameter = 5, maxdim = 513` case.
#
# This box is shared. BLAS threads default to a fixed, modest 8 (override with the
# `MPSTOOLKIT_BENCH_BLAS_THREADS` env var) rather than half the machine's cores: the old
# `max(Sys.CPU_THREADS ÷ 2, 1)` (32 threads on this 64-core box) oversubscribes under any
# concurrent load and corrupts the smaller cells -- measured a physically impossible 12.7x
# real-vs-complex ratio on a 600x600 `gemm` against a true flop ratio of 4x. The report header
# below prints the thread count and `/proc/loadavg` so a reader can tell whether the box was quiet;
# a loadavg above ~8 means the numbers are junk, regardless of what they say.
#
# Run: julia --project=. dev/bench_dmt.jl    (was ~15 minutes on a quiet 64-core box before the
# `elt` axis; that axis roughly doubles the cell count, so budget accordingly, and re-time on a
# quiet box rather than trusting this comment)

using ITensors, ITensorMPS, LinearAlgebra, MPSToolkit, Printf, Random

Random.seed!(1)
BLAS.set_num_threads(parse(Int, get(ENV, "MPSTOOLKIT_BENCH_BLAS_THREADS", "8")))
println("BLAS threads = ", BLAS.get_num_threads())
println("loadavg (1m 5m 15m, running/total procs, last pid) = ", strip(read("/proc/loadavg", String)))

# d => (nsites, bond), with d^(2 (bond - 1)) >= the largest chi below on both sides of the cut.
const GEOMETRY = Dict(2 => (14, 7), 3 => (10, 5), 4 => (8, 4))
const CHIS = (400, 800, 1600)

"""
A random operator-space MPS at element type `elt`, already gauged to its benchmark bond, with the
shape of the bond tensor the kernel will factorize. Returns `(psi, bond, rows, cols)`.
"""
function bench_state(d, chi; elt=ComplexF64)
  nsites, bond = GEOMETRY[d]
  sites = operator_siteinds(nsites; d=d)
  psi = random_mps(elt, sites; linkdims=chi)
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

"""
Bytes allocated by one DMT bond truncation on `psi`, via `@allocated` -- the same idiom
`test/test_dmt_kernel.jl` already uses to characterize this exact kernel call ("one bond step
allocates a bounded multiple of the bond tensor"). `:dense` truncation materializes and factorizes
the full `chi' x chi'` complement (see `dmt.jl:39-43`), so this is the number a future
real-arithmetic kernel is expected to roughly halve; captured here as the pre-change baseline.
Warms up once (compile only) before the measured call, exactly like `bench_bond` above.
"""
function bench_bond_memory(psi, bond, maxdim; factorize, truncation)
  step!(state) = MPSToolkit._dmt_bond_truncate!(state, bond; maxdim=maxdim, cutoff=1e-14,
    factorize=factorize, truncation=truncation, orthogonalize=false)
  step!(copy(psi))                                  # warm up: compile, do not measure
  return @allocated step!(copy(psi))
end

# Warm up every branch -- both factorize x truncation and both elt -- on a tiny state, so the
# first real cell does not pay for compilation.
let
  for elt in (Float64, ComplexF64)
    psi, bond, _, _ = bench_state(2, 20; elt=elt)
    for factorize in (:svd, :qr)
      bench_factorize(psi, bond; factorize=factorize, reps=1)
      for truncation in (:dense, :random)
        bench_bond(psi, bond, 18; factorize=factorize, truncation=truncation, reps=1)
      end
    end
    bench_bond_memory(psi, bond, 18; factorize=:qr, truncation=:dense)
  end
end

println("\n== 1. per-bond kernel cost: the factorize and truncation knobs ==\n")
@printf("%-11s %-3s %-6s %-12s %-7s %25s %25s %25s %17s\n",
  "", "", "", "", "", "factorization (s)", "step, :dense (s)", "step, :random (s)", "")
@printf("%-11s %-3s %-6s %-12s %-7s %8s %8s %7s %8s %8s %7s %8s %8s %7s %17s\n",
  "elt", "d", "chi", "bond tensor", "maxdim",
  "svd", "qr", "gain", "svd", "qr", "gain", "svd", "qr", "gain", "mem :dense,qr (MB)")
for elt in (Float64, ComplexF64), d in (2, 3, 4), chi in CHIS
  maxdim = 2 * d^2 + 180
  psi, bond, rows, cols = bench_state(d, chi; elt=elt)
  factor_svd = bench_factorize(psi, bond; factorize=:svd)
  factor_qr = bench_factorize(psi, bond; factorize=:qr)
  steps = map((:dense, :random)) do truncation
    return (bench_bond(psi, bond, maxdim; factorize=:svd, truncation=truncation),
      bench_bond(psi, bond, maxdim; factorize=:qr, truncation=truncation))
  end
  # Memory is measured once, under the shipped `:qr` factorize default (see the file header: Stage
  # A's factorize knob has no observable effect and is not user-facing), so this is the peak-memory
  # number a real `dmt_evolve!` call actually pays.
  mem_dense = bench_bond_memory(psi, bond, maxdim; factorize=:qr, truncation=:dense)
  @printf("%-11s %-3d %-6d %-12s %-7d %8.3f %8.3f %6.2fx %8.3f %8.3f %6.2fx %8.3f %8.3f %6.2fx %17.1f\n",
    elt, d, chi, "$(rows)x$(cols)", maxdim, factor_svd, factor_qr, factor_svd / factor_qr,
    steps[1][1], steps[1][2], steps[1][1] / steps[1][2],
    steps[2][1], steps[2][2], steps[2][1] / steps[2][2], mem_dense / 2^20)
  flush(stdout)
end

# ---------------------------------------------------------------------------------------------
# Tables 2 and 3: the `gate_maxdim` and `truncation` knobs.
#
# `gate_maxdim` caps the bond during raw gate application, BEFORE DMT sees it. The cap is applied
# by a plain SVD, which discards the smallest singular values -- exactly the error DMT exists to
# avoid -- so the question these tables answer is what the cap actually buys.
#
# The arithmetic worth having in mind while reading them: a two-site gate inflates the bond from
# `chi` to `min(chi d^2, cap)`, because the two-site block is `(chi d^2) x (d^2 chi)`. The old
# default `max(16 maxdim, 64)` therefore caps nothing at all whenever `d^2 <= 16`, i.e. for every
# `d <= 4`: the shipped behaviour at those local dimensions was ALREADY exact gate application.
# `gate_caps` below includes that old default explicitly so the claim is measured, not asserted.
# ---------------------------------------------------------------------------------------------

# d => (nsites, bond) for the split table: `bond` is far enough from both edges that
# `d^(2 min(bond, nsites - bond))` exceeds the largest inflated bond below (`d^2 * maxdim`).
const SPLIT_GEOMETRY = Dict(2 => (12, 6), 3 => (8, 4))
# d => nsites for the sweep table. Shorter, because a sweep pays for every bond; the edge bonds
# are clipped by `d^(2 b)` there, exactly as they are in a real run.
const SWEEP_GEOMETRY = Dict(2 => 10, 3 => 8)
const MAXDIM_OFFSETS = (100, 300)

"""
The `gate_maxdim` settings compared in tables 2 and 3, as `(value, label)` pairs: no cap at all,
the pre-Task-10 default, and a cap that really does bite.
"""
gate_caps(maxdim) = ((0, "0 (exact)"), (max(16 * maxdim, 64), "max(16 maxdim, 64) [old]"),
  (2 * maxdim, "2 maxdim"))

"""
The bond gate used by tables 2 and 3: nearest-neighbour `XX + ZZ` at `dt = 0.1`, in the operator
basis of local dimension `d`, narrowed to `elt`.

`operator_gate_from_hamiltonian` always computes internally at `ComplexF64` (`gates.jl` widens any
input via `Matrix{ComplexF64}(h)` before exponentiating), so the returned gate is real only if it
is narrowed here -- otherwise a state built at `elt = Float64` would still be gated by a complex
tensor and `product` would promote the whole bond back to `ComplexF64`, same trap as the melt in
`bench_sweep` below. The narrowing is exact, not approximate: for ANY unitary `A`, `G[a,b] =
tr(P_a (A P_b A'))` is real whenever the basis `P` is Hermitian (`P_a` and `A P_b A'` are both
Hermitian, and the trace of a product of two Hermitian matrices is always real), so a gate built in
a Hermitian operator basis is real regardless of `A`, not only for a real Hamiltonian. The assert
below is a tripwire, not a tolerance to tune: the measured residue across `d = 2, 3, 4` is
1e-18 .. 1.2e-16 (rounding noise from `exp` and the dense products), eight-plus orders below it.
"""
function bench_gate(d; elt=ComplexF64)
  sx = d == 2 ? elt[0 1; 1 0] / 2 : elt[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
  sz = d == 2 ? elt[1 0; 0 -1] / 2 : elt[1 0 0; 0 0 0; 0 0 -1]
  gate = operator_gate_from_hamiltonian(kron(sx, sx) + kron(sz, sz), 0.1; d=d)
  elt <: Real || return gate
  tol = 1e-8 * maximum(abs, gate)
  @assert maximum(abs, imag.(gate)) <= tol "bench_gate: operator_gate is not real to within " *
    "tolerance at d = $(d); narrowing to $(elt) would be lossy"
  return real.(gate)
end

"""
A trace-ful operator-space state at element type `elt` whose bonds sit at `chi`: the identity plus
a random component, which is what a melt looks like to the kernel once it has filled its bond
budget.
"""
function bench_operator_state(sites, chi; elt=ComplexF64)
  # `operator_basis_state`'s default `coefficient = 1.0` builds a real Float64 MPS regardless of
  # `elt` (see states.jl), so plain `add` already promotes correctly here: Float64 + elt gives elt
  # exactly (unlike `bench_sweep`'s melt below, which needs explicit narrowing -- see its comment).
  return add(operator_basis_state(sites, fill(1, length(sites))),
    0.3 * random_mps(elt, sites; linkdims=chi); maxdim=chi, cutoff=0.0)
end

"""
Best of `reps` for one gate application and the DMT bond truncation that follows it, timed
SEPARATELY. Returns `(gate_seconds, dmt_seconds, inflated_bond_dimension)`.

The gate is timed at whatever `gate_maxdim` says; the truncation is then timed on the bond the
gate left behind, so the two columns show where the cap moves work rather than only the total.
"""
function bench_split(d, maxdim, gate_maxdim; truncation, reps=2, elt=ComplexF64)
  nsites, bond = SPLIT_GEOMETRY[d]
  sites = operator_siteinds(nsites; d=d)
  base = bench_operator_state(sites, maxdim; elt=elt)
  gate = bench_gate(d; elt=elt)
  gate_best, dmt_best, inflated = Inf, Inf, 0
  for _ in 1:reps
    work = copy(base)
    orthogonalize!(work, bond)     # excluded: the same cost whatever the knobs say
    gate_best = min(gate_best,
      @elapsed tebd_evolve!(work, gate, bond; maxdim=gate_maxdim, cutoff=0.0))
    inflated = dim(linkind(work, bond))
    dmt_best = min(dmt_best, @elapsed MPSToolkit._dmt_bond_truncate!(work, bond; maxdim=maxdim,
      cutoff=1e-14, truncation=truncation, orthogonalize=true))
  end
  return gate_best, dmt_best, inflated
end

"""
Narrow an operator-space `MPS` from `ComplexF64` to `Float64` when `elt` calls for it; a no-op at
`elt = ComplexF64`.

`operator_product_state` and `operator_local_sum_state` (`states.jl`) always build `ComplexF64`
tensors, independent of the caller's operator matrices, so the melt `bench_sweep` builds from them
below is always `ComplexF64` regardless of `elt` -- and `add`ing it to an `elt`-typed random MPS
would silently promote the whole state back to `ComplexF64` at `elt = Float64` (`add(real, real) ->
Float64` but `add(complex, real) -> ComplexF64`). The narrowing here is EXACT, not approximate:
`identity_matrix` and `sz` are both real Hermitian, and the trace of a real symmetric operator
against any purely-imaginary operator-basis element is identically zero (a standard trace
identity: for antisymmetric real `A` and symmetric real `O`, `tr(AO) = tr((AO)^T) = tr(O A^T) =
-tr(AO)`, so `tr(AO) = 0`), so the discarded imaginary part is exactly `0.0`, not merely small --
confirmed numerically at dev time on the actual tensors this narrows.
"""
melt_as_elt(elt, psi::MPS) = elt <: Real ? MPS([real(t) for t in psi]) : psi

"""
Wall clock for `nstep` full forward-and-reverse DMT sweeps, from a domain-wall melt already at
its bond budget. This is the end-to-end number, re-gauging and all.
"""
function bench_sweep(d, maxdim, gate_maxdim; truncation, nstep=1, elt=ComplexF64)
  nsites = SWEEP_GEOMETRY[d]
  sites = operator_siteinds(nsites; d=d)
  identity_matrix = Matrix{elt}(I, d, d)
  sz = d == 2 ? elt[1 0; 0 -1] / 2 : elt[1 0 0; 0 0 0; 0 0 -1]
  melt = add(operator_product_state(sites, fill(identity_matrix, nsites)),
    operator_local_sum_state(sites, sz, [j <= nsites ÷ 2 ? -0.25 : 0.25 for j in 1:nsites]);
    maxdim=8, cutoff=0.0)
  melt = melt_as_elt(elt, melt)
  # Bring the state up to its bond budget first, so the timed sweep is the steady-state cost
  # rather than the cost of a nearly-product state that has not filled `maxdim` yet.
  rho = add(melt, 0.05 * random_mps(elt, sites; linkdims=maxdim); maxdim=maxdim, cutoff=0.0)
  schedule = collect(1:(nsites - 1))
  evo = DMTGateEvolution(bench_gate(d; elt=elt), 0.1; schedule=schedule,
    reverse_schedule=reverse(schedule), nstep=nstep, maxdim=maxdim, gate_maxdim=gate_maxdim,
    cutoff=1e-14, truncation=truncation, normalize=false)
  return @elapsed dmt_evolve!(rho, evo)
end

let d = 2, maxdim = 2 * d^2 + 20               # warm up both tables on a small cell, both elt
  for elt in (Float64, ComplexF64), truncation in (:dense, :random)
    bench_split(d, maxdim, 0; truncation=truncation, reps=1, elt=elt)
    bench_sweep(d, maxdim, 0; truncation=truncation, elt=elt)
  end
end

println("\n== 2. gate application vs DMT truncation, across gate_maxdim ==\n")
@printf("%-11s %-3s %-7s %-22s %-6s %-8s %8s %8s %9s\n",
  "elt", "d", "maxdim", "gate_maxdim", "bond", "trunc", "gate(s)", "dmt(s)", "total(s)")
for elt in (Float64, ComplexF64), d in (2, 3), offset in MAXDIM_OFFSETS
  maxdim = 2 * d^2 + offset
  for (cap, label) in gate_caps(maxdim), truncation in (:dense, :random)
    gate_seconds, dmt_seconds, inflated = bench_split(d, maxdim, cap; truncation=truncation, elt=elt)
    @printf("%-11s %-3d %-7d %-22s %-6d %-8s %8.3f %8.3f %9.3f\n",
      elt, d, maxdim, label, inflated, truncation, gate_seconds, dmt_seconds,
      gate_seconds + dmt_seconds)
    flush(stdout)
  end
end

println("\n== 3. full dmt_evolve! sweep, across gate_maxdim and truncation ==\n")
@printf("%-11s %-3s %-7s %-7s %-22s %10s %10s %7s\n",
  "elt", "d", "nsites", "maxdim", "gate_maxdim", ":dense(s)", ":random(s)", "gain")
for elt in (Float64, ComplexF64), d in (2, 3), offset in MAXDIM_OFFSETS
  maxdim = 2 * d^2 + offset
  for (cap, label) in gate_caps(maxdim)
    dense_seconds = bench_sweep(d, maxdim, cap; truncation=:dense, elt=elt)
    random_seconds = bench_sweep(d, maxdim, cap; truncation=:random, elt=elt)
    @printf("%-11s %-3d %-7d %-7d %-22s %10.3f %10.3f %6.2fx\n",
      elt, d, SWEEP_GEOMETRY[d], maxdim, label, dense_seconds, random_seconds,
      dense_seconds / random_seconds)
    flush(stdout)
  end
end
