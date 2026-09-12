using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Statistics

const NSITES = parse(Int, get(ENV, "MPSTOOLKIT_BOND_BENCH_NSITES", "8"))
const LINKDIM = parse(Int, get(ENV, "MPSTOOLKIT_BOND_BENCH_LINKDIM", "20"))
const MAXDIM = parse(Int, get(ENV, "MPSTOOLKIT_BOND_BENCH_MAXDIM", "12"))
const REPEATS = parse(Int, get(ENV, "MPSTOOLKIT_BOND_BENCH_REPEATS", "7"))
const BLAS_THREADS = parse(Int, get(ENV, "MPSTOOLKIT_BOND_BENCH_BLAS_THREADS", "1"))
const DURATION = parse(Float64, get(ENV, "MPSTOOLKIT_BOND_BENCH_DURATION", "0.05"))

BLAS.set_num_threads(BLAS_THREADS)

function git_output(args...)
  try
    return readchomp(`git -C $(pkgdir(MPSToolkit)) $(args)`)
  catch
    return "unavailable"
  end
end

function benchmark_backend(base, gate, bond, backend)
  warmup = copy(base)
  dmt_step!(warmup, gate, bond; maxdim=MAXDIM, cutoff=0.0, gate_backend=backend)
  times = Float64[]
  allocations = Int[]
  gc_times = Float64[]
  for _ in 1:REPEATS
    state = copy(base)
    GC.gc()
    sample = @timed dmt_step!(state, gate, bond;
      maxdim=MAXDIM, cutoff=0.0, gate_backend=backend)
    push!(times, sample.time)
    push!(allocations, sample.bytes)
    push!(gc_times, sample.gctime)
  end
  return (wall=median(times), spread=median(abs.(times .- median(times))),
    allocated=median(allocations), gc=median(gc_times))
end

function center_diagnostics(base, gate, bond)
  state = copy(base)
  orthogonalize!(state, bond)
  sites = [siteind(state, bond), siteind(state, bond + 1)]
  block = noprime(MPSToolkit._dense_local_operator(sites, gate) * state[bond] * state[bond + 1])
  previous = linkind(state, bond - 1)
  following = linkind(state, bond + 1)
  left = isnothing(previous) ? (sites[1],) : (previous, sites[1])
  right = isnothing(following) ? (sites[2],) : (sites[2], following)
  center = reshape(array(block, left..., right...), prod(dim.(left)), prod(dim.(right)))
  values = svdvals(center)
  rank_tolerance = isempty(values) ? 0.0 : maximum(size(center)) * eps(Float64) * values[1]
  return (size(center)..., count(>(rank_tolerance), values))
end

function profile_direct_stages(base, gate, bond)
  state = copy(base)
  gauge = @elapsed begin
    orthogonalize!(state, bond)
    MPSToolkit._direct_canonical_residual(state, bond)
  end
  sites = (siteind(state, bond), siteind(state, bond + 1))
  previous = linkind(state, bond - 1)
  following = linkind(state, bond + 1)
  left_inds = isnothing(previous) ? (sites[1],) : (previous, sites[1])
  right_inds = isnothing(following) ? (sites[2],) : (sites[2], following)
  center = nothing
  gate_center = @elapsed begin
    block = noprime(MPSToolkit._dense_local_operator(collect(sites), gate) *
      state[bond] * state[bond + 1])
    center = reshape(array(block, left_inds..., right_inds...), prod(dim.(left_inds)),
      prod(dim.(right_inds)))
  end
  protected_left = protected_right = nothing
  protected = @elapsed begin
    left_env = MPSToolkit._left_identity_environment(state, bond - 1)
    right_env = MPSToolkit._right_identity_environment(state, bond + 2)
    left_ops = MPSToolkit._preserved_operator_tensors([sites[1]], 2, nothing)
    right_ops = MPSToolkit._preserved_operator_tensors([sites[2]], 2, nothing)
    protected_left = conj(MPSToolkit._direct_protected_columns(left_env, left_inds,
      left_ops, eltype(center)))
    protected_right = MPSToolkit._direct_protected_columns(right_env, right_inds,
      right_ops, eltype(center))
  end
  ql = MPSToolkit._protected_basis(protected_left, eltype(center))
  qr = MPSToolkit._protected_basis(protected_right, eltype(center))
  q0 = MPSToolkit._unit_direction(protected_left[:, 1])
  r0 = MPSToolkit._unit_direction(protected_right[:, 1])
  a, b, _ = MPSToolkit._dmt_connector(center, q0, r0, eltype(center))
  ops = MPSToolkit._dmt_complement_ops(center, a, b, ql, qr)
  budget = max(MPSToolkit._dmt_complement_budget(MAXDIM, size(ql, 2), size(qr, 2)), 1)
  uc = sc = vc = nothing
  complement = @elapsed begin
    uc, sc, vc = MPSToolkit._truncated_svd(ops.mul, ops.adj, size(center)..., budget,
      eltype(center); mode=:dense, dense=ops.dense)
  end
  factor_left = hcat(a, ql, ops.BQRc, uc * Diagonal(sc))
  factor_right = hcat(conj(b), ops.QLtB', qr, vc)
  new_u = new_s = new_v = nothing
  refactor = @elapsed begin
    new_u, new_s, new_v = MPSToolkit._dmt_refactor(factor_left, factor_right, MAXDIM,
      MPSToolkit._dmt_refactor_tolerance(factor_left))
  end
  reconstruction = @elapsed begin
    new_link = Index(length(new_s), "Link,l=$(bond)")
    state[bond] = ITensor(reshape(new_u, dim.(left_inds)..., length(new_s)),
      left_inds..., new_link)
    state[bond + 1] = ITensor(reshape(Diagonal(new_s) * new_v', length(new_s),
      dim.(right_inds)...), dag(new_link), right_inds...)
  end
  return (gauge, gate_center, protected, complement, refactor, reconstruction)
end

function benchmark_model(name, hamiltonian)
  Random.seed!(20260910)
  sites = pauli_siteinds(NSITES)
  identity2 = Matrix{ComplexF64}(I, 2, 2)
  base = add(operator_product_state(sites, fill(identity2, NSITES)),
    0.1 * random_mps(sites; linkdims=LINKDIM); maxdim=LINKDIM, cutoff=0.0)
  bond = NSITES ÷ 2
  gate = pauli_gate_from_hamiltonian(hamiltonian, DURATION)
  original = copy(base)
  copies = [copy(base) for _ in 1:2]
  dmt_step!(copies[1], gate, bond; maxdim=MAXDIM, cutoff=0.0, gate_backend=:qr)
  norm(base - original) == 0 || error("benchmark state copies share mutable tensor data")
  nleft, nright, gated_rank = center_diagnostics(base, gate, bond)
  results = Dict(backend => benchmark_backend(base, gate, bond, backend)
    for backend in (:product, :qr, :fused, :direct))
  for backend in (:product, :qr, :fused, :direct)
    result = results[backend]
    println("| $(name) | `:$(backend)` | $(nleft)x$(nright) | $(gated_rank) | " *
      "$(round(1e3result.wall; digits=3)) ± $(round(1e3result.spread; digits=3)) ms | " *
      "$(round(result.allocated / 2^20; digits=3)) MiB | $(round(1e3result.gc; digits=3)) ms |")
  end
  profile_direct_stages(base, gate, bond)
  stages = [profile_direct_stages(base, gate, bond) for _ in 1:REPEATS]
  labels = ("gauge/canonical", "gate/center", "protected", "complement SVD",
    "refactor", "reconstruction")
  println("direct stages $(name): " * join(("$(labels[index])=$(round(1e3 * median(
    getindex.(stages, index)); digits=3)) ms" for index in eachindex(labels)), ", "))
  return results
end

println("# Generic two-site DMT backend benchmark")
println()
println("- Commit: `$(git_output("rev-parse", "HEAD"))`")
println("- Dirty worktree: `$(git_output("status", "--porcelain") == "" ? "false" : "true")`")
println("- Julia: `$(VERSION)`")
println("- CPU: `$(Sys.CPU_NAME)`")
println("- BLAS threads: `$(BLAS.get_num_threads())`")
println("- Input: N=$(NSITES), generated real operator MPS linkdim=$(LINKDIM), " *
  "center bond=$(NSITES ÷ 2), maxdim=$(MAXDIM), duration=$(DURATION), cutoff=0")
println("- Samples: $(REPEATS) separately copied states after one warm-up; medians reported")
println("- Allocation is cumulative `@allocated` volume for `dmt_step!`, not peak RSS")
println()
println("| model | backend | nL x nR | gated rank | median wall ± MAD | median allocation | median GC |")
println("| --- | --- | ---: | ---: | ---: | ---: | ---: |")

xxz = spinhalf_xyz_bond_hamiltonian(Jx=1.0, Jy=1.0, Jz=0.6)
mfi = spinhalf_mixed_field_ising_bond_hamiltonian(NSITES, NSITES ÷ 2;
  J=0.7, gx=1.1, gz=-0.3)
benchmark_model("XXZ", xxz)
benchmark_model("MFI bulk", mfi)
