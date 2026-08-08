using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Test

include("dmt_test_helpers.jl")

@testset "DMT at higher local dimension" begin
  Random.seed!(20260808)

  @testset "the fused protected block keeps the identity in column 1" begin
    # The kernel reads the trace direction off column 1 of the protected block, which is fused
    # from `radius` site indices by a `combiner`. `combiner` may fuse in any order; what the
    # kernel actually needs is the weaker property pinned here -- that whatever the order, the
    # all-identity multi-index still lands in column 1, so column 1 equals the identity-capped
    # contraction. If this ever stops holding, `q0`/`r0` must be built from an explicit identity
    # cap instead of by column index. The fused width is also asserted, because it is the
    # `d^(2 radius)` the budget arithmetic assumes.
    for d in (2, 3, 4), radius in (1, 2)
      sites = operator_siteinds(radius; d=d)
      link = Index(3, "Link,probe")
      block = random_itensor(ComplexF64, link, sites...)
      fuse = combiner(sites...)
      fused = matrix(block * fuse, link, combinedind(fuse))
      @test size(fused, 2) == d^(2 * radius)
      capped = block
      for site in sites
        capped *= MPSToolkit._identity_env(site)
      end
      @test fused[:, 1] ≈ vector(capped, link)
    end
  end

  @testset "diameter <= 3 preserved at d = 3 and d = 4" begin
    for d in (3, 4)
      nsites, chi = 6, 60
      maxdim = 2 * d^2 + 12
      sites = operator_siteinds(nsites; d=d)
      probes = diameter_probes(nsites, d, 3)
      for (label, noise) in (("Hermitian (real coefficients)", random_mps(sites; linkdims=chi)),
        ("non-Hermitian (complex)", random_mps(ComplexF64, sites; linkdims=chi)))
        rho = add(operator_basis_state(sites, fill(1, nsites)), 0.3 * noise;
          maxdim=chi + 1, cutoff=0.0)
        @test dim(linkind(rho, 3)) > maxdim          # the truncation really fires
        before = operator_expectation_profile(rho, probes; normalize=false)
        trace_before = operator_trace(rho)
        MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim=maxdim, cutoff=0.0)
        sweep_error = preservation_error(before,
          operator_expectation_profile(rho, probes; normalize=false))
        @info "d = $(d), $(label): $(length(probes)) probes, preservation error $(sweep_error)"
        @test dim(linkind(rho, 3)) <= maxdim
        @test sweep_error < 1e-11
        @test abs(operator_trace(rho) - trace_before) <= 1e-11 * abs(trace_before)
      end
    end
  end

  @testset "budget floor scales with d" begin
    # 2 d^(preserve_diameter - 1) + 1: 19 at d = 3 and 33 at d = 4 for diameter 3, 33 at d = 2
    # and 163 at d = 3 for diameter 5. The message must name the local dimension it applies to,
    # because that is the only way a caller learns which floor bit them.
    for (d, diameter, floor_value) in ((3, 3, 19), (4, 3, 33), (2, 5, 33), (3, 5, 163))
      sites = operator_siteinds(6; d=d)
      rho = random_mps(ComplexF64, sites; linkdims=20)
      failure = try
        MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim=floor_value - 1, cutoff=0.0,
          preserve_diameter=diameter)
        nothing
      catch exception
        exception
      end
      @test failure isa ArgumentError
      message = sprint(showerror, failure)
      @test occursin("local dimension d = $(d)", message)
      @test occursin("+ 1 = $(floor_value) ", message)
      @test MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim=floor_value, cutoff=0.0,
        preserve_diameter=diameter) === rho
    end
  end

  @testset "preserve_diameter = 5 preserves diameter 5 and stops at diameter 6" begin
    d, nsites, chi = 2, 9, 60
    maxdim = 2 * d^4 + 12      # 2 * 16 + 12 = 44
    sites = operator_siteinds(nsites; d=d)
    kept = diameter_probes(nsites, d, 5)
    rho = add(operator_basis_state(sites, fill(1, nsites)),
      0.3 * random_mps(ComplexF64, sites; linkdims=chi); maxdim=chi + 1, cutoff=0.0)
    reference = copy(rho)
    before = operator_expectation_profile(rho, kept; normalize=false)
    MPSToolkit._dmt_bond_truncate!(rho, 5; maxdim=maxdim, cutoff=0.0, preserve_diameter=5)
    after = operator_expectation_profile(rho, kept; normalize=false)
    wide = [k for k in eachindex(kept) if size(kept[k][2], 1) == d^5]
    sweep_error = preservation_error(before, after)
    # Measured on the width-5 probes alone as well, so a large-scale width-1 probe cannot mask a
    # width-5 failure through the shared sup-norm denominator.
    wide_error = preservation_error(before[wide], after[wide])
    @info "preserve_diameter = 5: $(length(kept)) probes error $(sweep_error); " *
          "$(length(wide)) width-5 probes error $(wide_error)"
    @test dim(linkind(rho, 5)) <= maxdim
    @test sweep_error < 1e-11
    @test wide_error < 1e-11

    # At preserve_diameter = 3 the same diameter-5 probes are NOT preserved: this pins that the
    # parameter does something rather than being cosmetic.
    narrow = copy(reference)
    MPSToolkit._dmt_bond_truncate!(narrow, 5; maxdim=maxdim, preserve_diameter=3, cutoff=0.0)
    narrow_error = preservation_error(
      before[wide], operator_expectation_profile(narrow, kept[wide]; normalize=false))
    @info "preserve_diameter = 3 on the same $(length(wide)) width-5 probes: error $(narrow_error)"
    @test narrow_error > 1e-8
  end

  @testset "preserve_diameter = 5 at d = 3" begin
    # The intersection the task exists to prove: wider radius AND higher local dimension at the
    # same time. The floor is 2 * 3^4 + 1 = 163, so the bond has to carry more than ~175 for the
    # truncation to fire at all; six sites is the shortest chain with such a bond (9^3 = 729 at
    # bond 3) that still keeps the width-5 window blocks small enough to measure, because those
    # windows then touch the chain edges.
    d, nsites, chi = 3, 6, 200
    maxdim = 2 * d^4 + 12      # 174, floor 163
    sites = operator_siteinds(nsites; d=d)
    kept = diameter_probes(nsites, d, 5)
    rho = add(operator_basis_state(sites, fill(1, nsites)),
      0.3 * random_mps(ComplexF64, sites; linkdims=chi); maxdim=chi + 1, cutoff=0.0)
    @test dim(linkind(rho, 3)) > maxdim
    reference = copy(rho)
    before = operator_expectation_profile(rho, kept; normalize=false)
    MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim=maxdim, cutoff=0.0, preserve_diameter=5)
    after = operator_expectation_profile(rho, kept; normalize=false)
    wide = [k for k in eachindex(kept) if size(kept[k][2], 1) == d^5]
    sweep_error = preservation_error(before, after)
    wide_error = preservation_error(before[wide], after[wide])
    @info "d = 3, preserve_diameter = 5: $(length(kept)) probes error $(sweep_error); " *
          "$(length(wide)) width-5 probes error $(wide_error)"
    @test dim(linkind(rho, 3)) <= maxdim
    @test sweep_error < 1e-11
    # Measured on the width-5 probes alone as well, so a large-scale width-1 probe cannot mask a
    # width-5 failure through the shared sup-norm denominator.
    @test wide_error < 1e-11

    narrow = copy(reference)
    MPSToolkit._dmt_bond_truncate!(narrow, 3; maxdim=maxdim, cutoff=0.0, preserve_diameter=3)
    narrow_error = preservation_error(
      before[wide], operator_expectation_profile(narrow, kept[wide]; normalize=false))
    @info "d = 3, preserve_diameter = 3 on the same width-5 probes: error $(narrow_error)"
    @test narrow_error > 1e-8
  end

  @testset "the cached environments match the rebuild at radius 2" begin
    # `_dmt_bond_truncate!` stops its environments at `bond - radius` / `bond + 1 + radius`, so
    # the cache invalidation window widens with the radius. `_DMT_VERIFY_ENVS[]` compares every
    # cached environment against a from-scratch rebuild and throws on a mismatch, which is the
    # direct check that the widened window is wide enough. The state starts above `maxdim` on the
    # central bonds so that the truncation -- and therefore the environment fetch -- actually
    # runs; a state that stays under budget short-circuits before the cache is ever touched.
    previous = MPSToolkit._DMT_VERIFY_ENVS[]
    MPSToolkit._DMT_VERIFY_ENVS[] = true
    try
      for (d, nsites, diameter, maxdim) in ((2, 8, 5, 44), (3, 6, 3, 30))
        sites = operator_siteinds(nsites; d=d)
        sz = Matrix{ComplexF64}(Diagonal(collect(Float64, (d - 1):-2:(1 - d))))
        identity_matrix = Matrix{ComplexF64}(I, d, d)
        h = kron(sz, sz) + 0.5 * (kron(sz, identity_matrix) + kron(identity_matrix, sz))
        gate = operator_gate_from_hamiltonian(h, 0.2; d=d)
        rho = add(operator_basis_state(sites, fill(1, nsites)),
          0.3 * random_mps(ComplexF64, sites; linkdims=60); maxdim=61, cutoff=0.0)
        # At least three central bonds start over budget, so the truncation -- and with it the
        # radius-`radius` environment fetch -- runs at bonds far enough from both edges that
        # neither `left_count` nor `right_count` is clamped.
        @test count(b -> dim(linkind(rho, b)) > maxdim, 1:(nsites - 1)) >= 3
        schedule = collect(1:(nsites - 1))
        evo = DMTGateEvolution(gate, 0.2; schedule=schedule, reverse_schedule=reverse(schedule),
          nstep=2, maxdim=maxdim, cutoff=1e-14, preserve_diameter=diameter, normalize=false)
        @test dmt_evolve!(rho, evo) === rho
        @test maximum(dim(linkind(rho, b)) for b in 1:(nsites - 1)) <= maxdim
      end
    finally
      MPSToolkit._DMT_VERIFY_ENVS[] = previous
    end
  end

  @testset "spin-1 Heisenberg melt conserves total S^z end to end" begin
    d, nsites, dt = 3, 8, 0.1
    maxdim = 2 * d^2 + 22      # 40
    sites = operator_siteinds(nsites; d=d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    sx = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sy = ComplexF64[0 -im 0; im 0 -im; 0 im 0] / sqrt(2)
    h = kron(sx, sx) + kron(sy, sy) + kron(sz, sz)
    gate = operator_gate_from_hamiltonian(h, dt; d=d)
    identity_matrix = Matrix{ComplexF64}(I, d, d)
    rho = add(operator_product_state(sites, fill(identity_matrix, nsites)),
      operator_local_sum_state(sites, sz, [j <= nsites ÷ 2 ? -0.25 : 0.25 for j in 1:nsites]);
      maxdim=8, cutoff=0.0)
    profile(state) = real.(operator_expectation_profile(
      state, [(x, sz) for x in 1:nsites]; normalize=false))
    charge(state) = sum(profile(state))
    initial_charge = charge(rho)
    # `normalize = false` reports the literal trace, so a single site carries
    # 0.25 * tr((S^z)^2) * d^(nsites - 1) = 1093.5 here. The conserved quantity must therefore be
    # judged RELATIVE to that scale: an absolute 1e-9 would be a 9e-13 relative demand, below what
    # 10 sweeps of 14 gates hold in double precision. Measured worst relative drift is 1.4e-12,
    # while the same circuit on 6 sites with a budget so large that no truncation fires (maxdim
    # 800, max bond 9^3) drifts 4.9e-13 relative -- MORE than the truncated 6-site run's 2.0e-13.
    # The residue is Trotter-circuit roundoff, not DMT.
    scale = maximum(abs, profile(rho))
    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule=schedule, reverse_schedule=reverse(schedule),
      nstep=1, maxdim=maxdim, cutoff=1e-14, normalize=false)
    worst = 0.0
    for _ in 1:10
      dmt_evolve!(rho, evo)
      worst = max(worst, abs(charge(rho) - initial_charge))
      @test abs(charge(rho) - initial_charge) < 1e-11 * scale
      @test maximum(dim(linkind(rho, b)) for b in 1:(nsites - 1)) <= maxdim
    end
    @info "spin-1 melt: per-site charge scale $(scale), worst drift over 10 sweeps $(worst) " *
          "(relative $(worst / scale))"
  end

  @testset "exact-diagonalization oracle at d = 3" begin
    # With a budget large enough that no truncation fires, DMT evolution must equal dense
    # Liouvillian evolution of the same Trotter circuit.
    d, nsites, dt, nstep = 3, 4, 0.05, 3
    sites = operator_siteinds(nsites; d=d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    sx = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sy = ComplexF64[0 -im 0; im 0 -im; 0 im 0] / sqrt(2)
    identity_matrix = Matrix{ComplexF64}(I, d, d)
    # The field term is load-bearing. Without it, h = SxSx + SzSz leaves every probe below
    # exactly zero except S^z on the site that carries the initial operator, and an oracle that
    # only compares zeros to zeros would pass against almost any kernel.
    h = kron(sx, sx) + kron(sz, sz) + 0.5 * (kron(sz, identity_matrix) + kron(identity_matrix, sz))
    gate = operator_gate_from_hamiltonian(h, dt; d=d)
    rho0 = add(operator_product_state(sites, fill(identity_matrix, nsites)),
      operator_product_state(sites, [j == 2 ? sz : identity_matrix for j in 1:nsites]);
      maxdim=4, cutoff=0.0)

    dense = kron(kron(kron(identity_matrix, identity_matrix), identity_matrix), identity_matrix) +
            kron(kron(kron(identity_matrix, sz), identity_matrix), identity_matrix)
    u2 = exp(-im * dt * Matrix{ComplexF64}(h))
    function embed(u, start)
      left = start == 1 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^(start - 1), d^(start - 1))
      right_sites = nsites - start - 1
      right = right_sites == 0 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^right_sites, d^right_sites)
      return kron(kron(left, u), right)
    end

    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule=schedule, reverse_schedule=reverse(schedule),
      nstep=1, maxdim=400, cutoff=1e-14, normalize=false)
    rho = copy(rho0)
    for _ in 1:nstep
      dmt_evolve!(rho, evo)
      # `dmt_evolve!` runs the forward schedule and then the reverse schedule, so the dense
      # reference must apply the same gates in the same order.
      for start in schedule
        big = embed(u2, start)
        dense = big * dense * big'
      end
      for start in reverse(schedule)
        big = embed(u2, start)
        dense = big * dense * big'
      end
    end
    # 9^2 is the largest bond a 4-site operator-space chain can carry: nothing was truncated,
    # so any disagreement below is a bug rather than a truncation error.
    @test maximum(dim(linkind(rho, b)) for b in 1:(nsites - 1)) == d^(2 * (nsites ÷ 2))

    nontrivial = 0
    for (start, op) in ((1, sz), (2, sz), (3, sz), (4, sz), (1, kron(sz, sz)), (1, kron(sx, sx)),
      (2, kron(sy, sy)), (3, kron(sz, sz)))
      span = size(op, 1) == d ? 1 : 2
      left = start == 1 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^(start - 1), d^(start - 1))
      right_sites = nsites - start - span + 1
      right = right_sites == 0 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^right_sites, d^right_sites)
      padded = kron(kron(left, ComplexF64.(op)), right)
      reference = tr(dense * padded)
      abs(reference) > 1e-4 && (nontrivial += 1)
      @test operator_expectation(rho, op, start; normalize=false) ≈ reference atol = 1e-8
    end
    # Guard against the oracle silently degenerating into a sweep of 0 ≈ 0 comparisons.
    @test nontrivial >= 6
  end
end
