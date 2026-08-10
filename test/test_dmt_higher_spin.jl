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

  # `preservation_error` divides by the sup-norm of the whole `before` profile, so a family with
  # large expectation values (the dense random probes) would otherwise flatter a family with
  # small ones (the structured basis products). Every family is therefore scored against its own
  # denominator, including the widest width on its own -- that is the width the radius is
  # actually on trial for.
  function preservation_families(probes, structured_count, d, diameter)
    spans = [probe_span(op, d) for (_, op) in probes]
    return ["all" => collect(eachindex(probes)),
      "structured products" => collect(1:structured_count),
      "random dense" => collect((structured_count + 1):length(probes)),
      "widest (span $(diameter))" => [k for k in eachindex(probes) if spans[k] == diameter]]
  end

  function assert_preserved(label, probes, structured_count, before, after, d, diameter)
    for (name, index) in preservation_families(probes, structured_count, d, diameter)
      isempty(index) && continue
      value = preservation_error(before[index], after[index])
      @info "$(label): $(name), $(length(index)) probes, preservation error $(value)"
      @test value < 1e-11
    end
    return nothing
  end

  @testset "diameter <= 3 preserved at d = 3 and d = 4" begin
    for d in (3, 4)
      nsites, chi = 6, 60
      maxdim = 2 * d^2 + 12
      sites = operator_siteinds(nsites; d=d)
      # Complete onsite basis at widths 1-2, plus dense random probes per window. A kernel that
      # protected only a few basis directions instead of the whole d^(2 radius) subspace would
      # pass a 4-element capped sweep at d > 2 and fail this one.
      structured = diameter_probes(nsites, d, 3; nrandom=0)
      probes = vcat(structured, random_window_probes(nsites, d, 1:3))
      # `elt` is the element type each arm is supposed to run in, asserted BEFORE and AFTER the
      # truncation; see the same loop in `test_dmt_preservation.jl` for why the "after" assertion
      # is the load-bearing one. `d > 2` is where it matters most: the promoted copy `NDTensors`
      # materializes across a mixed seam scales with the bond, and the bond scales as `d^2`.
      for (label, noise, elt) in (
        ("Hermitian (real coefficients)", random_mps(sites; linkdims=chi), Float64),
        ("non-Hermitian (complex)", random_mps(ComplexF64, sites; linkdims=chi), ComplexF64))
        rho = add(operator_basis_state(sites, fill(1, nsites)), 0.3 * noise;
          maxdim=chi + 1, cutoff=0.0)
        @test mps_eltype(rho) === elt
        @test dim(linkind(rho, 3)) > maxdim          # the truncation really fires
        before = operator_expectation_profile(rho, probes; normalize=false)
        trace_before = operator_trace(rho)
        MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim=maxdim, cutoff=0.0)
        @test mps_eltype(rho) === elt
        after = operator_expectation_profile(rho, probes; normalize=false)
        assert_preserved("d = $(d), $(label)", probes, length(structured), before, after, d, 3)
        @test dim(linkind(rho, 3)) <= maxdim
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

  @testset "preserve_diameter = 5 preserves diameter 5" begin
    for (d, nsites, bond, chi, maxdim) in ((2, 9, 5, 60, 2 * 2^4 + 12),      # 44, floor 33
      (3, 6, 3, 200, 2 * 3^4 + 12))                                          # 174, floor 163
      # The d = 3 row is the intersection the task exists to prove: wider radius AND higher
      # local dimension at once. Its floor is 2 * 3^4 + 1 = 163, so the bond has to carry more
      # than 174 for the truncation to fire; six sites is the shortest chain with such a bond
      # (9^3 = 729 at bond 3) whose width-5 measurement windows still touch the chain edges and
      # so stay small enough to contract.
      sites = operator_siteinds(nsites; d=d)
      structured = diameter_probes(nsites, d, 5; nrandom=0)
      kept = vcat(structured, random_window_probes(nsites, d, 1:5))
      rho = add(operator_basis_state(sites, fill(1, nsites)),
        0.3 * random_mps(ComplexF64, sites; linkdims=chi); maxdim=chi + 1, cutoff=0.0)
      @test dim(linkind(rho, bond)) > maxdim        # the truncation really fires
      reference = copy(rho)
      before = operator_expectation_profile(rho, kept; normalize=false)
      MPSToolkit._dmt_bond_truncate!(rho, bond; maxdim=maxdim, cutoff=0.0, preserve_diameter=5)
      after = operator_expectation_profile(rho, kept; normalize=false)
      assert_preserved("d = $(d), preserve_diameter = 5", kept, length(structured),
        before, after, d, 5)
      @test dim(linkind(rho, bond)) <= maxdim

      # At preserve_diameter = 3 the same diameter-5 probes are NOT preserved: this pins that
      # the parameter does something rather than being cosmetic.
      wide = [k for k in eachindex(kept) if probe_span(kept[k][2], d) == 5]
      narrow = copy(reference)
      MPSToolkit._dmt_bond_truncate!(narrow, bond; maxdim=maxdim, preserve_diameter=3, cutoff=0.0)
      narrow_error = preservation_error(
        before[wide], operator_expectation_profile(narrow, kept[wide]; normalize=false))
      @info "d = $(d), preserve_diameter = 3 on the same $(length(wide)) width-5 probes: " *
            "error $(narrow_error)"
      @test narrow_error > 1e-8
    end
  end

  @testset "the guarantee stops exactly one width later" begin
    # `preserve_diameter = 2 radius + 1` is a sharp edge, not a soft one. At width 2 radius + 2 a
    # window can put more than `radius` sites on BOTH sides of the cut for the first time, and
    # exactly those windows must break -- while the same width on a window that still fits one
    # side inside the radius must stay exact. Random dense probes only: they carry generic
    # weight on every basis product of the window, which is what makes the break unambiguous.
    for (d, nsites, bond, diameter, chi, maxdim) in ((3, 6, 3, 3, 60, 30), (4, 5, 3, 3, 60, 44),
      (2, 9, 5, 5, 60, 44), (3, 6, 3, 5, 200, 174))
      radius = (diameter - 1) ÷ 2
      width = diameter + 1
      sites = operator_siteinds(nsites; d=d)
      probes = random_window_probes(nsites, d, width:width; nrandom=12)
      covered = [k for k in eachindex(probes) if guarantee_covers(probes[k][1], width, bond, radius)]
      uncovered = [k for k in eachindex(probes) if !guarantee_covers(probes[k][1], width, bond, radius)]
      @test !isempty(uncovered)
      rho = add(operator_basis_state(sites, fill(1, nsites)),
        0.3 * random_mps(ComplexF64, sites; linkdims=chi); maxdim=chi + 1, cutoff=0.0)
      @test dim(linkind(rho, bond)) > maxdim
      before = operator_expectation_profile(rho, probes; normalize=false)
      MPSToolkit._dmt_bond_truncate!(rho, bond; maxdim=maxdim, cutoff=0.0,
        preserve_diameter=diameter)
      after = operator_expectation_profile(rho, probes; normalize=false)
      broken = preservation_error(before[uncovered], after[uncovered])
      @info "d = $(d), preserve_diameter = $(diameter): width-$(width) probes beyond the " *
            "guarantee ($(length(uncovered)) of $(length(probes))) error $(broken)"
      @test broken > 1e-4
      if !isempty(covered)
        held = preservation_error(before[covered], after[covered])
        @info "d = $(d), preserve_diameter = $(diameter): width-$(width) probes still inside " *
              "the guarantee ($(length(covered))) error $(held)"
        @test held < 1e-11
      end
    end
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

  # One spin-1 Heisenberg melt: returns the worst relative total-S^z drift over `nsweep` sweeps
  # and the largest bond dimension reached.
  function spin1_melt_drift(nsites, maxdim, nsweep)
    d, dt = 3, 0.1
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
    scale = maximum(abs, profile(rho))
    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule=schedule, reverse_schedule=reverse(schedule),
      nstep=1, maxdim=maxdim, cutoff=1e-14, normalize=false)
    worst = 0.0
    reached = 0
    for _ in 1:nsweep
      dmt_evolve!(rho, evo)
      worst = max(worst, abs(charge(rho) - initial_charge))
      reached = max(reached, maximum(dim(linkind(rho, b)) for b in 1:(nsites - 1)))
      @test maximum(dim(linkind(rho, b)) for b in 1:(nsites - 1)) <= maxdim
    end
    return (relative=worst / scale, absolute=worst, scale=scale, maxlink=reached)
  end

  @testset "spin-1 Heisenberg melt conserves total S^z end to end" begin
    # `normalize = false` reports the literal trace, so a single site carries
    # 0.25 * tr((S^z)^2) * d^(nsites - 1) = 1093.5 at 8 sites. The conserved quantity has to be
    # judged RELATIVE to that scale: an absolute 1e-9 would be a 9e-13 relative demand, below
    # what 10 sweeps of 14 gates hold in double precision.
    main = spin1_melt_drift(8, 2 * 3^2 + 22, 10)     # maxdim 40
    @info "spin-1 melt, 8 sites: per-site charge scale $(main.scale), worst drift over 10 " *
          "sweeps $(main.absolute) (relative $(main.relative)), max bond $(main.maxlink)"
    @test main.relative < 1e-11

    # The residue above is roundoff in the Trotter circuit and the repeated factorizations, NOT
    # charge leaking through the truncation -- and that claim is asserted, not merely asserted in
    # a comment. The same circuit on 6 sites is run twice: once truncating at maxdim 40, once
    # with a budget past the 9^3 = 729 ceiling of a 6-site chain, so `_dmt_bond_truncate!`
    # short-circuits at every bond and no DMT truncation happens at all. Roundoff alone already
    # produces the observed drift, so if DMT ever does start shedding charge the truncating run
    # will pull away from the untruncated one and this comparison fails.
    truncating = spin1_melt_drift(6, 40, 10)
    untruncated = spin1_melt_drift(6, 800, 10)
    @info "spin-1 melt, 6 sites: truncating (max bond $(truncating.maxlink)) relative drift " *
          "$(truncating.relative); untruncated (max bond $(untruncated.maxlink)) relative drift " *
          "$(untruncated.relative); ratio $(truncating.relative / untruncated.relative)"
    @test truncating.maxlink <= 40
    @test untruncated.maxlink == 3^(2 * (6 ÷ 2))     # 729: the budget never bound
    @test truncating.relative < 1e-11
    @test untruncated.relative < 1e-11
    @test truncating.relative <= 3 * untruncated.relative
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
