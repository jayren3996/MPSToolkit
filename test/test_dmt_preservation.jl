using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Test

include("dmt_test_helpers.jl")

@testset "DMT preserves local observables exactly" begin
  Random.seed!(20260808)

  @testset "diameter <= 3 preserved to machine precision (d = 2)" begin
    nsites, chi, maxdim = 7, 40, 24
    sites = operator_siteinds(nsites; d=2)
    probes = diameter_probes(nsites, 2, 3)
    # `elt` is the element type each arm is supposed to run in, asserted BEFORE and AFTER the
    # truncation. The "after" assertion is the load-bearing one: `random_mps(sites; ...)` (no
    # element type argument) is `Float64`, but a kernel that coerced its own intermediates to
    # `ComplexF64` returned a MIXED `[Float64, ComplexF64]` chain -- so this arm silently
    # exercised the complex path and the A/B against the complex arm below tested nothing.
    for (label, noise, elt) in (
      ("Hermitian (real coefficients)", random_mps(sites; linkdims=chi), Float64),
      ("non-Hermitian (complex)", random_mps(ComplexF64, sites; linkdims=chi), ComplexF64))
      rho = add(operator_basis_state(sites, fill(1, nsites)), 0.3 * noise;
        maxdim=chi + 1, cutoff=0.0)
      @test mps_eltype(rho) === elt
      before = operator_expectation_profile(rho, probes; normalize=false)
      trace_before = operator_trace(rho)
      MPSToolkit._dmt_bond_truncate!(rho, 4; maxdim=maxdim, cutoff=0.0)
      @test mps_eltype(rho) === elt
      after = operator_expectation_profile(rho, probes; normalize=false)
      @test dim(linkind(rho, 4)) <= maxdim
      @test preservation_error(before, after) < 1e-11
      @test abs(operator_trace(rho) - trace_before) <= 1e-11 * abs(trace_before)
    end
  end

  # One d = 2 XXZ domain-wall melt: returns the worst total-S^z drift over `nsweep` sweeps,
  # relative to the per-site charge scale, and the largest bond dimension reached. The total is
  # ~0 for a symmetric wall, so the drift has to be judged against the per-site scale.
  function melt_drift(nsites, maxdim, nsweep)
    dt = 0.1
    sites = operator_siteinds(nsites; d=2)
    gate = operator_gate_from_hamiltonian(
      spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=1.0, Jz=1.0), dt; d=2)
    rho = add(operator_basis_state(sites, fill(1, nsites)),
      0.25 * pauli_domain_wall_state(sites; kink=nsites ÷ 2); maxdim=8, cutoff=0.0)
    z = ComplexF64[1 0; 0 -1]
    profile(state) = real.(operator_expectation_profile(
      state, [(x, z) for x in 1:nsites]; normalize=false))
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

  @testset "the guarantee survives a full sweep" begin
    nsites, maxdim = 10, 24
    truncating = melt_drift(nsites, maxdim, 15)
    # A conserved total is a WEAK probe of the guarantee: it is exactly what stayed at 1e-6 while
    # individual diameter-3 observables were 26% wrong in the defect this kernel was rebuilt to
    # fix. Asserted here only in company with the control below.
    @test truncating.absolute < 1e-10
    @test truncating.maxlink <= maxdim

    # The control is what gives the testset its name. The same circuit is run with a budget past
    # the 4^5 = 1024 ceiling of a 10-site operator-space chain, so `_dmt_bond_truncate!`
    # short-circuits at every bond and no DMT truncation happens at all -- roundoff in the Trotter
    # circuit and the repeated factorizations is then the only source of drift. If DMT ever does
    # start shedding charge through the truncation, the truncating run pulls away from the
    # untruncated one and this comparison fails, where the absolute bound above would not.
    untruncated = melt_drift(nsites, 4^(nsites ÷ 2), 15)
    @info "d = 2 melt, $(nsites) sites: per-site charge scale $(truncating.scale); truncating " *
          "(max bond $(truncating.maxlink)) relative drift $(truncating.relative); untruncated " *
          "(max bond $(untruncated.maxlink)) relative drift $(untruncated.relative); ratio " *
          "$(truncating.relative / untruncated.relative)"
    # The control never truncated: the bond saturated at the full 4^5 the chain admits, which is
    # the budget it was given, so every bond short-circuited on `dim(link) <= maxdim`.
    @test untruncated.maxlink == 4^(nsites ÷ 2)
    @test truncating.relative < 1e-11
    @test untruncated.relative < 1e-11
    @test truncating.relative <= 3 * untruncated.relative
  end

  @testset "purity monotone: DMT keeps Z >= 1 where plain SVD truncation does not" begin
    nsites, chi, maxdim = 8, 40, 20
    sites = operator_siteinds(nsites; d=2)
    rho = add(operator_basis_state(sites, fill(1, nsites)),
      0.3 * random_mps(sites; linkdims=chi); maxdim=chi + 1, cutoff=0.0)
    # Z = tr(rho) / sqrt(tr(rho^2)); in an orthonormal operator basis tr(rho^2) = <psi|psi>.
    purity_ratio(state) = real(operator_trace(state)) / sqrt(real(inner(state, state)))
    dmt_state = copy(rho)
    MPSToolkit._dmt_bond_truncate!(dmt_state, 4; maxdim=maxdim, cutoff=0.0)
    svd_state = copy(rho)
    orthogonalize!(svd_state, 4)
    u, s, v = svd(svd_state[4], (linkind(svd_state, 3), siteind(svd_state, 4));
      maxdim=maxdim, cutoff=0.0)
    svd_state[4] = u
    svd_state[5] = s * v * svd_state[5]
    @test purity_ratio(dmt_state) >= purity_ratio(svd_state)
  end

  @testset "both sweep directions produce the same state" begin
    # The :R and :L branches differ only in which tensor absorbs the singular values, so the
    # orthogonality centre lands where the next step expects it (bond + 1 for :R, bond for :L).
    # The physical state must be identical.
    nsites, chi, maxdim = 8, 40, 20
    sites = operator_siteinds(nsites; d=2)
    base = add(operator_basis_state(sites, fill(1, nsites)),
      0.3 * random_mps(ComplexF64, sites; linkdims=chi); maxdim=chi + 1, cutoff=0.0)
    for bond in (3, 4)
      right_state = copy(base)
      MPSToolkit._dmt_bond_truncate!(right_state, bond; maxdim=maxdim, cutoff=0.0, direction=:R)
      left_state = copy(base)
      MPSToolkit._dmt_bond_truncate!(left_state, bond; maxdim=maxdim, cutoff=0.0, direction=:L)
      @test abs(inner(right_state, left_state)) / (norm(right_state) * norm(left_state)) ≈ 1.0 atol =
        1e-10

      # ... and the singular values land on the tensor that puts the orthogonality centre where
      # the next step expects it: at bond + 1 for :R (so psi[bond] is a left isometry), at bond
      # for :L (so psi[bond + 1] is a right isometry).
      new_link = linkind(right_state, bond)
      gram = right_state[bond] * dag(prime(right_state[bond], new_link))
      @test matrix(gram, new_link, new_link') ≈ Matrix{ComplexF64}(I, dim(new_link), dim(new_link)) atol =
        1e-10
      new_link = linkind(left_state, bond)
      gram = left_state[bond + 1] * dag(prime(left_state[bond + 1], new_link))
      @test matrix(gram, new_link, new_link') ≈ Matrix{ComplexF64}(I, dim(new_link), dim(new_link)) atol =
        1e-10
    end
  end

  @testset "budget validation names the local dimension" begin
    sites = operator_siteinds(6; d=2)
    rho = random_mps(ComplexF64, sites; linkdims=20)
    # 2 d^2 + 1 = 9 for d = 2. Assert the MESSAGE, not just the type: the whole point of this
    # error is to tell a caller which floor applies at their local dimension and that `maxdim`
    # changed meaning, and a bare `@test_throws ArgumentError` would pass on any of them.
    failure = try
      MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim=8, cutoff=0.0)
      nothing
    catch exception
      exception
    end
    @test failure isa ArgumentError
    message = sprint(showerror, failure)
    @test occursin("local dimension d = 2", message)
    @test occursin("maxdim >= 2 d^(preserve_diameter - 1) + 1 = 9", message)
    @test occursin("total bond dimension, inclusive of the protected block", message)

    @test MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim=9, cutoff=0.0) === rho
    @test_throws ArgumentError DMTOptions(maxdim=30, preserve_diameter=4)   # must be odd
  end

  @testset "the guarantee survives the cutoff, not just cutoff = 0" begin
    # Every other preservation cell in this suite runs at `cutoff = 0.0` while the default
    # everywhere in the package is `1e-12` -- so the guarantee used to be verified only at a
    # setting nobody runs. It did not hold at the others. `cutoff` was applied by `_dmt_refactor`
    # to the REASSEMBLED bond matrix, whose top-k singular directions need not contain the
    # protected row and column spaces, and in operator space the identity dominates `sigma_1`: on
    # a signal `eps` above an infinite-temperature background, every cutoff above `eps` discarded
    # the protected block wholesale. Measured pre-fix, and the failure is total rather than
    # graceful: at d = 2, eps = 1e-6, bond 4, maxdim = 14, `cutoff = 1e-6` collapsed the bond to
    # `chi_out = 1` -- the pure identity -- for a preservation error of 0.98, against 1.5e-7 at
    # `cutoff = 0`. At d = 3 it was chi_out = 1 and 0.76. Even at eps = 1e-3 the guarantee broke
    # by `cutoff = 1e-4` (0.37 at d = 2, 0.21 at d = 3). The cutoff now binds on the complement
    # alone (`_dmt_complement_keep`).
    #
    # The claim under test is cutoff *invariance*, not one fixed bound. The reachable floor is set
    # by the signal-to-background ratio -- an absolute 1e-16 roundoff on `sigma_1` is a relative
    # `1e-16 / eps` on the signal -- so at eps = 1e-6 it is ~1e-7 and no implementation beats it.
    # What must hold is that turning the cutoff up does not move it.
    for d in (2, 3)
      nsites, chi = 8, 40
      maxdim = 2 * d^2 + 6
      sites = operator_siteinds(nsites; d=d)
      # Structured probes only. A dense random probe has an O(1) identity overlap, so its
      # expectation does not scale with `eps` and it would put a constant ~0.8 in
      # `preservation_error`'s denominator -- masking the destruction of a 1e-6 signal by
      # measuring it against a background that survives. Every structured probe carries at least
      # one traceless factor, so the whole profile scales with the signal.
      probes = diameter_probes(nsites, d, 3; nrandom=0)
      covered = [k for (k, (start, op)) in enumerate(probes)
                 if guarantee_covers(start, probe_span(op, d), 4, 1)]
      for eps in (1e-3, 1e-6)
        Random.seed!(9182)
        base = add(operator_basis_state(sites, fill(1, nsites)),
          eps * random_mps(sites; linkdims=chi); maxdim=chi, cutoff=0.0)
        # Without this the cell is vacuous: `_dmt_bond_truncate!` returns untouched when the bond
        # already fits, and every cutoff then trivially "preserves" everything.
        @test dim(linkind(base, 4)) > maxdim
        before = operator_expectation_profile(base, probes)
        # The probes see the signal and nothing else: a structured probe's expectation is
        # proportional to `eps` and below it -- two decades down at d = 2, three at d = 3, where
        # the same random state spreads over 9 basis directions per site instead of 4.
        @test 1e-4 * eps < maximum(abs, before[covered]) < eps

        function truncated_error(cutoff)
          rho = copy(base)
          MPSToolkit._dmt_bond_truncate!(rho, 4; maxdim=maxdim, cutoff=cutoff)
          @test dim(linkind(rho, 4)) <= maxdim
          return preservation_error(before[covered],
                                    operator_expectation_profile(rho, probes)[covered])
        end
        floor_error = truncated_error(0.0)
        # At a resolvable signal the guarantee is machine precision outright; at eps = 1e-6 the
        # conditioning floor above applies and only invariance is assertable.
        eps == 1e-3 && @test floor_error < 1e-11
        for cutoff in (1e-12, 1e-8, 1e-6, 1e-4)
          # A factor of 10 against a pre-fix excursion of 6 to 7 orders of magnitude. `max` keeps
          # a `floor_error` that happens to land at exactly zero from making this unsatisfiable.
          @test truncated_error(cutoff) <= 10 * max(floor_error, 1e-14)
        end
      end
    end
  end
end
