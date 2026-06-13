using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Test

# Dense exact-diagonalization oracles for the PXP transport infrastructure. Conventions match
# the source: local basis index 1 = |0> (ground), 2 = |1> (excited), kron factor order = chain
# site order (site 1 = most significant bit), normalized Pauli strings sigma/sqrt(2) per site.

const _PXP_ID2 = Matrix{ComplexF64}(I, 2, 2)
const _PXP_X = ComplexF64[0 1; 1 0]
const _PXP_GROUND = ComplexF64[1 0; 0 0]   # |0><0|
const _PXP_EXCITED = ComplexF64[0 0; 0 1]  # |1><1|

# Embed a dense local operator starting at `start` into the full 2^nsites space.
function _embed_term(op, start, nsites)
  span = round(Int, log2(size(op, 1)))
  left = start == 1 ? Matrix{ComplexF64}(I, 1, 1) : foldl(kron, fill(_PXP_ID2, start - 1))
  right_count = nsites - (start + span - 1)
  right = right_count == 0 ? Matrix{ComplexF64}(I, 1, 1) : foldl(kron, fill(_PXP_ID2, right_count))
  return kron(left, ComplexF64.(op), right)
end

function _pxp_dense(nsites; omega=1.0)
  return sum(
    _embed_term(pxp_term_hamiltonian(nsites, j; omega=omega), first(pxp_term_support(nsites, j)), nsites)
    for j in 1:nsites
  )
end

# Diagonal dense P_G = prod_j (1 - n_j n_{j+1}); site 1 = most significant bit.
function _pxp_projector_dense(nsites)
  d = 2^nsites
  diagvals = ones(Float64, d)
  for k in 0:(d - 1)
    bits = reverse(digits(k; base=2, pad=nsites))   # bits[1] = site 1
    if any(bits[i] == 1 && bits[i + 1] == 1 for i in 1:(nsites - 1))
      diagvals[k + 1] = 0.0
    end
  end
  return Matrix{ComplexF64}(Diagonal(diagvals))
end

# Contract an MPO into a dense matrix with big-endian site ordering (matching _embed_term).
function _mpo_dense(op::MPO, sites)
  contracted = op[1]
  for n in 2:length(op)
    contracted *= op[n]
  end
  arr = Array(contracted, prime.(sites)..., sites...)
  n = length(sites)
  perm = vcat(reverse(1:n), reverse((n + 1):(2n)))
  d = prod(dim.(sites))
  return reshape(permutedims(arr, perm), d, d)
end

const _PAULI_NORM = [m / sqrt(2) for m in values(pauli_matrices())]

_pauli_string(labels) = foldl(kron, (_PAULI_NORM[l] for l in labels))

# Vectorized Pauli coefficients of a dense operator: c_alpha = tr(P_alpha' * dense).
function _dense_pauli_amplitude(dense, labels)
  return tr(_pauli_string(labels)' * dense)
end

# Random single-site operator-space product MPS together with its dense operator.
function _random_pauli_product(sites, rng_offset::Int)
  tensors = ITensor[]
  dense = Matrix{ComplexF64}(I, 1, 1)
  for (n, site) in enumerate(sites)
    coeffs = ComplexF64[sin(0.7 * n + rng_offset) + 0.3im * cos(1.3 * n),
                        cos(0.4 * n) - 0.2im * sin(n + rng_offset),
                        0.5 * sin(2.1 * n) + 0.1im,
                        cos(0.9 * n + 2 * rng_offset)]
    tensor = ITensor(site)
    for alpha in 1:4
      tensor[site => alpha] = coeffs[alpha]
    end
    push!(tensors, tensor)
    dense = kron(dense, sum(coeffs[alpha] * _PAULI_NORM[alpha] for alpha in 1:4))
  end
  return MPS(tensors), dense
end

@testset "PXP model helpers" begin
  @testset "pxp_term_support" begin
    @test pxp_term_support(6, 1) == 1:2
    @test pxp_term_support(6, 2) == 1:3
    @test pxp_term_support(6, 5) == 4:6
    @test pxp_term_support(6, 6) == 5:6
    @test pxp_term_support(2, 1) == 1:2
    @test pxp_term_support(2, 2) == 1:2
    @test_throws ArgumentError pxp_term_support(1, 1)
    @test_throws ArgumentError pxp_term_support(4, 0)
    @test_throws ArgumentError pxp_term_support(4, 5)
  end

  @testset "pxp_term_hamiltonian dense blocks" begin
    @test pxp_term_hamiltonian(6, 1) ≈ kron(_PXP_X, _PXP_GROUND)
    @test pxp_term_hamiltonian(6, 6) ≈ kron(_PXP_GROUND, _PXP_X)
    @test pxp_term_hamiltonian(6, 3) ≈ kron(_PXP_GROUND, _PXP_X, _PXP_GROUND)
    @test pxp_term_hamiltonian(6, 3; omega=0.7) ≈ 0.7 * kron(_PXP_GROUND, _PXP_X, _PXP_GROUND)
    @test_throws ArgumentError pxp_term_hamiltonian(1, 1)
    @test_throws ArgumentError pxp_term_hamiltonian(4, 5)
  end

  @testset "embedded terms sum to the open-chain PXP Hamiltonian" begin
    nsites = 5
    h = _pxp_dense(nsites)
    # Independent reference built directly from definitions.
    reference = zeros(ComplexF64, 2^nsites, 2^nsites)
    reference += _embed_term(kron(_PXP_X, _PXP_GROUND), 1, nsites)
    for j in 2:(nsites - 1)
      reference += _embed_term(kron(_PXP_GROUND, _PXP_X, _PXP_GROUND), j - 1, nsites)
    end
    reference += _embed_term(kron(_PXP_GROUND, _PXP_X), nsites - 1, nsites)
    @test h ≈ reference
    @test h ≈ h'
  end

  @testset "pxp_constraint_mpo matches the dense projector" begin
    for nsites in 2:6
      sites = siteinds("S=1/2", nsites)
      projector = pxp_constraint_mpo(sites)
      dense = _mpo_dense(projector, sites)
      @test dense ≈ _pxp_projector_dense(nsites)
      @test dense * dense ≈ dense
    end
    @test_throws ArgumentError pxp_constraint_mpo(siteinds("S=1", 3))
  end

  @testset "every PXP term commutes with the constraint projector" begin
    nsites = 5
    projector = _pxp_projector_dense(nsites)
    for j in 1:nsites
      term = _embed_term(
        pxp_term_hamiltonian(nsites, j),
        first(pxp_term_support(nsites, j)),
        nsites,
      )
      @test norm(term * projector - projector * term) ≈ 0 atol = 1e-12
    end
    h = _pxp_dense(nsites)
    @test norm(h * projector - projector * h) ≈ 0 atol = 1e-12
  end
end

@testset "Pauli vectorization converters" begin
  @testset "pauli_state_from_mpo reproduces tr(P_alpha' O)" begin
    # (a) Hermitian diagonal MPO: the PXP constraint projector on four sites.
    nsites = 4
    phys = siteinds("S=1/2", nsites)
    psites = pauli_siteinds(nsites)
    projector = pxp_constraint_mpo(phys)
    dense = _mpo_dense(projector, phys)
    vectorized = pauli_state_from_mpo(projector, psites)
    @test length(vectorized) == nsites
    @test maxlinkdim(vectorized) <= 2
    for labels in Iterators.product(ntuple(_ -> 1:4, nsites)...)
      amplitude = inner(pauli_basis_state(psites, collect(labels)), vectorized)
      @test amplitude ≈ _dense_pauli_amplitude(dense, labels) atol = 1e-12
    end

    # (b) Non-Hermitian MPO with complex structure.
    os = OpSum()
    os += "Sz", 1
    os += 0.7, "S+", 2
    os += 0.3im, "Sy", 1, "Sx", 2
    phys2 = siteinds("S=1/2", 2)
    psites2 = pauli_siteinds(2)
    mpo = MPO(os, phys2)
    dense2 = _mpo_dense(mpo, phys2)
    vectorized2 = pauli_state_from_mpo(mpo, psites2)
    for labels in Iterators.product(1:4, 1:4)
      amplitude = inner(pauli_basis_state(psites2, collect(labels)), vectorized2)
      @test amplitude ≈ _dense_pauli_amplitude(dense2, labels) atol = 1e-12
    end

    @test_throws ArgumentError pauli_state_from_mpo(projector, pauli_siteinds(3))
  end

  @testset "pauli_superoperator_mpo implements rho -> M rho M'" begin
    # Full dense superoperator oracle on two sites with a non-Hermitian MPO.
    os = OpSum()
    os += "Sz", 1
    os += 0.7, "S+", 2
    os += 0.4, "Sx", 1, "Sz", 2
    phys = siteinds("S=1/2", 2)
    psites = pauli_siteinds(2)
    mpo = MPO(os, phys)
    dense = _mpo_dense(mpo, phys)
    superop = pauli_superoperator_mpo(mpo, psites)
    for in_labels in Iterators.product(1:4, 1:4)
      transformed = apply(superop, pauli_basis_state(psites, collect(in_labels)); cutoff=0.0)
      expected_dense = dense * _pauli_string(in_labels) * dense'
      for out_labels in Iterators.product(1:4, 1:4)
        amplitude = inner(pauli_basis_state(psites, collect(out_labels)), transformed)
        @test amplitude ≈ _dense_pauli_amplitude(expected_dense, out_labels) atol = 1e-10
      end
    end
  end

  @testset "pauli_pxp_constraint_state" begin
    # tr(P_G) counts allowed configurations: Fibonacci F(nsites + 2).
    fibonacci = Dict(1 => 2, 2 => 3, 3 => 5, 4 => 8, 5 => 13)
    for (nsites, count) in fibonacci
      psites = pauli_siteinds(nsites)
      state = pauli_pxp_constraint_state(psites)
      @test maxlinkdim(state) <= 2
      @test pauli_trace(state) ≈ count atol = 1e-10
    end
  end

  @testset "pauli_pxp_constraint_projector acts as P rho P" begin
    nsites = 3
    psites = pauli_siteinds(nsites)
    superop = pauli_pxp_constraint_projector(psites)
    dense_projector = _pxp_projector_dense(nsites)

    rho, dense_rho = _random_pauli_product(psites, 1)
    projected = apply(superop, rho; cutoff=0.0)
    expected_dense = dense_projector * dense_rho * dense_projector
    for labels in Iterators.product(ntuple(_ -> 1:4, nsites)...)
      amplitude = inner(pauli_basis_state(psites, collect(labels)), projected)
      @test amplitude ≈ _dense_pauli_amplitude(expected_dense, labels) atol = 1e-10
    end

    # Idempotent: applying twice equals applying once.
    twice = apply(superop, projected; cutoff=0.0)
    @test inner(twice, twice) + inner(projected, projected) - 2 * real(inner(twice, projected)) ≈ 0 atol = 1e-10

    # Fixes the vectorized constraint projector itself.
    constraint_state = pauli_pxp_constraint_state(psites)
    fixed = apply(superop, constraint_state; cutoff=0.0)
    difference = inner(fixed, fixed) + inner(constraint_state, constraint_state) - 2 * real(inner(fixed, constraint_state))
    @test difference ≈ 0 atol = 1e-10

    # Annihilates an operator supported entirely on blockade-violating configurations.
    blocked_labels = [(1, 1), (1, 4), (4, 1), (4, 4)]
    blocked_signs = Dict(1 => 1.0, 4 => -1.0)
    psites2 = pauli_siteinds(2)
    violating = sum(
      0.5 * blocked_signs[a] * blocked_signs[b] * pauli_basis_state(psites2, [a, b])
      for (a, b) in blocked_labels
    )  # |11><11| = (I - Z)/2 ⊗ (I - Z)/2 in normalized Pauli coordinates
    superop2 = pauli_pxp_constraint_projector(psites2)
    annihilated = apply(superop2, violating; cutoff=0.0)
    @test norm(annihilated) ≈ 0 atol = 1e-10
  end
end

@testset "operator-space expectations" begin
  nsites = 4
  psites = pauli_siteinds(nsites)
  rho, dense_rho = _random_pauli_product(psites, 3)

  # Make the operator Hermitian so normalized expectations are real diagnostics, but keep the
  # state generic (non-Hermitian) to exercise the full complex path.
  _random_hermitian(span, seed) = begin
    raw = ComplexF64[
      sin(0.31 * seed * (i + 2j)) + 1im * cos(0.17 * seed * (2i + j)) for i in 1:(2^span), j in 1:(2^span)
    ]
    raw + raw'
  end

  @testset "pauli_trace matches the dense trace" begin
    @test pauli_trace(rho) ≈ tr(dense_rho) atol = 1e-10
  end

  @testset "pauli_expectation matches dense traces" begin
    cases = [(1, _random_hermitian(2, 1)), (2, _random_hermitian(2, 2)), (3, _random_hermitian(1, 3)),
             (1, _random_hermitian(3, 4)), (2, _random_hermitian(3, 5)), (4, _random_hermitian(1, 6))]
    for (start, op) in cases
      embedded = _embed_term(op, start, nsites)
      expected = tr(dense_rho * embedded) / tr(dense_rho)
      @test pauli_expectation(rho, op, start) ≈ expected atol = 1e-10
      @test pauli_expectation(rho, op, start; normalize=false) ≈ tr(dense_rho * embedded) atol = 1e-10
    end
    @test_throws ArgumentError pauli_expectation(rho, _random_hermitian(3, 7), 3)
    @test_throws ArgumentError pauli_expectation(rho, ones(ComplexF64, 3, 3), 1)
  end

  @testset "pauli_expectation_profile equals per-term expectations" begin
    terms = [(first(pxp_term_support(nsites, j)), pxp_term_hamiltonian(nsites, j)) for j in 1:nsites]
    profile = pauli_expectation_profile(rho, terms)
    for (k, (start, op)) in enumerate(terms)
      @test profile[k] ≈ pauli_expectation(rho, op, start) atol = 1e-12
    end
    # Unsorted input keeps input ordering of outputs.
    shuffled = reverse(terms)
    reversed_profile = pauli_expectation_profile(rho, shuffled)
    @test reversed_profile ≈ reverse(profile) atol = 1e-12
    # Unnormalized variant scales by the trace.
    raw_profile = pauli_expectation_profile(rho, terms; normalize=false)
    @test raw_profile ≈ profile .* tr(dense_rho) atol = 1e-9
  end

  @testset "Hermitian state gives real energies" begin
    psites2 = pauli_siteinds(2)
    hermitian_rho = sum(
      0.5 * coeff * pauli_basis_state(psites2, collect(labels)) for
      (coeff, labels) in zip((1.0, 0.3, 0.3, 0.7), ((1, 1), (1, 4), (4, 1), (4, 4)))
    )
    value = pauli_expectation(hermitian_rho, Matrix(pxp_term_hamiltonian(2, 1)), 1)
    @test abs(imag(value)) < 1e-12
  end

  @testset "error and edge paths" begin
    @test isempty(pauli_expectation_profile(rho, Tuple{Int,Matrix{ComplexF64}}[]))
    @test_throws ArgumentError pauli_expectation_profile(rho, [(0, _random_hermitian(1, 1))])
    @test_throws ArgumentError pauli_expectation_profile(rho, [(nsites, _random_hermitian(2, 1))])
    @test_throws ArgumentError pauli_expectation_profile(rho, [(1, ones(ComplexF64, 2, 3))])
    # numerics-1: a numerically-negligible (not exactly zero) trace is rejected under
    # normalize=true, but allowed under normalize=false.
    leaky = add(
      pauli_basis_state(psites, [4, 1, 1, 1]),
      pauli_basis_state(psites, [1, 1, 1, 1]; coefficient=1e-13);
      maxdim=8,
      cutoff=0.0,
    )
    @test_throws ArgumentError pauli_expectation_profile(leaky, [(1, _random_hermitian(1, 1))])
    @test length(pauli_expectation_profile(leaky, [(1, _random_hermitian(1, 1))]; normalize=false)) == 1
    # test-7: non-Pauli (non-dimension-4) sites are rejected by trace and profile.
    bad = MPS(siteinds("S=1", 3), n -> "Up")
    @test_throws ArgumentError pauli_trace(bad)
    @test_throws ArgumentError pauli_expectation_profile(bad, [(1, ComplexF64[1.0 0.0; 0.0 -1.0])])
  end
end

# Dense K = sum_j w_j h_j for the PXP term list.
function _dense_weighted_pxp(nsites, weights; omega=1.0)
  return sum(
    weights[j] * _embed_term(
      pxp_term_hamiltonian(nsites, j; omega=omega),
      first(pxp_term_support(nsites, j)),
      nsites,
    ) for j in 1:nsites
  )
end

function _dense_pxp_profile(rho_dense, nsites)
  return [
    real(
      tr(rho_dense * _embed_term(
        pxp_term_hamiltonian(nsites, j),
        first(pxp_term_support(nsites, j)),
        nsites,
      )) / tr(rho_dense),
    ) for j in 1:nsites
  ]
end

_pxp_terms(nsites; omega=1.0) =
  [(first(pxp_term_support(nsites, j)), pxp_term_hamiltonian(nsites, j; omega=omega)) for j in 1:nsites]

@testset "imaginary-time thermal preparation" begin
  @testset "pauli_gate_from_imaginary_time" begin
    h = Matrix(pxp_term_hamiltonian(4, 2))
    dbeta = 0.37
    @test pauli_gate_from_imaginary_time(h, dbeta) ≈ pauli_gate(exp(-(dbeta / 2) * h)) atol = 1e-12
    @test pauli_gate_from_imaginary_time(h, 0.0) ≈ Matrix{ComplexF64}(I, 4^3, 4^3) atol = 1e-12
    @test_throws ArgumentError pauli_gate_from_imaginary_time(ones(2, 3), 0.1)
    @test_throws ArgumentError pauli_gate_from_imaginary_time([0.0 1.0; 0.0 0.0], 0.1)
  end

  nsites = 4
  beta = 0.6
  psites = pauli_siteinds(nsites)
  terms = _pxp_terms(nsites)
  projector_dense = _pxp_projector_dense(nsites)

  function _mps_profile(state)
    return real.(pauli_expectation_profile(state, terms))
  end

  @testset "uniform-beta constrained Gibbs state matches ED" begin
    weights = fill(beta, nsites)
    k_dense = _dense_weighted_pxp(nsites, weights)
    rho_dense = exp(-Matrix(k_dense) / 2) * projector_dense * exp(-Matrix(k_dense) / 2)
    expected = _dense_pxp_profile(rho_dense, nsites)

    errors = Float64[]
    for nsteps in (3, 20)
      state = pauli_gibbs_state(
        psites,
        terms,
        weights;
        nsteps=nsteps,
        maxdim=256,
        cutoff=0.0,
        initial_state=pauli_pxp_constraint_state(psites),
      )
      push!(errors, maximum(abs.(_mps_profile(state) - expected)))
    end
    @test errors[end] < 1e-3
    @test errors[1] > errors[end]   # Trotter convergence with nsteps
  end

  @testset "domain-wall weights match the factorized dense state" begin
    weights = [beta, 0.0, 0.0, -beta]
    k_dense = _dense_weighted_pxp(nsites, weights)
    rho_dense = exp(-Matrix(k_dense) / 2) * projector_dense * exp(-Matrix(k_dense) / 2)
    expected = _dense_pxp_profile(rho_dense, nsites)

    state = pauli_gibbs_state(
      psites,
      terms,
      weights;
      nsteps=16,
      maxdim=256,
      cutoff=0.0,
      initial_state=pauli_pxp_constraint_state(psites),
    )
    @test maximum(abs.(_mps_profile(state) - expected)) < 1e-3
  end

  @testset "identity default initial state" begin
    weights = [beta, 0.0, 0.0, -beta]
    k_dense = _dense_weighted_pxp(nsites, weights)
    rho_dense = exp(-Matrix(k_dense))   # e^{-K/2} * I * e^{-K/2}
    expected = _dense_pxp_profile(rho_dense, nsites)

    state = pauli_gibbs_state(psites, terms, weights; nsteps=16, maxdim=256, cutoff=0.0)
    @test maximum(abs.(_mps_profile(state) - expected)) < 1e-3
  end

  @testset "argument validation" begin
    @test_throws ArgumentError pauli_gibbs_state(psites, terms, fill(beta, nsites - 1))
    @test_throws ArgumentError pauli_gibbs_state(psites, terms, fill(beta, nsites); nsteps=0)
    # nstep is accepted as an alias for nsteps (consistency with the evolution drivers).
    @test_throws ArgumentError pauli_gibbs_state(psites, terms, fill(beta, nsites); nstep=0)
    @test pauli_gibbs_state(psites, terms, fill(beta, nsites); nstep=2) isa MPS
  end
end

@testset "constrained DMT evolution" begin
  nsites = 6
  beta = 0.4
  dt = 0.05
  nstep = 5     # one sweep advances 2*dt -> t_total = 0.5
  psites = pauli_siteinds(nsites)
  terms = _pxp_terms(nsites)
  # Wall between sites 3 and 4: terms fully inside a half keep +/- beta, straddlers get 0.
  weights = [beta, beta, 0.0, 0.0, -beta, -beta]

  rho = pauli_gibbs_state(
    psites,
    terms,
    weights;
    nsteps=12,
    maxdim=256,
    cutoff=0.0,
    initial_state=pauli_pxp_constraint_state(psites),
  )

  # Dense reference: rho(t) = e^{-iHt} rho_0 e^{+iHt}.
  k_dense = _dense_weighted_pxp(nsites, weights)
  rho_dense = exp(-Matrix(k_dense) / 2) * _pxp_projector_dense(nsites) * exp(-Matrix(k_dense) / 2)
  h_dense = _pxp_dense(nsites)
  t_total = 2 * dt * nstep
  u_dense = exp(-1im * Matrix(h_dense) * t_total)
  rho_dense_t = u_dense * rho_dense * u_dense'
  expected_profile = _dense_pxp_profile(rho_dense_t, nsites)

  gates = [pauli_gate_from_hamiltonian(h, dt) for (_, h) in terms]
  schedule = [start for (start, _) in terms]
  evo = DMTGateEvolution(
    gates,
    dt;
    schedule=schedule,
    reverse_schedule=reverse(schedule),
    nstep=nstep,
    maxdim=64,
    cutoff=1e-12,
    gate_maxdim=256,
    connector_buffer=4,
  )
  projector = pauli_pxp_constraint_projector(psites)

  evolved = copy(rho)
  constrained_dmt_evolve!(evolved, evo, projector; project_every=2)
  profile = real.(pauli_expectation_profile(evolved, terms))
  @test maximum(abs.(profile - expected_profile)) < 5e-3

  # The state stays in the constrained sector: P rho P = rho up to numerical error.
  projected = apply(projector, evolved; maxdim=256, cutoff=0.0)
  overlap = inner(projected, projected) + inner(evolved, evolved) - 2 * real(inner(projected, evolved))
  @test overlap ≈ 0 atol = 1e-8

  # Total energy is conserved by the evolution (and matches the dense value).
  expected_total = real(tr(rho_dense * h_dense) / tr(rho_dense))
  @test sum(profile) ≈ expected_total atol = 5e-3

  @test_throws ArgumentError constrained_dmt_evolve!(evolved, evo, projector; project_every=0)
  @test_throws ArgumentError constrained_dmt_evolve!(evolved, evo, pauli_pxp_constraint_projector(pauli_siteinds(4)))

  # test-5: a non-default projector budget still matches dense ED at exact bond dimension.
  evolved2 = copy(rho)
  constrained_dmt_evolve!(evolved2, evo, projector; project_every=2, projector_maxdim=256, projector_cutoff=1e-14)
  @test maximum(abs.(real.(pauli_expectation_profile(evolved2, terms)) - expected_profile)) < 5e-3
  # project_every is irrelevant at exact bond dimension (no truncation leakage to remove).
  every1 = copy(rho)
  constrained_dmt_evolve!(every1, evo, projector; project_every=1)
  everyN = copy(rho)
  constrained_dmt_evolve!(everyN, evo, projector; project_every=nstep)
  @test maximum(abs.(real.(pauli_expectation_profile(every1, terms)) - real.(pauli_expectation_profile(everyN, terms)))) < 1e-6
  # projector_cutoff is validated.
  @test_throws ArgumentError constrained_dmt_evolve!(copy(rho), evo, projector; projector_cutoff=-1e-12)
end

@testset "energy-correlator protocol matches dense ED (normalize=false)" begin
  # End-to-end check of the pxp_energy_correlator.jl protocol: a traceless, sector-projected
  # energy density O(0) = P_G h_center P_G, Heisenberg-evolved with normalize=false, measured
  # as the unnormalized profile C(x) = tr(O h_x). At N=6 with maxdim = 4^3 = 64 the
  # operator-space evolution is exact (no truncation), so it matches dense ED.
  nsites = 6
  center = nsites ÷ 2
  dt = 0.05
  nstep = 4                       # t_total = 2*dt*nstep = 0.4
  psites = pauli_siteinds(nsites)
  terms = _pxp_terms(nsites)
  projector = pauli_pxp_constraint_projector(psites)

  # O(0): vectorize the center PXP term, sector-project, and HS-normalize (as the script does).
  phys = siteinds("S=1/2", nsites)
  os = OpSum()
  os += "ProjUp", center - 1, "X", center, "ProjUp", center + 1
  center_mpo = MPO(os, phys)
  O = pauli_state_from_mpo(center_mpo, psites)
  O = apply(projector, O; maxdim=256, cutoff=0.0)
  normalize!(O)

  # Dense O(0) from the SAME operator: o0 = P_G h_center P_G / ||.||_HS (HS norm = the
  # vectorized MPS norm that normalize! divides by).
  pg = _pxp_projector_dense(nsites)
  o0_raw = pg * _mpo_dense(center_mpo, phys) * pg
  o0 = o0_raw / sqrt(real(tr(o0_raw' * o0_raw)))
  h_dense = [_embed_term(h, start, nsites) for (start, h) in terms]

  # O is traceless => normalize=true must be rejected (numerics-1); normalize=false gives
  # tr(O h_x) exactly, matching the dense correlator with no extra factors. Tight at t=0.
  @test_throws ArgumentError pauli_expectation_profile(O, terms)
  profile0 = real.(pauli_expectation_profile(O, terms; normalize=false))
  @test profile0 ≈ [real(tr(o0 * hx)) for hx in h_dense] atol = 1e-8

  # Heisenberg-evolve at exact bond dimension with normalize=false.
  gates = [pauli_gate_from_hamiltonian(h, dt) for (_, h) in terms]
  schedule = [start for (start, _) in terms]
  evo = DMTGateEvolution(
    gates,
    dt;
    schedule=schedule,
    reverse_schedule=reverse(schedule),
    nstep=nstep,
    maxdim=64,
    cutoff=0.0,
    gate_maxdim=256,
    connector_buffer=8,
  )
  constrained_dmt_evolve!(O, evo, projector; project_every=1, normalize=false)
  profile_t = real.(pauli_expectation_profile(O, terms; normalize=false))

  # Dense reference uses the exact propagator; the MPS uses a 2nd-order Trotter sweep, so the
  # site-by-site match is at Trotter tolerance (a sign/normalization/sector bug would be O(1)).
  u = exp(-1im * Matrix(_pxp_dense(nsites)) * (2 * dt * nstep))
  o_t = u * o0 * u'
  @test maximum(abs.(profile_t - [real(tr(o_t * hx)) for hx in h_dense])) < 5e-3

  # Conserved total: sum_x C(x,t) = tr(H O(t)) = tr(H O(0)) (Trotter tolerance).
  @test sum(profile_t) ≈ sum(real(tr(o0 * hx)) for hx in h_dense) atol = 5e-3
end

@testset "constrained evolution with normalize=false preserves absolute scales" begin
  nsites = 6
  psites = pauli_siteinds(nsites)
  terms = _pxp_terms(nsites)
  weights = [0.4, 0.4, 0.0, 0.0, -0.4, -0.4]
  rho = pauli_gibbs_state(
    psites,
    terms,
    weights;
    nsteps=8,
    maxdim=256,
    cutoff=0.0,
    initial_state=pauli_pxp_constraint_state(psites),
  )
  gates = [pauli_gate_from_hamiltonian(h, 0.05) for (_, h) in terms]
  schedule = [start for (start, _) in terms]
  evo = DMTGateEvolution(
    gates,
    0.05;
    schedule=schedule,
    reverse_schedule=reverse(schedule),
    nstep=4,
    maxdim=64,
    cutoff=1e-12,
    gate_maxdim=256,
    connector_buffer=4,
  )
  projector = pauli_pxp_constraint_projector(psites)

  # Unitary evolution preserves tr(rho); at N=6 with maxdim=64 the operator space is exact,
  # the state stays in-sector, and with normalize=false nothing rescales the operator, so
  # the absolute trace must be conserved through sweeps + projections.
  trace0 = pauli_trace(rho)
  evolved = copy(rho)
  constrained_dmt_evolve!(evolved, evo, projector; project_every=2, normalize=false)
  # Conservation up to the cutoff-level truncation of the projector applications (the
  # normalize=true path would instead rescale the trace by O(1) over a long run).
  @test pauli_trace(evolved) ≈ trace0 atol = 1e-6 * abs(trace0)

  # The default normalized path rescales the state but leaves ratio observables unchanged:
  # both paths must produce the same energy profile.
  normalized = copy(rho)
  constrained_dmt_evolve!(normalized, evo, projector; project_every=2)
  profile_free = real.(pauli_expectation_profile(evolved, terms))
  profile_normalized = real.(pauli_expectation_profile(normalized, terms))
  @test maximum(abs.(profile_free - profile_normalized)) < 1e-8

  # dmt_evolve! honors normalize=false as well.
  unnormalized = copy(rho)
  dmt_evolve!(unnormalized, evo; normalize=false)
  @test pauli_trace(unnormalized) ≈ trace0 atol = 1e-6 * abs(trace0)
end
