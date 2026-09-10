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

function _pxp_mu_term(nsites, center; omega=1.0, mu=0.0)
  support = pxp_term_support(nsites, center)
  center_offset = center - first(support) + 1
  factors = [offset == center_offset ? mu * _PXP_EXCITED : _PXP_ID2
    for offset in 1:length(support)]
  return pxp_term_hamiltonian(nsites, center; omega=omega) + foldl(kron, factors)
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

  @testset "controlled MPO gates equal dense PXP gates" begin
    for nsites in 2:8
      sites = pauli_siteinds(nsites)
      for center in unique((1, cld(nsites, 2), nsites)), duration in (-0.07, 0.11)
        gate = pauli_pxp_controlled_gate(sites, center, duration; omega=0.8)
        support = pxp_term_support(nsites, center)
        local_sites = collect(sites[support])
        dense = pauli_gate_from_hamiltonian(
          pxp_term_hamiltonian(nsites, center; omega=0.8), duration)
        @test _mpo_dense(gate.mpo, local_sites) ≈ dense atol=1e-12
        @test maxlinkdim(gate.mpo) <= 4
        @test mapreduce(eltype, promote_type, ITensorMPS.data(gate.mpo)) === Float64
        @test gate.start == first(support)
        @test gate.span == length(support)
      end
    end
  end


  @testset "controlled PXP plus chemical-potential gates" begin
    for nsites in 2:8
      sites = pauli_siteinds(nsites)
      for center in unique((1, cld(nsites, 2), nsites)), duration in (-0.07, 0.11), mu in (-0.4, 0.3)
        gate = pauli_pxp_controlled_gate(sites, center, duration; omega=0.8, mu=mu)
        support = pxp_term_support(nsites, center)
        dense = pauli_gate_from_hamiltonian(
          _pxp_mu_term(nsites, center; omega=0.8, mu=mu), duration)
        @test _mpo_dense(gate.mpo, collect(sites[support])) ≈ dense atol=2e-12
        @test maxlinkdim(gate.mpo) <= 4
        @test gate.mu == mu
      end
    end
    sites = pauli_siteinds(5)
    onsite_only = pxp_s4_dmt_plan(sites, 0.1; omega=0.0, mu=-0.7, maxdim=16)
    @test all(gate.omega == 0.0 && gate.mu == -0.7 for gate in onsite_only.evolution.gate)
    @test occursin("mu=-0.7", onsite_only.signature)

    nsites = 4
    omega, mu, total_time = 0.8, -0.35, 0.4
    h = sum(_embed_term(_pxp_mu_term(nsites, center; omega=omega, mu=mu),
      first(pxp_term_support(nsites, center)), nsites) for center in 1:nsites)
    for (builder, threshold) in ((pxp_s2_dmt_plan, 3.5), (pxp_s4_dmt_plan, 10.0))
      errors = Float64[]
      for nsteps in (builder === pxp_s2_dmt_plan ? (2, 4, 8) : (1, 2, 4))
        plan = builder(pauli_siteinds(nsites), total_time / nsteps;
          omega=omega, mu=mu, maxdim=16)
        step = Matrix{ComplexF64}(I, 2^nsites, 2^nsites)
        for entry in plan.entries
          local_gate = exp(-im * entry.duration *
            _pxp_mu_term(nsites, entry.center; omega=omega, mu=mu))
          step = _embed_term(local_gate, entry.start, nsites) * step
        end
        push!(errors, norm(step^nsteps - exp(-im * total_time * h)))
      end
      @test errors[1] / errors[2] > threshold
      @test errors[2] / errors[3] > threshold
    end
    left_odd = _embed_term(_pxp_mu_term(nsites, 1; omega=omega, mu=mu), 1, nsites)
    right_odd = _embed_term(_pxp_mu_term(nsites, 3; omega=omega, mu=mu), 2, nsites)
    @test norm(left_odd * right_odd - right_odd * left_odd) <= 1e-13
    @test_throws ArgumentError pauli_pxp_controlled_gate(sites, 2, Inf; mu=mu)
  end

  @testset "S2 parity plan metadata and order" begin
    for nsites in 2:8
      plan = pxp_s2_dmt_plan(pauli_siteinds(nsites), 0.2; omega=0.8,
        maxdim=12, normalize=false)
      @test plan.scheme == :S2
      @test plan.tau == 0.2
      @test plan.evolution.gate_backend == :controlled
      @test isempty(plan.evolution.reverse_schedule)
      @test [entry.id for entry in plan.entries] == collect(eachindex(plan.entries))
      @test all(entry.direction == :R for entry in plan.entries)
      @test pxp_dmt_attempted_updates(plan) == 3 * (nsites - 1)
      @test [entry.center for entry in plan.entries if entry.layer == 1] ==
        collect(1:2:nsites)
      @test [entry.center for entry in plan.entries if entry.layer == 2] ==
        collect(2:2:nsites)
      @test all(entry.duration == 0.1 for entry in plan.entries if entry.layer != 2)
      @test all(entry.duration == 0.2 for entry in plan.entries if entry.layer == 2)
    end

    nsites = 4
    total_time = 0.4
    exact = exp(-1im * total_time * _pxp_dense(nsites; omega=0.8))
    errors = Float64[]
    for nsteps in (2, 4, 8)
      tau = total_time / nsteps
      plan = pxp_s2_dmt_plan(pauli_siteinds(nsites), tau; omega=0.8, maxdim=16)
      step = Matrix{ComplexF64}(I, 2^nsites, 2^nsites)
      for entry in plan.entries
        local_gate = exp(-1im * entry.duration *
          pxp_term_hamiltonian(nsites, entry.center; omega=0.8))
        step = _embed_term(local_gate, entry.start, nsites) * step
      end
      push!(errors, norm(step^nsteps - exact))
    end
    @test errors[1] > errors[2] > errors[3]
    @test errors[1] / errors[2] > 3.5
    @test errors[2] / errors[3] > 3.5
    @test_throws ArgumentError pxp_s2_dmt_plan(pauli_siteinds(4), 0.0)
  end

  @testset "S4 Yoshida plan metadata and order" begin
    for nsites in 2:8
      plan = pxp_s4_dmt_plan(pauli_siteinds(nsites), 0.2; omega=0.8,
        maxdim=12, normalize=false)
      @test plan.scheme == :S4
      @test isempty(plan.evolution.reverse_schedule)
      @test pxp_dmt_attempted_updates(plan) == 7 * (nsites - 1)
      @test all(entry.direction == :R for entry in plan.entries)
      @test any(entry.duration < 0 for entry in plan.entries)
      @test [entry.center for entry in plan.entries if entry.layer == 1] ==
        collect(1:2:nsites)
      @test [entry.center for entry in plan.entries if entry.layer == 2] ==
        collect(2:2:nsites)
      odd_duration = sum(entry.duration for entry in plan.entries if entry.center == 1)
      even_center = 2
      even_duration = sum(entry.duration for entry in plan.entries
        if entry.center == even_center)
      @test odd_duration ≈ plan.tau atol=1e-14
      @test even_duration ≈ plan.tau atol=1e-14
    end

    nsites = 4
    total_time = 0.4
    exact = exp(-1im * total_time * _pxp_dense(nsites; omega=0.8))
    errors = Float64[]
    for nsteps in (1, 2, 4)
      tau = total_time / nsteps
      plan = pxp_s4_dmt_plan(pauli_siteinds(nsites), tau; omega=0.8, maxdim=16)
      step = Matrix{ComplexF64}(I, 2^nsites, 2^nsites)
      for entry in plan.entries
        local_gate = exp(-1im * entry.duration *
          pxp_term_hamiltonian(nsites, entry.center; omega=0.8))
        step = _embed_term(local_gate, entry.start, nsites) * step
      end
      push!(errors, norm(step^nsteps - exact))
    end
    @test errors[1] > errors[2] > errors[3]
    @test errors[1] / errors[2] > 10
    @test errors[2] / errors[3] > 10
  end

  @testset "parity-layer MPOs equal sequential controlled gates" begin
    for nsites in 2:6, parity in (:odd, :even)
      duration = isodd(nsites) ? -0.07 : 0.11
      sites = pauli_siteinds(nsites)
      layer = pauli_pxp_parity_layer(sites, parity, duration; omega=0.8)
      @test maxlinkdim(layer.mpo) <= 4
      @test mapreduce(eltype, promote_type, ITensorMPS.data(layer.mpo)) === Float64
      state, _ = _random_pauli_product(sites, nsites + (parity === :odd ? 1 : 2))
      via_layer = apply(layer.mpo, state; cutoff=0.0)
      via_gates = copy(state)
      centers = parity === :odd ? (1:2:nsites) : (2:2:nsites)
      for center in centers
        gate = pauli_pxp_controlled_gate(sites, center, duration; omega=0.8)
        via_gates = apply(MPSToolkit._embed_local_mpo(gate.mpo, sites, gate.start),
          via_gates; cutoff=0.0)
      end
      @test abs(inner(via_layer, via_gates)) /
        (norm(via_layer) * norm(via_gates)) ≈ 1.0 atol=1e-11
    end

    for scheme in (:S2, :S4)
      plan = pxp_layer_dmt_plan(pauli_siteinds(6), 0.1; scheme=scheme, maxdim=12,
        normalize=false)
      expected_layers = scheme === :S2 ? 3 : 7
      @test length(plan.layers) == expected_layers
      @test pxp_dmt_attempted_updates(plan) == expected_layers * 5
    end


    @test_throws ArgumentError pxp_layer_dmt_plan(pauli_siteinds(4), 0.1; nstep=0)
    @test_throws ArgumentError pxp_layer_dmt_plan(pauli_siteinds(4), 0.1; nstep=-1)
    sites = pauli_siteinds(4)
    right = pxp_layer_dmt_plan(sites, 0.1; maxdim=16, direction=:R)
    left = pxp_layer_dmt_plan(sites, 0.1; maxdim=16, direction=:L)
    @test right.direction == :R
    @test left.direction == :L
    @test right.signature != left.signature
    @test right.signature != pxp_layer_dmt_plan(sites, 0.1;
      maxdim=12, direction=:R).signature

    state, _ = _random_pauli_product(sites, 19)
    via_right = copy(state)
    via_left = copy(state)
    dmt_evolve!(via_right, right; normalize=false)
    dmt_evolve!(via_left, left; normalize=false)
    @test norm(via_right - via_left) <= 1e-9 * norm(via_right)

    @test_throws ArgumentError pxp_layer_dmt_plan(sites, 0.1; maxdim=1, normalize=true)
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
      @test mapreduce(eltype, promote_type, ITensorMPS.data(state)) === Float64
      legacy = pauli_state_from_mpo(
        pxp_constraint_mpo(MPSToolkit._pxp_physical_sites(nsites)), psites)
      @test norm(state) ≈ norm(legacy) atol=1e-12
      @test abs(inner(state, legacy)) / (norm(state) * norm(legacy)) ≈ 1.0 atol=1e-12
    end
  end

  @testset "pauli_pxp_constraint_projector acts as P rho P" begin
    nsites = 3
    psites = pauli_siteinds(nsites)
    superop = pauli_pxp_constraint_projector(psites)
    dense_projector = _pxp_projector_dense(nsites)
    @test mapreduce(eltype, promote_type, ITensorMPS.data(superop)) === Float64
    legacy_superop = pauli_superoperator_mpo(
      pxp_constraint_mpo(MPSToolkit._pxp_physical_sites(nsites)), psites)
    @test _mpo_dense(superop, psites) ≈ _mpo_dense(legacy_superop, psites) atol=1e-12

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
    # Regression: start beyond the chain length used to reach siteind(rho, start) while
    # computing span, raising a raw BoundsError instead of this ArgumentError.
    @test_throws ArgumentError pauli_expectation_profile(rho, [(length(rho) + 1, _random_hermitian(1, 1))])
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

  @testset "exact gate backends track at finite chi" begin
    product_evo = DMTGateEvolution(gates, dt; schedule=schedule,
      reverse_schedule=reverse(schedule), nstep=2, maxdim=12, cutoff=1e-12,
      gate_maxdim=0, gate_backend=:product, normalize=false)
    qr_evo = DMTGateEvolution(gates, dt; schedule=schedule,
      reverse_schedule=reverse(schedule), nstep=2, maxdim=12, cutoff=1e-12,
      gate_maxdim=0, gate_backend=:qr, normalize=false)
    via_product = copy(rho)
    via_qr = copy(rho)
    dmt_evolve!(via_product, product_evo)
    dmt_evolve!(via_qr, qr_evo)
    @test maximum(abs.(pauli_expectation_profile(via_product, terms) .-
      pauli_expectation_profile(via_qr, terms))) < 1e-6
    @test abs(inner(via_product, via_qr)) / (norm(via_product) * norm(via_qr)) > 1 - 1e-6
  end

  @testset "controlled gate schedule matches dense QR evolution" begin
    controlled_gates = [pauli_pxp_controlled_gate(psites, center, dt)
      for center in 1:nsites]
    dense_evo = DMTGateEvolution(gates, dt; schedule=schedule,
      reverse_schedule=reverse(schedule), maxdim=64, cutoff=0.0,
      gate_maxdim=0, gate_backend=:qr, normalize=false)
    controlled_evo = DMTGateEvolution(controlled_gates, dt; schedule=schedule,
      reverse_schedule=reverse(schedule), maxdim=64, cutoff=0.0,
      gate_maxdim=0, gate_backend=:controlled, normalize=false)
    via_dense = copy(rho)
    via_controlled = copy(rho)
    dmt_evolve!(via_dense, dense_evo)
    dmt_evolve!(via_controlled, controlled_evo)
    @test abs(inner(via_dense, via_controlled)) /
      (norm(via_dense) * norm(via_controlled)) ≈ 1.0 atol=1e-10
    @test pauli_expectation_profile(via_dense, terms) ≈
      pauli_expectation_profile(via_controlled, terms) atol=1e-10
    @test MPSToolkit._mps_eltype(via_controlled) === Float64

    plan = pxp_s2_dmt_plan(psites, dt; nstep=2, maxdim=64, cutoff=0.0,
      normalize=false)
    planned = copy(rho)
    run_state = ConstrainedDMTRunState()
    constrained_dmt_evolve!(planned, plan, projector; project_every=2,
      run_state=run_state, normalize=false)
    @test run_state.completed_steps == 2
    @test run_state.physical_time ≈ 2 * dt
    @test run_state.projection_count == 1
    @test MPSToolkit._mps_eltype(planned) === Float64

    gatewise_plan = pxp_s2_dmt_plan(psites, dt; maxdim=64, cutoff=0.0,
      normalize=false)
    layerwise_plan = pxp_layer_dmt_plan(psites, dt; scheme=:S2, maxdim=64,
      cutoff=0.0, normalize=false)
    gatewise = copy(rho)
    layerwise = copy(rho)
    dmt_evolve!(gatewise, gatewise_plan)
    dmt_evolve!(layerwise, layerwise_plan)
    @test abs(inner(gatewise, layerwise)) / (norm(gatewise) * norm(layerwise)) ≈
      1.0 atol=1e-10
    @test pauli_expectation_profile(gatewise, terms) ≈
      pauli_expectation_profile(layerwise, terms) atol=1e-10
    @test MPSToolkit._mps_eltype(layerwise) === Float64
  end

  @testset "persistent projection cadence is independent of call boundaries" begin
    make_evo(nstep) = DMTGateEvolution(
      gates, dt; schedule=schedule, reverse_schedule=reverse(schedule), nstep=nstep,
      maxdim=64, cutoff=0.0, gate_maxdim=256, normalize=false)

    whole = copy(rho)
    whole_state = ConstrainedDMTRunState()
    constrained_dmt_evolve!(whole, make_evo(10), projector; project_every=3,
      normalize=false, run_state=whole_state)

    split = copy(rho)
    split_state = ConstrainedDMTRunState()
    for _ in 1:10
      constrained_dmt_evolve!(split, make_evo(1), projector; project_every=3,
        normalize=false, run_state=split_state)
    end
    @test whole_state.completed_steps == split_state.completed_steps == 10
    @test whole_state.physical_time ≈ split_state.physical_time ≈ 20 * dt
    @test whole_state.steps_since_projection == split_state.steps_since_projection == 1
    @test whole_state.projection_count == split_state.projection_count == 3
    @test pauli_expectation_profile(whole, terms) ≈ pauli_expectation_profile(split, terms)

    normalized_whole = copy(rho)
    normalized_split = copy(rho)
    normalized_whole[1] *= 2
    normalized_split[1] *= 2
    normalized_whole_state = ConstrainedDMTRunState()
    normalized_split_state = ConstrainedDMTRunState()
    constrained_dmt_evolve!(normalized_whole, make_evo(4), projector; project_every=3,
      normalize=true, run_state=normalized_whole_state)
    for _ in 1:4
      constrained_dmt_evolve!(normalized_split, make_evo(1), projector; project_every=3,
        normalize=true, run_state=normalized_split_state)
    end
    @test pauli_expectation_profile(normalized_whole, terms) ≈
      pauli_expectation_profile(normalized_split, terms)

    dir = mktempdir()
    restarted = copy(rho)
    restart_state = ConstrainedDMTRunState()
    constrained_dmt_evolve!(restarted, make_evo(7), projector; project_every=3,
      normalize=false, run_state=restart_state)
    path = dmt_checkpoint_save(dmt_checkpoint_path("cadence", restart_state.physical_time;
      dir=dir), restarted, Dict(:project_every => 3, :dt => dt); time=restart_state.physical_time,
      sweep=restart_state.completed_steps, run_state=restart_state)
    checkpoint = dmt_checkpoint_load(path)
    restarted = checkpoint.state
    restart_state = dmt_checkpoint_run_state(checkpoint)
    constrained_dmt_evolve!(restarted, make_evo(3), projector; project_every=3,
      normalize=false, run_state=restart_state)
    @test restart_state.completed_steps == 10
    @test restart_state.projection_count == 3
    @test pauli_expectation_profile(whole, terms) ≈ pauli_expectation_profile(restarted, terms)

    constrained_dmt_evolve!(restarted, make_evo(1), projector; project_every=3,
      normalize=false, run_state=restart_state, final_project=true)
    @test restart_state.steps_since_projection == 0
    @test restart_state.projection_count == 4
  end

  @testset "chunk copies every DMT option" begin
    selected = [Matrix{Float64}(I, 2, 2)]
    configured = DMTGateEvolution(gates, dt; schedule=schedule,
      reverse_schedule=reverse(schedule), nstep=2, maxdim=12, cutoff=2e-10,
      gate_maxdim=19, preserve_diameter=3, preserve_operators=selected,
      truncation=:random, normalize=false)
    copied = MPSToolkit._copy_dmt_evolution(configured; nstep=1)
    for field in fieldnames(DMTGateEvolution)
      field == :nstep && continue
      @test getfield(copied, field) == getfield(configured, field)
    end
    @test copied.nstep == 1
  end

  @testset "direct leakage contraction agrees with an exact residual" begin
    projected = apply(projector, rho; maxdim=256, cutoff=0.0)
    residual_squared = real(inner(rho, rho) + inner(projected, projected) -
      2 * inner(rho, projected)) / real(inner(rho, rho))
    @test constraint_leakage_squared(rho, projector) ≈ residual_squared atol=1e-10
  end
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

@testset "constrained_dmt_evolve! honors the evo.normalize field (A2-1 regression)" begin
  # Audit follow-up: the constrained driver's `normalize` keyword must DEFAULT to evo.normalize,
  # not a hardcoded `true`. Previously a normalize=false evolution called WITHOUT an explicit
  # normalize keyword was silently re-normalized (the field was ignored), re-introducing the
  # trace-inflation footgun the field exists to prevent. At N=6 the operator space is exact at
  # maxdim=4^3=64, so unitary evolution + in-sector projection preserve the HS norm; scaling the
  # initial state to norm 2 makes a normalize=false run (norm stays 2) cleanly distinct from a
  # normalize=true run (norm reset to 1), independent of any truncation magnitude.
  nsites = 6
  psites = pauli_siteinds(nsites)
  terms = _pxp_terms(nsites)
  weights = [0.4, 0.4, 0.0, 0.0, -0.4, -0.4]
  rho = pauli_gibbs_state(
    psites, terms, weights;
    nsteps=8, maxdim=256, cutoff=0.0, initial_state=pauli_pxp_constraint_state(psites),
  )
  rho[1] = 2.0 * rho[1]                       # norm(rho) = 2, distinctly != 1
  gates = [pauli_gate_from_hamiltonian(h, 0.05) for (_, h) in terms]
  schedule = [start for (start, _) in terms]
  projector = pauli_pxp_constraint_projector(psites)
  # normalize=false carried ONLY in the field (no call-site keyword below); exact bond dim.
  evo = DMTGateEvolution(
    gates, 0.05;
    schedule=schedule, reverse_schedule=reverse(schedule),
    nstep=4, maxdim=64, cutoff=1e-12, gate_maxdim=256, normalize=false,
  )

  honored = copy(rho)
  constrained_dmt_evolve!(honored, evo, projector; project_every=2)                        # no kwarg
  explicit_false = copy(rho)
  constrained_dmt_evolve!(explicit_false, evo, projector; project_every=2, normalize=false)
  explicit_true = copy(rho)
  constrained_dmt_evolve!(explicit_true, evo, projector; project_every=2, normalize=true)

  @test norm(honored) ≈ norm(explicit_false) atol = 1e-10   # field honored: no-kwarg ≡ explicit false
  @test isapprox(norm(explicit_true), 1.0; atol = 1e-6)     # explicit true renormalizes (kwarg overrides)
  @test isapprox(norm(honored), 2.0; rtol = 1e-3)           # field false: the norm-2 scale is preserved
  @test !isapprox(norm(honored), 1.0; rtol = 1e-2)          # the regression: default no longer forces -> 1
end
