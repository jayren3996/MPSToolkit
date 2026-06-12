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
