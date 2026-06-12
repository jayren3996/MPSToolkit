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
