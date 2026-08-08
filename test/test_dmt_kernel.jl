using LinearAlgebra
using MPSToolkit
using Random
using Test

@testset "DMT truncation primitives" begin
  Random.seed!(20260808)

  @testset "_protected_basis is orthonormal and spans the input" begin
    for (chi, k) in ((40, 4), (40, 9), (12, 16))
      protected = randn(ComplexF64, chi, k)
      basis = MPSToolkit._protected_basis(protected)
      @test size(basis, 1) == chi
      @test size(basis, 2) == min(chi, k)
      @test basis' * basis ≈ Matrix{ComplexF64}(I, size(basis, 2), size(basis, 2)) atol = 1e-12
      if k <= chi
        @test norm(protected - basis * (basis' * protected)) < 1e-10
      end
    end
  end

  @testset "_protected_basis stays orthonormal on a rank-deficient block" begin
    # `randn` inputs are full rank and do not exercise the degenerate case. This function's
    # predecessor (`_complete_orthonormal_basis`) shipped a real orthogonality bug that random
    # inputs did not catch, so duplicate / scaled / zero columns are tested explicitly: the
    # basis must stay orthonormal AND still span every input column.
    chi = 24
    base = randn(ComplexF64, chi, 2)
    protected = hcat(base[:, 1], base[:, 1], base[:, 2], zeros(ComplexF64, chi), 3.0 * base[:, 1])
    @test rank(protected; atol = 1e-10) == 2      # genuinely rank deficient
    basis = MPSToolkit._protected_basis(protected)
    @test size(basis) == (chi, 5)
    @test basis' * basis ≈ Matrix{ComplexF64}(I, 5, 5) atol = 1e-12
    @test norm(protected - basis * (basis' * protected)) < 1e-10
  end

  @testset "the connector is rank one and annihilated by B" begin
    for bond_matrix in (Diagonal(ComplexF64.(sort(rand(30); rev = true))),
                        UpperTriangular(randn(ComplexF64, 30, 30)))
      chi = size(bond_matrix, 1)
      q0 = normalize(randn(ComplexF64, chi))
      r0 = normalize(randn(ComplexF64, chi))
      a, b, has_connector = MPSToolkit._dmt_connector(bond_matrix, q0, r0)
      @test has_connector
      # b is the row vector q0' S, not S' q0
      @test transpose(b) ≈ q0' * Matrix(bond_matrix) atol = 1e-10
      connector = a * transpose(b)
      @test rank(connector; atol = 1e-10) == 1
      bmat = Matrix(bond_matrix) - connector
      # B annihilates the identity directions: this is why the connector needs no extra budget.
      @test norm(bmat * r0) < 1e-10
      @test norm(q0' * bmat) < 1e-10
    end
  end

  @testset "a negligible trace overlap disables the connector rather than blowing up" begin
    chi = 20
    bond_matrix = Diagonal(ComplexF64.(sort(rand(chi); rev = true)))
    q0 = zeros(ComplexF64, chi); q0[1] = 1
    r0 = zeros(ComplexF64, chi); r0[2] = 1     # orthogonal directions => q0' S r0 == 0
    a, b, has_connector = MPSToolkit._dmt_connector(bond_matrix, q0, r0)
    @test !has_connector
    @test all(iszero, a)
  end

  @testset "the complement operator is doubly orthogonal to the protected subspaces" begin
    chi, k = 40, 9
    bond_matrix = Diagonal(ComplexF64.(sort(rand(chi); rev = true)))
    QL = MPSToolkit._protected_basis(randn(ComplexF64, chi, k))
    QR = MPSToolkit._protected_basis(randn(ComplexF64, chi, k))
    q0 = QL[:, 1]
    r0 = QR[:, 1]
    a, b, _ = MPSToolkit._dmt_connector(bond_matrix, q0, r0)
    ops = MPSToolkit._dmt_complement_ops(bond_matrix, a, b, QL, QR)
    dense = ops.mul(Matrix{ComplexF64}(I, chi, chi))
    @test norm(QL' * dense) < 1e-10        # rows orthogonal to the left protected space
    @test norm(dense * QR) < 1e-10         # columns orthogonal to the right protected space
    # adj really is the adjoint
    x = randn(ComplexF64, chi, 3)
    y = randn(ComplexF64, chi, 3)
    @test tr(y' * ops.mul(x)) ≈ tr((ops.adj(y))' * x) atol = 1e-10
  end

  @testset "randomized truncated SVD matches the dense one" begin
    chi, rank_target = 200, 20
    base = randn(ComplexF64, chi, chi)
    mul(x) = base * x
    adj(x) = base' * x
    u_dense, s_dense, v_dense = MPSToolkit._truncated_svd(mul, adj, chi, rank_target; mode = :dense)
    u_rand, s_rand, v_rand = MPSToolkit._truncated_svd(mul, adj, chi, rank_target; mode = :random)
    @test length(s_dense) == rank_target
    @test length(s_rand) == rank_target
    @test s_rand ≈ s_dense rtol = 0.05
    dense_error = norm(base - u_dense * Diagonal(s_dense) * v_dense')
    rand_error = norm(base - u_rand * Diagonal(s_rand) * v_rand')
    @test rand_error <= 1.05 * dense_error
  end

  @testset "_dmt_refactor reproduces the dense SVD of a low-rank product" begin
    chi, r = 300, 25
    f = randn(ComplexF64, chi, r)
    g = randn(ComplexF64, chi, r)
    u, s, v = MPSToolkit._dmt_refactor(f, g, 1000, 1e-14)
    @test size(u, 2) == length(s) == size(v, 2)
    @test length(s) <= r
    @test u * Diagonal(s) * v' ≈ f * g' atol = 1e-9
    @test u' * u ≈ Matrix{ComplexF64}(I, length(s), length(s)) atol = 1e-10
    # maxdim clips
    u2, s2, v2 = MPSToolkit._dmt_refactor(f, g, 10, 0.0)
    @test length(s2) == 10
  end
end
