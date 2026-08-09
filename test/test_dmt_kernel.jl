using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Test

include("dmt_test_helpers.jl")

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

@testset "QR and SVD bond factorization give the same truncation" begin
  # The kernel needs orthonormal bases on both sides of the cut and the bond matrix expressed in
  # them -- never the Schmidt form -- so `psi[bond] = Q R` and `psi[bond] = U S V'` must produce
  # the SAME truncated operator, not merely two states that both satisfy the guarantee.
  Random.seed!(4242)
  for d in (2, 3)
    nsites, chi, bond = 6, 50, 3
    maxdim = 2 * d^2 + 10
    sites = operator_siteinds(nsites; d = d)
    base = add(operator_basis_state(sites, fill(1, nsites)),
      0.3 * random_mps(ComplexF64, sites; linkdims = chi); maxdim = chi + 1, cutoff = 0.0)
    # Two paths that both short-circuit and return the input unchanged would satisfy every
    # assertion below, so pin that the truncation actually fires.
    @test dim(linkind(base, bond)) > maxdim
    # A moderate probe set: this test compares the two factorizations against each other, while
    # the exhaustive sweep of the guarantee itself lives in `test_dmt_preservation.jl` and
    # `test_dmt_higher_spin.jl` (both of which now run the `:qr` default).
    probes = diameter_probes(nsites, d, 3; full_basis_width = 1, cap = 2, nrandom = 3)
    before = operator_expectation_profile(base, probes; normalize = false)

    via_svd = copy(base)
    MPSToolkit._dmt_bond_truncate!(via_svd, bond; maxdim = maxdim, cutoff = 1e-14,
      factorize = :svd)
    via_qr = copy(base)
    MPSToolkit._dmt_bond_truncate!(via_qr, bond; maxdim = maxdim, cutoff = 1e-14,
      factorize = :qr)

    @test dim(linkind(via_qr, bond)) == dim(linkind(via_svd, bond))
    @test dim(linkind(via_qr, bond)) <= maxdim
    @test abs(inner(via_svd, via_qr)) / (norm(via_svd) * norm(via_qr)) ≈ 1.0 atol = 1e-9
    @test norm(via_qr) ≈ norm(via_svd) rtol = 1e-12        # the overlap above fixes only direction
    # The overlap is a single number. Compare the whole diameter-3 observable profile as well,
    # so a discrepancy that happens to be small in the state's own direction cannot hide in it.
    svd_after = operator_expectation_profile(via_svd, probes; normalize = false)
    qr_after = operator_expectation_profile(via_qr, probes; normalize = false)
    @test preservation_error(svd_after, qr_after) < 1e-12
    # ... and the guarantee itself holds equally well on both paths.
    @test preservation_error(before, svd_after) < 1e-11
    @test preservation_error(before, qr_after) < 1e-11
  end
end

@testset "a wider-than-tall bond tensor falls back to the SVD" begin
  # `qr` returns a trapezoidal `R` when the bond tensor has fewer rows than columns, and the
  # square-bond kernel downstream cannot use it. Such a bond is rank deficient by construction --
  # here the left of the cut spans only 2 * d^2 = 8 directions while the cut itself carries 32 --
  # and the `:qr` request must quietly fall back rather than throw.
  d, nsites, bond = 2, 5, 2
  sites = operator_siteinds(nsites; d = d)
  links = [Index(k, "Link,l=$(j)") for (j, k) in enumerate((2, 40, 8, 4))]
  tensors = [random_itensor(ComplexF64, sites[1], links[1])]
  for site in 2:(nsites - 1)
    push!(tensors, random_itensor(ComplexF64, dag(links[site - 1]), sites[site], links[site]))
  end
  push!(tensors, random_itensor(ComplexF64, dag(links[nsites - 1]), sites[nsites]))
  base = MPS(tensors)
  maxdim = 2 * d^2 + 4

  # Pin that this input really does reach the rectangular branch: gauged to the bond, the left
  # dimension must be smaller than the bond dimension, and the truncation must fire.
  gauged = copy(base)
  orthogonalize!(gauged, bond)
  @test dim(linkind(gauged, bond - 1)) * d^2 < dim(linkind(gauged, bond))
  @test dim(linkind(gauged, bond)) > maxdim

  via_qr = copy(base)
  MPSToolkit._dmt_bond_truncate!(via_qr, bond; maxdim = maxdim, cutoff = 1e-14, factorize = :qr)
  via_svd = copy(base)
  MPSToolkit._dmt_bond_truncate!(via_svd, bond; maxdim = maxdim, cutoff = 1e-14, factorize = :svd)
  @test dim(linkind(via_qr, bond)) == dim(linkind(via_svd, bond))
  @test dim(linkind(via_qr, bond)) <= maxdim
  @test abs(inner(via_svd, via_qr)) / (norm(via_svd) * norm(via_qr)) ≈ 1.0 atol = 1e-9
  @test operator_trace(via_qr) ≈ operator_trace(base) rtol = 1e-11
end

@testset "the complement machinery never forms a chi x chi matrix" begin
  # This is the promise of the basis-free kernel, pinned where it is unambiguous: at chi = 6000 a
  # single dense chi x chi ComplexF64 matrix is 576 MB, while the matrix-free path only ever
  # touches chi x (rank + oversample) blocks -- measured 250 MB, and falling as 1/chi relative to
  # the bound (0.46x here, 0.68x at chi = 4000). `:dense` is excluded by construction: it
  # materializes the operator and calls LAPACK on it, so the guard is on `:random`.
  chi, rank_target = 6000, 40
  bond_matrix = Diagonal(ComplexF64.(sort(rand(chi); rev = true)))
  QL = MPSToolkit._protected_basis(randn(ComplexF64, chi, 4))
  QR = MPSToolkit._protected_basis(randn(ComplexF64, chi, 4))
  a, b, _ = MPSToolkit._dmt_connector(bond_matrix, QL[:, 1], QR[:, 1])
  function complement_svd()
    ops = MPSToolkit._dmt_complement_ops(bond_matrix, a, b, QL, QR)
    return MPSToolkit._truncated_svd(ops.mul, ops.adj, chi, rank_target; mode = :random)
  end
  complement_svd()                                  # warm up: compile, do not measure
  bytes = @allocated complement_svd()
  @info "complement machinery at chi = $(chi), rank $(rank_target): $(bytes) bytes " *
        "($(round(bytes / (16 * chi^2); digits = 2)) x a dense chi x chi ComplexF64 matrix)"
  @test bytes < 16 * chi^2
end

@testset "one bond step allocates a bounded multiple of the bond tensor" begin
  # A bond step cannot allocate less than the bond tensor it is handed: that tensor is
  # (chi d^2) x chi, i.e. `d^2` chi x chi matrices already, and its Q factor is as many again. So
  # the honest coarse bound is a multiple of *that*, not of chi^2 -- the sharp "never a chi x chi
  # matrix" statement is the testset above, where it can be made without competing against the
  # input. Measured on Julia 1.12 / ITensors 0.9, both truncation modes land at 17-19x the bond
  # tensor, against the 36x bound here.
  d, nsites, bond, chi = 2, 10, 5, 200
  maxdim = 2 * d^2 + 60
  sites = operator_siteinds(nsites; d = d)
  psi = random_mps(ComplexF64, sites; linkdims = chi)
  normalize!(psi)
  # Gauging is not part of the working set under test and would dominate it, so it happens once,
  # here, and the kernel runs with `orthogonalize = false`.
  orthogonalize!(psi, bond)
  @test dim(linkind(psi, bond)) > maxdim            # the truncation really fires
  bond_tensor_bytes = 16 * d^2 * chi^2
  for truncation in (:dense, :random)
    step!(state) = MPSToolkit._dmt_bond_truncate!(state, bond; maxdim = maxdim, cutoff = 1e-14,
      truncation = truncation, orthogonalize = false)
    step!(copy(psi))                                # warm up: compile, do not measure
    bytes = @allocated step!(copy(psi))
    @info "one bond step at d = $(d), chi = $(chi), maxdim = $(maxdim), $(truncation): " *
          "$(bytes) bytes ($(round(bytes / bond_tensor_bytes; digits = 2)) x the bond tensor)"
    @test bytes < 36 * bond_tensor_bytes
  end
end

@testset "an unknown factorization is rejected" begin
  sites = operator_siteinds(4; d = 2)
  psi = random_mps(ComplexF64, sites; linkdims = 20)
  # Rejected even though this bond is under budget and short-circuits before factorizing.
  @test_throws ArgumentError MPSToolkit._dmt_bond_truncate!(
    psi, 2; maxdim = 400, cutoff = 0.0, factorize = :lu)
end
