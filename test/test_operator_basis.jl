using ITensors
using LinearAlgebra
using MPSToolkit
using Test

@testset "generic onsite operator basis" begin
  @testset "d = 2 reproduces the normalized Pauli basis exactly" begin
    basis = operator_basis_matrices(2)
    expected = [m / sqrt(2) for m in values(pauli_matrices())]   # (I, X, Y, Z) / sqrt(2)
    @test length(basis) == 4
    for k in 1:4
      @test basis[k] ≈ expected[k] atol = 1e-15
    end
  end

  @testset "orthonormal, Hermitian, identity first" begin
    for d in 2:5
      basis = operator_basis_matrices(d)
      @test length(basis) == d^2
      @test basis[1] ≈ Matrix{ComplexF64}(I, d, d) / sqrt(d) atol = 1e-14
      gram = [tr(basis[i]' * basis[j]) for i in eachindex(basis), j in eachindex(basis)]
      @test gram ≈ Matrix{ComplexF64}(I, d^2, d^2) atol = 1e-12
      for m in basis
        @test m ≈ m' atol = 1e-14
      end
      # every element other than the identity is traceless
      for k in 2:length(basis)
        @test abs(tr(basis[k])) < 1e-14
      end
    end
  end

  @testset "the cache returns the same matrices, not fresh copies" begin
    @test operator_basis_matrices(3) === operator_basis_matrices(3)
  end

  @testset "local_dimension inverts the site dimension" begin
    @test local_dimension(4) == 2
    @test local_dimension(9) == 3
    @test local_dimension(Index(16, "OperatorSpace,n=1")) == 4
    @test_throws ArgumentError local_dimension(8)     # not a perfect square
    @test_throws ArgumentError local_dimension(1)     # d = 1 is not a spin space
    @test_throws ArgumentError operator_basis_matrices(1)
  end
end
