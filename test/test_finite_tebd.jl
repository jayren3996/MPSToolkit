using ITensors
using ITensorMPS
using LinearAlgebra

@testset "finite tebd evolve" begin
  s = siteinds("S=1/2", 4)
  psi = MPS(s, n -> "↑")
  gate = [1.0 0.0 0.0 0.0;
          0.0 1.0 0.0 0.0;
          0.0 0.0 1.0 0.0;
          0.0 0.0 0.0 1.0]
  evo = LocalGateEvolution(gate, 0.01; schedule=[1, 2, 3], nstep=1, maxdim=4, cutoff=1e-12)

  evolve!(psi, evo)

  @test length(psi) == 4
  @test maxlinkdim(psi) <= 1
end

@testset "finite one-site tebd evolve" begin
  s = siteinds("S=1/2", 4)
  psi = MPS(s, n -> "Up")
  gate = ComplexF64[0 1; 1 0]
  evo = LocalGateEvolution(gate, 0.1; schedule=[2], nstep=1, maxdim=4, cutoff=1e-12)

  evolve!(psi, evo)

  sz = expect(psi, "Sz")
  @test sz[1] ≈ 0.5 atol = 1e-10
  @test sz[2] ≈ -0.5 atol = 1e-10
  @test sz[3] ≈ 0.5 atol = 1e-10
  @test sz[4] ≈ 0.5 atol = 1e-10
end

@testset "finite three-site tebd evolve" begin
  s = siteinds("S=1/2", 5)
  psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  gate = diagm(ones(8))
  evo = LocalGateEvolution(gate, 0.1; schedule=[1, 2, 3], nstep=1, maxdim=8, cutoff=1e-12)

  evolve!(psi, evo)

  @test length(psi) == 5
  @test maxlinkdim(psi) == 1
end

@testset "finite tebd supports per-bond gates" begin
  s = siteinds("S=1/2", 5)
  psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
  gates = [diagm(ones(8)), diagm(ones(8))]
  evo = LocalGateEvolution(gates, 0.1; schedule=[1, 2], nstep=1, maxdim=8, cutoff=1e-12)

  evolve!(psi, evo)

  @test length(psi) == 5
  @test maxlinkdim(psi) == 1
end

@testset "finite tebd supports periodic boundary gates" begin
  s = siteinds("S=1/2", 4)
  psi = MPS(s, n -> "Up")
  x = ComplexF64[0 1; 1 0]
  gate = kron(x, x)
  evo = LocalGateEvolution(gate, 0.1; schedule=[4], nstep=1, maxdim=8, cutoff=1e-12)

  evolve!(psi, evo)

  sz = expect(psi, "Sz")
  @test sz[1] ≈ -0.5 atol = 1e-10
  @test sz[2] ≈ 0.5 atol = 1e-10
  @test sz[3] ≈ 0.5 atol = 1e-10
  @test sz[4] ≈ -0.5 atol = 1e-10
end

@testset "local gates from Hamiltonians" begin
  h = ComplexF64[1 0; 0 -1]
  gate = local_gates_from_hamiltonians(h, 0.2)

  @test gate ≈ exp(-0.2im * h)
end

@testset "local gates from Hamiltonians supports custom mapping" begin
  hs = [reshape(1:4, 2, 2), reshape(5:8, 2, 2)]
  gates = local_gates_from_hamiltonians(hs, 0.1; map_hamiltonian=(h, step_dt) -> step_dt .* h)

  @test gates[1] == 0.1 .* hs[1]
  @test gates[2] == 0.1 .* hs[2]
end

@testset "tebd evolution from Hamiltonians" begin
  s = siteinds("S=1/2", 4)
  psi = MPS(s, n -> "Up")
  h = ComplexF64[0 1; 1 0]
  evo = tebd_evolution_from_hamiltonians(h, pi / 2; schedule=[2], nstep=1, maxdim=4, cutoff=1e-12)

  evolve!(psi, evo)

  sz = expect(psi, "Sz")
  @test sz[1] ≈ 0.5 atol = 1e-10
  @test sz[2] ≈ -0.5 atol = 1e-10
  @test sz[3] ≈ 0.5 atol = 1e-10
  @test sz[4] ≈ 0.5 atol = 1e-10
end

@testset "Strang schedule helper" begin
  schedule, weights = tebd_strang_schedule(6)
  @test schedule == [1, 3, 5, 2, 4, 1, 3, 5]
  @test weights == [0.5, 0.5, 0.5, 1.0, 1.0, 0.5, 0.5, 0.5]
end

@testset "Strang evolution helper" begin
  s = siteinds("S=1/2", 4)
  psi = MPS(s, n -> "Up")
  evo = tebd_strang_evolution(4, 0.3; local_hamiltonian=(bond, weight) -> zeros(ComplexF64, 2, 2), maxdim=4, cutoff=1e-12)
  evolve!(psi, evo)
  sz = expect(psi, "Sz")
  @test evo.schedule == [1, 3, 2, 1, 3]
  @test all(isapprox(value, 0.5; atol=1e-10) for value in sz)
end

@testset "finite energy density" begin
  s = siteinds("S=1/2", 4)
  psi = MPS(s, n -> "↑")
  id2 = [1.0 0.0 0.0 0.0;
         0.0 1.0 0.0 0.0;
         0.0 0.0 1.0 0.0;
         0.0 0.0 0.0 1.0]

  @test energy_density(psi, id2) ≈ 1.0 atol = 1e-10
end
