using ITensors
using ITensorMPS

function _transverse_field_mpo(sites; h=1.0)
  opsum = OpSum()
  for j in 1:length(sites)
    opsum += h, "Sx", j
  end
  return MPO(opsum, sites)
end

function _pxp_mpo(sites)
  nsites_total = length(sites)
  opsum = OpSum()
  opsum += 1.0, "Sx", 1, "ProjUp", 2
  for j in 2:(nsites_total - 1)
    opsum += 1.0, "ProjUp", j - 1, "Sx", j, "ProjUp", j + 1
  end
  opsum += 1.0, "ProjUp", nsites_total - 1, "Sx", nsites_total
  return MPO(opsum, sites)
end

function _longitudinal_field_mpo(sites; h=1.0)
  opsum = OpSum()
  for j in 1:length(sites)
    opsum += h, "Sz", j
  end
  return MPO(opsum, sites)
end

function _tilted_field_mpo(sites; hz=1.0, hx=0.3)
  opsum = OpSum()
  for j in 1:length(sites)
    opsum += hz, "Sz", j
    opsum += hx, "Sx", j
  end
  return MPO(opsum, sites)
end

@testset "finite tdvp evolve" begin
  sites = siteinds("S=1/2", 4)
  hamiltonian = _transverse_field_mpo(sites)

  @testset "zero time leaves state unchanged" begin
    psi = MPS(sites, n -> "Up")
    evo = TDVPEvolution(hamiltonian, 0.0; nsteps=1, updater_backend="exponentiate")
    evolve!(psi, evo)
    @test expect(psi, "Sz")[1] ≈ 0.5 atol = 1e-10
    @test maxlinkdim(psi) == 1
  end

  @testset "real evolution changes observables" begin
    psi = MPS(sites, n -> "Up")
    evo = TDVPEvolution(hamiltonian, -0.1im; time_step=-0.05im, nsteps=2, normalize=true, solver_kwargs=(maxdim=16, cutoff=1e-12))
    evolve!(psi, evo)
    @test expect(psi, "Sz")[1] < 0.5
  end

  @testset "mpo energy density matches direct expectation value" begin
    psi = MPS(sites, n -> "Up")
    expected = real(inner(psi', hamiltonian, psi)) / length(sites)
    @test energy_density(psi, hamiltonian) ≈ expected atol = 1e-10
  end

  @testset "mpo target energy correction can move energy toward zero" begin
    sites = siteinds("S=1/2", 6)
    psi = MPS(sites, n -> "Up")
    hamiltonian = _tilted_field_mpo(sites; hz=1.0, hx=0.3)
    evolution = TDVPEvolution(
      hamiltonian,
      0.0im;
      time_step=0.0im,
      nsteps=1,
      normalize=true,
      solver_kwargs=(maxdim=32, cutoff=1e-10),
    )
    truncation = BondDimTruncation(2; cutoff=1e-10)
    target = EnergyTarget(0.0; operator=hamiltonian, tol=1e-8, alpha=0.1, maxstep=5)
    initial_energy = energy_density(psi, hamiltonian)

    scarfinder_step!(psi, evolution, truncation; target_energy=target)

    @test abs(energy_density(psi, hamiltonian)) < abs(initial_energy)
  end
end
