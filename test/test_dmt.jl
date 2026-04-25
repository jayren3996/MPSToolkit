using ITensors
using ITensorMPS
using LinearAlgebra
using Test

function _dmt_test_state(nsites)
  sites = pauli_siteinds(nsites)
  labels = [isodd(n) ? "Z" : "I" for n in 1:nsites]
  state = pauli_basis_state(sites, labels)
  normalize!(state)
  return sites, state
end

function _identity_gate(window)
  dim = 4^window
  return Matrix{ComplexF64}(I, dim, dim)
end

function _manual_dmt_sweep!(psi, gates; maxdim, cutoff, gate_maxdim, connector_buffer)
  for (bond, gate) in enumerate(gates)
    dmt_step!(
      psi,
      gate,
      bond;
      maxdim=maxdim,
      cutoff=cutoff,
      direction=:R,
      gate_maxdim=gate_maxdim,
      connector_buffer=connector_buffer,
    )
  end
  for (offset, gate) in enumerate(reverse(gates))
    bond = length(gates) - offset + 1
    dmt_step!(
      psi,
      gate,
      bond;
      maxdim=maxdim,
      cutoff=cutoff,
      direction=:L,
      gate_maxdim=gate_maxdim,
      connector_buffer=connector_buffer,
    )
  end
  normalize!(psi)
  return psi
end

@testset "operator-space DMT" begin
  @test isdefined(MPSToolkit, :dmt_step!)
  @test isdefined(MPSToolkit, :dmt_evolve!)
  @test isdefined(MPSToolkit, :DMTGateEvolution)
  @test MPSToolkit.OperatorSpace.dmt_step! === MPSToolkit.dmt_step!
  @test MPSToolkit.OperatorSpace.dmt_evolve! === MPSToolkit.dmt_evolve!
  @test MPSToolkit.OperatorSpace.DMTGateEvolution === MPSToolkit.DMTGateEvolution

  @testset "identity DMT step preserves a product operator" begin
    _, psi = _dmt_test_state(4)
    reference = copy(psi)

    dmt_step!(psi, _identity_gate(2), 2; maxdim=8, cutoff=1e-12, direction=:R, gate_maxdim=64)

    @test inner(reference, psi) ≈ 1.0 atol = 1e-10
  end

  @testset "scheduled DMT evolution matches explicit two-site sweep" begin
    _, manual = _dmt_test_state(5)
    scheduled = copy(manual)
    x = ComplexF64[0 1; 1 0]
    zz = ComplexF64[1 0; 0 -1]
    gate = pauli_gate(exp(-0.1im * kron(x, zz)))
    gates = [gate, gate, gate, gate]

    _manual_dmt_sweep!(manual, gates; maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=8)
    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=4, cutoff=1e-12, gate_maxdim=64)
    dmt_evolve!(scheduled, evo)

    @test inner(manual, scheduled) ≈ 1.0 atol = 1e-8
  end

  @testset "scheduled DMT evolution matches explicit three-site sweep" begin
    _, manual = _dmt_test_state(6)
    scheduled = copy(manual)
    x = ComplexF64[0 1; 1 0]
    z = ComplexF64[1 0; 0 -1]
    gate = pauli_gate(exp(-0.05im * kron(kron(x, z), x)))
    gates = [gate, gate, gate, gate]

    _manual_dmt_sweep!(manual, gates; maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=8)
    evo = DMTGateEvolution(gates, 0.05; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=4, cutoff=1e-12, gate_maxdim=64)
    dmt_evolve!(scheduled, evo)

    @test inner(manual, scheduled) ≈ 1.0 atol = 1e-8
  end

  @testset "reverse DMT sweep uses gates associated with original bonds" begin
    sites = pauli_siteinds(4)
    state_a = pauli_basis_state(sites, [2, 1, 1, 1])
    state_b = pauli_basis_state(sites, [1, 1, 2, 1])
    manual = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))
    scheduled = copy(manual)

    gate1 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 2.0 for n in 1:16])
    gate2 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 0.5 for n in 1:16])
    gate3 = Matrix{ComplexF64}(I, 16, 16)
    gates = [gate1, gate2, gate3]

    _manual_dmt_sweep!(manual, gates; maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 3], reverse_schedule=[3, 2, 1], maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    dmt_evolve!(scheduled, evo)

    @test abs(inner(manual, scheduled)) ≈ 1.0 atol = 1e-10
  end

  @testset "reverse DMT sweep maps repeated forward schedule entries by position" begin
    sites = pauli_siteinds(3)
    state_a = pauli_basis_state(sites, [2, 1, 1])
    state_b = pauli_basis_state(sites, [1, 2, 1])
    manual = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))
    scheduled = copy(manual)

    gate1 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 2.0 for n in 1:16])
    gate2 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 0.5 for n in 1:16])
    gate3 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 1.5 for n in 1:16])
    gates = [gate1, gate2, gate3]

    for (bond, gate) in zip([1, 2, 1], gates)
      dmt_step!(manual, gate, bond; maxdim=8, cutoff=1e-12, direction=:R, gate_maxdim=16, connector_buffer=0)
    end
    for (bond, gate) in zip([1, 2, 1], reverse(gates))
      dmt_step!(manual, gate, bond; maxdim=8, cutoff=1e-12, direction=:L, gate_maxdim=16, connector_buffer=0)
    end
    normalize!(manual)

    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 1], reverse_schedule=[1, 2, 1], maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    dmt_evolve!(scheduled, evo)

    @test abs(inner(manual, scheduled)) ≈ 1.0 atol = 1e-10
    ambiguous = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 1], reverse_schedule=[2, 1, 1], maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    @test_throws ArgumentError dmt_evolve!(copy(scheduled), ambiguous)
  end

  @testset "reverse DMT sweep maps callable gates by forward position" begin
    sites = pauli_siteinds(4)
    state_a = pauli_basis_state(sites, [2, 1, 1, 1])
    state_b = pauli_basis_state(sites, [1, 1, 2, 1])
    manual = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))
    scheduled = copy(manual)

    gate1 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 2.0 for n in 1:16])
    gate2 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 0.5 for n in 1:16])
    gate3 = Matrix{ComplexF64}(I, 16, 16)
    gates = [gate1, gate2, gate3]

    _manual_dmt_sweep!(manual, gates; maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    evo = DMTGateEvolution((bond, index) -> gates[index], 0.1; schedule=[1, 2, 3], reverse_schedule=[3, 2, 1], maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    dmt_evolve!(scheduled, evo)

    @test abs(inner(manual, scheduled)) ≈ 1.0 atol = 1e-10
  end
end
