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

function _link_dims(psi)
  return [dim(linkind(psi, bond)) for bond in 1:(length(psi) - 1)]
end

function _dense_pauli_coefficients(psi)
  sites = [siteind(psi, n) for n in 1:length(psi)]
  coeffs = ComplexF64[]
  for labels in Iterators.product(ntuple(_ -> 1:4, length(psi))...)
    push!(coeffs, inner(pauli_basis_state(sites, collect(labels)), psi))
  end
  return coeffs
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
    @test dmt_step!(psi, _identity_gate(2), 2; maxdim=8, cutoff=1e-12, direction=:R, gate_maxdim=64) === psi
  end

  @testset "DMT truncates an enlarged first bond" begin
    sites = pauli_siteinds(3)
    terms = [
      pauli_basis_state(sites, [1, 1, 1]),
      pauli_basis_state(sites, [2, 2, 1]),
      pauli_basis_state(sites, [3, 3, 1]),
    ]
    psi = normalize(add(terms...; maxdim=8, cutoff=1e-14))

    @test dim(linkind(psi, 1)) > 1
    dmt_step!(psi, _identity_gate(2), 1; maxdim=1, cutoff=1e-12, direction=:R, gate_maxdim=8, connector_buffer=0)

    @test dim(linkind(psi, 1)) <= 1
    @test isfinite(real(inner(psi, psi)))
  end

  @testset "DMT validates connector buffer budget" begin
    @test_throws ArgumentError DMTGateEvolution(_identity_gate(2), 0.1; schedule=[1], maxdim=2, connector_buffer=3)

    _, psi = _dmt_test_state(3)
    @test_throws ArgumentError dmt_step!(psi, _identity_gate(2), 1; maxdim=2, connector_buffer=3)
  end

  @testset "DMT completes protected bases beyond local Pauli dimension" begin
    sites = pauli_siteinds(6)
    psi = random_mps(sites; linkdims=20)
    normalize!(psi)

    @test dim(linkind(psi, 3)) == 20
    dmt_step!(psi, _identity_gate(2), 3; maxdim=12, cutoff=1e-12, gate_maxdim=40, connector_buffer=8)

    @test dim(linkind(psi, 3)) == 12
  end

  @testset "DMT preserves identity and local Pauli data under truncation" begin
    sites = pauli_siteinds(5)
    terms = [
      pauli_basis_state(sites, fill(1, 5); coefficient=1.0),
      pauli_basis_state(sites, [1, 4, 1, 1, 1]; coefficient=0.3),
      pauli_basis_state(sites, fill(2, 5); coefficient=0.2),
    ]
    psi = add(terms...; maxdim=8, cutoff=1e-14)
    probes = [fill(1, 5), [1, 4, 1, 1, 1], [1, 1, 4, 1, 1], [1, 1, 1, 4, 1]]
    before = [inner(pauli_basis_state(sites, labels), psi) for labels in probes]

    dmt_step!(psi, _identity_gate(3), 2; maxdim=1, cutoff=1e-12, direction=:R, gate_maxdim=8, connector_buffer=0)

    after = [inner(pauli_basis_state(sites, labels), psi) for labels in probes]
    @test after ≈ before atol = 1e-12
  end

  @testset "complex DMT projection uses adjoint orthonormal bases" begin
    left_protected = ComplexF64[
      1 1+im
      im 2
      0.5-im -0.25
      0.1 0.3im
    ]
    right_protected = ComplexF64[
      1-im 0.2
      0.7 0.4im
      im 1
      0.3+0.2im -0.1
    ]
    left_basis = MPSToolkit._complete_orthonormal_basis(left_protected, 4)
    right_basis = MPSToolkit._complete_orthonormal_basis(right_protected, 4)
    singular_data = ComplexF64[
      1.0 0.2im 0.1 0.3
      -0.4im 0.8 0.2+0.1im 0.0
      0.1im -0.2 0.5 0.3im
      0.0 0.1 0.2im 0.2
    ]

    reduced = left_basis' * singular_data * right_basis
    repaired = left_basis * reduced * right_basis'

    @test left_basis' * left_basis ≈ Matrix{ComplexF64}(I, 4, 4) atol = 1e-12
    @test right_basis' * right_basis ≈ Matrix{ComplexF64}(I, 4, 4) atol = 1e-12
    @test repaired ≈ singular_data atol = 1e-12
  end

  @testset "DMT truncates every internal bond of wider update windows" begin
    for span in (1, 2, 3, 4, 5)
      nsites = span + 3
      sites = pauli_siteinds(nsites)
      psi = random_mps(sites; linkdims=16)
      normalize!(psi)
      start = 2

      dmt_step!(psi, _identity_gate(span), start; maxdim=10, cutoff=1e-12, direction=:R, gate_maxdim=32, connector_buffer=4)

      for bond in start:(start + span - 2)
        @test dim(linkind(psi, bond)) <= 10
      end
    end
  end

  @testset "invalid DMT calls do not mutate the state" begin
    sites = pauli_siteinds(5)
    psi = random_mps(sites; linkdims=8)
    normalize!(psi)
    reference = copy(psi)
    reference_dims = _link_dims(reference)

    @test_throws ArgumentError dmt_step!(psi, _identity_gate(3), 2; direction=:bad, maxdim=4, gate_maxdim=16)
    @test _link_dims(psi) == reference_dims
    @test inner(reference, psi) ≈ inner(reference, reference) atol = 1e-12

    @test_throws ArgumentError dmt_step!(psi, _identity_gate(3), 4; direction=:R, maxdim=4, gate_maxdim=16)
    @test _link_dims(psi) == reference_dims
    @test inner(reference, psi) ≈ inner(reference, reference) atol = 1e-12

    @test_throws ArgumentError dmt_step!(psi, _identity_gate(2), length(psi); direction=:R, maxdim=4, gate_maxdim=16)
    @test _link_dims(psi) == reference_dims
    @test inner(reference, psi) ≈ inner(reference, reference) atol = 1e-12

    bad_sites = [Index(3, "NotPauli,n=$n") for n in 1:4]
    bad_psi = random_mps(bad_sites; linkdims=4)
    bad_reference_dims = _link_dims(bad_psi)
    @test_throws ArgumentError dmt_step!(bad_psi, Matrix{ComplexF64}(I, 9, 9), 2; maxdim=2, gate_maxdim=4)
    @test _link_dims(bad_psi) == bad_reference_dims
  end

  @testset "scheduled DMT evolution matches explicit two-site sweep" begin
    _, manual = _dmt_test_state(5)
    scheduled = copy(manual)
    x = ComplexF64[0 1; 1 0]
    zz = ComplexF64[1 0; 0 -1]
    gate = pauli_gate(exp(-0.1im * kron(x, zz)))
    gates = [gate, gate, gate, gate]

    _manual_dmt_sweep!(manual, gates; maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
    @test dmt_evolve!(scheduled, evo) === scheduled

    @test inner(manual, scheduled) ≈ 1.0 atol = 1e-8
  end

  @testset "scheduled DMT evolution matches explicit three-site sweep" begin
    _, manual = _dmt_test_state(6)
    scheduled = copy(manual)
    x = ComplexF64[0 1; 1 0]
    z = ComplexF64[1 0; 0 -1]
    gate = pauli_gate(exp(-0.05im * kron(kron(x, z), x)))
    gates = [gate, gate, gate, gate]

    _manual_dmt_sweep!(manual, gates; maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
    evo = DMTGateEvolution(gates, 0.05; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
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

    sx = ComplexF64[0 1; 1 0]
    sy = ComplexF64[0 -im; im 0]
    sz = ComplexF64[1 0; 0 -1]
    gate1 = pauli_gate(exp(-0.07im * (kron(sx, sz) + 0.3 * kron(sz, sx))))
    gate2 = pauli_gate(exp(-0.11im * (kron(sy, sx) + 0.2 * kron(sz, sz))))
    gate3 = pauli_gate(exp(-0.05im * (kron(sz, sy) + 0.4 * kron(sx, sx))))
    gates = [gate1, gate2, gate3]

    _manual_dmt_sweep!(manual, gates; maxdim=32, cutoff=1e-12, gate_maxdim=64, connector_buffer=0)
    evo = DMTGateEvolution((bond, index) -> gates[index], 0.1; schedule=[1, 2, 3], reverse_schedule=[3, 2, 1], maxdim=32, cutoff=1e-12, gate_maxdim=64, connector_buffer=0)
    dmt_evolve!(scheduled, evo)

    @test _dense_pauli_coefficients(scheduled) ≈ _dense_pauli_coefficients(manual) atol = 1e-10
  end

  @testset "custom repeated reverse DMT schedule rejects ambiguous callable indices" begin
    sites = pauli_siteinds(3)
    state_a = pauli_basis_state(sites, [2, 1, 1])
    state_b = pauli_basis_state(sites, [1, 2, 1])
    scheduled = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))

    gate = Diagonal(ComplexF64[isodd(n) ? 1.0 : 1.2 for n in 1:16])
    forward_schedule = [1, 2, 1]
    reverse_schedule = [2, 1, 1]

    evo = DMTGateEvolution((bond, index) -> gate, 0.1; schedule=forward_schedule, reverse_schedule=reverse_schedule, maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    @test_throws ArgumentError dmt_evolve!(scheduled, evo)
  end
end
