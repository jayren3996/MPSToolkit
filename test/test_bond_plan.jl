using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Test

include("dmt_test_helpers.jl")

const _BOND_I2 = Matrix{ComplexF64}(I, 2, 2)

function _bond_embed(op, bond, nsites)
  left = bond == 1 ? Matrix{ComplexF64}(I, 1, 1) :
    foldl(kron, fill(_BOND_I2, bond - 1))
  right_count = nsites - bond - 1
  right = right_count == 0 ? Matrix{ComplexF64}(I, 1, 1) :
    foldl(kron, fill(_BOND_I2, right_count))
  return kron(left, Matrix{ComplexF64}(op), right)
end


@testset "direct-center nonuniform and selected-operator paths" begin
  Random.seed!(20260911)
  sites = pauli_siteinds(5)
  base = random_mps(ComplexF64, sites; linkdims=[3, 7, 5, 2])
  identity = operator_product_state(sites, fill(_BOND_I2, 5))
  base = add(identity, 0.2base; maxdim=8, cutoff=0.0)
  h = spinhalf_xyz_bond_hamiltonian(Jx=0.8, Jy=-0.3, Jz=0.5)
  gate = pauli_gate_from_hamiltonian(h, -0.07)
  z = ComplexF64[1 0; 0 -1]
  exact = copy(base)
  tebd_evolve!(exact, gate, 2; maxdim=0, cutoff=0.0)
  expected = [pauli_trace(exact), pauli_expectation(exact, z, 2; normalize=false),
    pauli_expectation(exact, z, 3; normalize=false)]
  for direction in (:R, :L), cache_mode in (:none, :cache)
    actual = copy(base)
    cache = cache_mode === :cache ? MPSToolkit._DMTEnvCache(actual) : nothing
    dmt_step!(actual, gate, 2; maxdim=7, cutoff=1e-8, direction=direction,
      preserve_operators=[z], gate_backend=:direct, cache=cache)
    measured = [pauli_trace(actual), pauli_expectation(actual, z, 2; normalize=false),
      pauli_expectation(actual, z, 3; normalize=false)]
    @test norm(measured - expected) <= 1e-11 * max(norm(expected), 1)
    @test dim(linkind(actual, 2)) <= 7
    @test mapreduce(eltype, promote_type, ITensorMPS.data(actual)) === ComplexF64
  end

  qutrit_sites = operator_siteinds(4; d=3)
  qutrit = add(operator_basis_state(qutrit_sites, fill(1, 4)),
    0.1random_mps(qutrit_sites; linkdims=[4, 24, 5]); maxdim=24, cutoff=0.0)
  diagonal_h = Diagonal(collect(range(-1.0, 1.0; length=9)))
  qutrit_gate = operator_gate_from_hamiltonian(diagonal_h, 0.05; d=3)
  before = operator_trace(qutrit)
  dmt_step!(qutrit, qutrit_gate, 2; maxdim=19, cutoff=0.0, gate_backend=:direct)
  @test operator_trace(qutrit) ≈ before atol=1e-11 * max(abs(before), 1)
  @test dim(linkind(qutrit, 2)) <= 19
end

function _bond_plan_dense_step(plan, terms, nsites)
  step = Matrix{ComplexF64}(I, 2^nsites, 2^nsites)
  for entry in plan.entries
    gate = exp(-im * entry.duration * terms[entry.bond])
    step = _bond_embed(gate, entry.bond, nsites) * step
  end
  return step
end

function _bond_dense_coefficients(state)
  sites = [siteind(state, site) for site in 1:length(state)]
  return vec(ComplexF64[inner(pauli_basis_state(sites, collect(labels)), state)
    for labels in Iterators.product(ntuple(_ -> 1:4, length(state))...)])
end

function _bond_embed_superoperator(gate, bond, nsites)
  identity4 = Matrix{Float64}(I, 4, 4)
  left = bond == 1 ? Matrix{Float64}(I, 1, 1) : foldl(kron, fill(identity4, bond - 1))
  right_count = nsites - bond - 1
  right = right_count == 0 ? Matrix{Float64}(I, 1, 1) :
    foldl(kron, fill(identity4, right_count))
  return kron(left, gate, right)
end

@testset "generic nearest-neighbor DMT plans" begin
  for model in (:mfi, :xxz)
    nsites = 4
    terms = if model === :mfi
      [spinhalf_mixed_field_ising_bond_hamiltonian(nsites, bond;
        J=0.7, gx=1.1, gz=-0.3) for bond in 1:(nsites - 1)]
    else
      fill(spinhalf_xyz_bond_hamiltonian(Jx=1.0, Jy=1.0, Jz=0.6), nsites - 1)
    end
    h = sum(_bond_embed(terms[bond], bond, nsites) for bond in 1:(nsites - 1))
    total_time = 0.4
    for (order, threshold) in ((2, 3.5), (4, 10.0))
      errors = Float64[]
      for nsteps in (order == 2 ? (2, 4, 8) : (1, 2, 4))
        tau = total_time / nsteps
        plan = bond_dmt_plan(pauli_siteinds(nsites), terms, tau; order=order,
          maxdim=16, normalize=false)
        push!(errors, norm(_bond_plan_dense_step(plan, terms, nsites)^nsteps -
          exp(-im * total_time * h)))
        @test any(entry.duration < 0 for entry in plan.entries) == (order == 4)
        @test all(entry.direction == :R for entry in plan.entries)
      end
      @test errors[1] > errors[2] > errors[3]
      @test errors[1] / errors[2] > threshold
      @test errors[2] / errors[3] > threshold
    end
  end

  sites = pauli_siteinds(6)
  bulk = spinhalf_xyz_bond_hamiltonian(Jx=1.0, Jy=1.0, Jz=0.5)
  terms = fill(bulk, 5)
  defaults = bond_dmt_plan(sites, terms, 0.1; maxdim=16)
  @test defaults.order == 4
  @test defaults.composition == :merged7
  @test defaults.options.gate_backend == :qr
  @test defaults.options.truncation == :dense
  merged = bond_dmt_plan(sites, terms, 0.1; order=4, composition=:merged7,
    maxdim=16, direction=:R)
  unmerged = bond_dmt_plan(sites, terms, 0.1; order=4, composition=:yoshida9,
    maxdim=16, direction=:L)
  @test bond_dmt_attempted_updates(merged) == 4 * 3 + 3 * 2
  @test bond_dmt_attempted_updates(unmerged) == 6 * 3 + 3 * 2
  @test merged.signature != unmerged.signature
  @test all(entry.direction == :L for entry in unmerged.entries)
  # Homogeneous terms and repeated signed durations share the same compiled gate object.
  matching = findall(entry -> entry.duration == merged.entries[1].duration, merged.entries)
  @test length(matching) > 1
  @test all(merged.gates[index] === merged.gates[first(matching)] for index in matching)

  @test_throws ArgumentError bond_dmt_plan(sites, terms, 0.0)
  @test_throws ArgumentError bond_dmt_plan(sites, terms, 0.1; nstep=0)
  @test_throws ArgumentError bond_dmt_plan(sites, terms, 0.1; order=3)
  @test_throws ArgumentError bond_dmt_plan(sites, terms[1:4], 0.1)
  bad = copy(terms)
  bad[1] = copy(bad[1]); bad[1][1, 2] = 1im
  @test_throws ArgumentError bond_dmt_plan(sites, bad, 0.1)

  @test_throws ArgumentError bond_dmt_plan(sites, terms, 0.1; maxdim=1)
end

@testset "compiled bond-plan gates and executor match dense operator space" begin
  Random.seed!(20260910)
  nsites = 4
  sites = pauli_siteinds(nsites)
  terms = [spinhalf_mixed_field_ising_bond_hamiltonian(nsites, bond;
    J=0.7, gx=1.1, gz=-0.3) for bond in 1:(nsites - 1)]
  initial = random_mps(ComplexF64, sites; linkdims=4)
  normalize!(initial)
  for (order, composition, direction) in
      ((2, :merged7, :R), (4, :merged7, :R), (4, :yoshida9, :L))
    plan = bond_dmt_plan(sites, terms, 0.04; order=order, composition=composition,
      direction=direction, nstep=2, maxdim=16, normalize=false)
    dense = _bond_dense_coefficients(initial)
    for _ in 1:plan.nstep, entry in plan.entries
      gate = pauli_gate_from_hamiltonian(terms[entry.bond], entry.duration)
      dense = _bond_embed_superoperator(gate, entry.bond, nsites) * dense
    end
    evolved = copy(initial)
    dmt_evolve!(evolved, plan)
    actual = _bond_dense_coefficients(evolved)
    @test norm(actual - dense) <= 5e-13 * norm(dense)
  end
end

@testset "fused two-site gate and DMT" begin
  sites = pauli_siteinds(4)
  h = spinhalf_xyz_bond_hamiltonian(Jx=0.9, Jy=0.7, Jz=-0.2)
  gate = pauli_gate_from_hamiltonian(h, 0.08)
  base = add(operator_product_state(sites, fill(_BOND_I2, 4)),
    0.1 * random_mps(sites; linkdims=12); maxdim=12, cutoff=0.0)
  @test dim(linkind(base, 2)) > 10
  via_product = copy(base)
  via_fused = copy(base)
  dmt_step!(via_product, gate, 2; maxdim=10, cutoff=0.0,
    gate_backend=:product, direction=:R)
  dmt_step!(via_fused, gate, 2; maxdim=10, cutoff=0.0,
    gate_backend=:fused, direction=:R)
  @test maxlinkdim(via_fused) <= 12
  @test dim(linkind(via_fused, 2)) <= 10
  @test abs(inner(via_product, via_fused)) /
    (norm(via_product) * norm(via_fused)) > 1 - 1e-10
  @test norm(via_product) ≈ norm(via_fused) rtol=1e-10

  exact_product = copy(base)
  exact_fused = copy(base)
  dmt_step!(exact_product, gate, 2; maxdim=64, cutoff=0.0,
    gate_backend=:product, direction=:L)
  dmt_step!(exact_fused, gate, 2; maxdim=64, cutoff=0.0,
    gate_backend=:fused, direction=:L)
  @test norm(exact_product - exact_fused) <= 1e-10 * norm(exact_product)

  @test_throws ArgumentError dmt_step!(copy(base), Matrix{Float64}(I, 4, 4), 2;
    maxdim=10, gate_backend=:fused)
  fused_plan = bond_dmt_plan(sites, fill(h, 3), 0.05; order=2, maxdim=10,
    gate_backend=:fused, normalize=false)
  @test fused_plan.options.gate_backend == :fused
  evolved = copy(base)
  dmt_evolve!(evolved, fused_plan)
  @test dim(linkind(evolved, 2)) <= 10
  direct_plan = bond_dmt_plan(sites, fill(h, 3), 0.05; order=2, maxdim=10,
    gate_backend=:direct, normalize=false)
  @test direct_plan.options.gate_backend == :direct
  @test direct_plan.signature != fused_plan.signature
  direct_evolved = copy(base)
  dmt_evolve!(direct_evolved, direct_plan)
  @test maxlinkdim(direct_evolved) <= 10
  @test_throws ArgumentError bond_dmt_plan(sites, fill(h, 3), 0.05; order=2,
    maxdim=10, gate_backend=:direct, truncation=:random)

  @test MPSToolkit._operator_basis_operators(0, 2) == [Matrix{ComplexF64}(I, 1, 1)]
  for preserve_operators in (nothing, Matrix{ComplexF64}[]), direction in (:R, :L)
    radius_zero = copy(base)
    trace_before = pauli_trace(radius_zero)
    dmt_step!(radius_zero, gate, 2; maxdim=3, cutoff=0.0,
      preserve_diameter=1, preserve_operators=preserve_operators,
      gate_backend=:fused, direction=direction)
    @test dim(linkind(radius_zero, 2)) <= 3
    @test pauli_trace(radius_zero) ≈ trace_before atol=1e-10 * max(abs(trace_before), 1)
    radius_zero_direct = copy(base)
    dmt_step!(radius_zero_direct, gate, 2; maxdim=3, cutoff=0.0,
      preserve_diameter=1, preserve_operators=preserve_operators,
      gate_backend=:direct, direction=direction)
    @test dim(linkind(radius_zero_direct, 2)) <= 3
    @test pauli_trace(radius_zero_direct) ≈ trace_before atol=1e-10 * max(abs(trace_before), 1)
  end

  probes = diameter_probes(4, 2, 3; full_basis_width=3, nrandom=1,
    rng=MersenneTwister(20260910))
  exact = copy(base)
  tebd_evolve!(exact, gate, 2; maxdim=0, cutoff=0.0)
  expected = operator_expectation_profile(exact, probes; normalize=false)
  covered = [index for (index, (start, op)) in enumerate(probes)
    if guarantee_covers(start, probe_span(op, 2), 2, 1)]
  for backend in (:product, :qr, :fused, :direct), direction in (:R, :L)
    truncated = copy(base)
    dmt_step!(truncated, gate, 2; maxdim=10, cutoff=0.0,
      gate_backend=backend, direction=direction)
    @test dim(linkind(truncated, 2)) <= 10
    actual = operator_expectation_profile(truncated, probes; normalize=false)
    @test preservation_error(expected[covered], actual[covered]) < 1e-11
  end

  low_rank_z = ComplexF64[1 0; 0 -1]
  low_rank = operator_product_state(sites, [_BOND_I2, low_rank_z, _BOND_I2, _BOND_I2])
  low_rank_exact = copy(low_rank)
  tebd_evolve!(low_rank_exact, gate, 2; maxdim=0, cutoff=0.0)
  low_rank_direct = copy(low_rank)
  dmt_step!(low_rank_direct, gate, 2; maxdim=10, cutoff=0.0, gate_backend=:direct)
  @test _bond_dense_coefficients(low_rank_direct) ≈
    _bond_dense_coefficients(low_rank_exact) atol=1e-11
end
