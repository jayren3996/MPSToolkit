@testset "configuration types" begin
  evo = LocalGateEvolution(rand(ComplexF64, 4, 4), 0.01; nstep=10)
  dmt_evo = DMTGateEvolution(rand(ComplexF64, 4, 4), 0.01; schedule=[1], maxdim=6)
  tdvp = TDVPEvolution(rand(2, 2), -0.1im; time_step=-0.05im, nsteps=2, nsweeps=2, reverse_step=false, updater_backend="exponentiate", normalize=false, solver_kwargs=(cutoff=1e-9,))
  trunc = BondDimTruncation(8; cutoff=1e-10)
  target = EnergyTarget(0.0; operator=rand(4, 4), tol=1e-6, alpha=0.1, maxstep=50)
  selector = EntropySelector()
  fidelity_selector = FidelitySelector()
  selection_context = SelectionContext()

  @test evo.dt == 0.01
  @test evo.nstep == 10
  @test dmt_evo.dt == 0.01
  @test dmt_evo.maxdim == 6
  @test dmt_evo.reverse_schedule == [1]
  @test tdvp.t == -0.1im
  @test tdvp.time_step == -0.05im
  @test tdvp.nsteps == 2
  @test tdvp.nsweeps == 2
  @test tdvp.reverse_step == false
  @test tdvp.updater_backend == "exponentiate"
  @test tdvp.normalize == false
  @test tdvp.solver_kwargs.cutoff == 1e-9
  @test trunc.maxdim == 8
  @test trunc.cutoff == 1e-10
  @test target.target == 0.0
  @test size(target.operator) == (4, 4)
  @test target.tol == 1e-6
  @test selector isa EntropySelector
  @test fidelity_selector isa FidelitySelector
  @test isnothing(selection_context.reference_state)
end

mutable struct DummyState
  value::Int
end

mutable struct TraceState
  trace::Vector{Symbol}
end

mutable struct RefineState
  value::Int
end

mutable struct FidelityState
  value::Int
end

mutable struct StepCountState
  value::Int
end

MPSToolkit.evolve!(psi::DummyState, evo::LocalGateEvolution) = (psi.value += 1; psi)
MPSToolkit.project!(psi::DummyState, trunc::BondDimTruncation) = (psi.value *= 2; psi)
MPSToolkit.energy_density(psi::DummyState, op; kwargs...) = psi.value

MPSToolkit.evolve!(psi::TraceState, evo::LocalGateEvolution) = (push!(psi.trace, :evolve); psi)
MPSToolkit.project!(psi::TraceState, trunc::BondDimTruncation) = (push!(psi.trace, :project); psi)
MPSToolkit.match_energy!(psi::TraceState, evo::LocalGateEvolution, trunc::BondDimTruncation, target::EnergyTarget) = (push!(psi.trace, :target); psi)
MPSToolkit.trajectory_refine!(psi::TraceState, evo::LocalGateEvolution, trunc::BondDimTruncation, selector; kwargs...) = (push!(psi.trace, :refine); psi)

MPSToolkit.evolve!(psi::RefineState, evo::LocalGateEvolution) = (psi.value += 1; psi)
MPSToolkit.project!(psi::RefineState, trunc::BondDimTruncation) = (psi.value *= 2; psi)
MPSToolkit.match_energy!(psi::RefineState, evo::LocalGateEvolution, trunc::BondDimTruncation, target::EnergyTarget) = (psi.value += 100; psi)
MPSToolkit.score(selector::EntropySelector, psi::RefineState, context::SelectionContext) = -psi.value
MPSToolkit._assign_state!(psi::RefineState, updated::RefineState) = (psi.value = updated.value; psi)

MPSToolkit.evolve!(psi::FidelityState, evo::LocalGateEvolution) = (psi.value += 1; psi)
MPSToolkit.project!(psi::FidelityState, trunc::BondDimTruncation) = psi
MPSToolkit.score(selector::FidelitySelector, psi::FidelityState, context::SelectionContext) = abs(psi.value - context.reference_state)
MPSToolkit._assign_state!(psi::FidelityState, updated::FidelityState) = (psi.value = updated.value; psi)

MPSToolkit.evolve!(psi::StepCountState, evo::LocalGateEvolution) = (psi.value += evo.nstep; psi)
MPSToolkit.evolve!(psi::StepCountState, evo::DMTGateEvolution) = (psi.value += evo.nstep; psi)
MPSToolkit.evolve!(psi::StepCountState, evo::TDVPEvolution) = (psi.value += something(evo.nsteps, evo.nsweeps); psi)
MPSToolkit.project!(psi::StepCountState, trunc::BondDimTruncation) = psi

@testset "feature namespaces" begin
  @test MPSToolkit.Evolution.LocalGateEvolution === LocalGateEvolution
  @test MPSToolkit.Evolution.DMTGateEvolution === DMTGateEvolution
  @test MPSToolkit.ScarFinder.BondDimTruncation === BondDimTruncation
  @test MPSToolkit.Observables.bond_entropy === bond_entropy
  @test MPSToolkit.Bases.pauli_basis === pauli_basis
  @test MPSToolkit.Evolution.local_gates_from_hamiltonians === local_gates_from_hamiltonians
  @test MPSToolkit.Evolution.tebd_evolution_from_hamiltonians === tebd_evolution_from_hamiltonians
  @test MPSToolkit.Evolution.tebd_strang_schedule === tebd_strang_schedule
  @test MPSToolkit.Evolution.tebd_strang_evolution === tebd_strang_evolution
  @test MPSToolkit.Models.spinhalf_tfim_bond_hamiltonian === spinhalf_tfim_bond_hamiltonian
  @test MPSToolkit.OperatorSpace.pauli_daoe_projector === pauli_daoe_projector
  @test MPSToolkit.OperatorSpace.fdaoe_projector === fdaoe_projector
  @test MPSToolkit.OperatorSpace.pauli_basis_state === pauli_basis_state
  @test MPSToolkit.OperatorSpace.pauli_total_sz_state === pauli_total_sz_state
  @test MPSToolkit.OperatorSpace.pauli_lindblad_generator === pauli_lindblad_generator
  @test MPSToolkit.OperatorSpace.pauli_gate_from_lindbladian === pauli_gate_from_lindbladian
  @test MPSToolkit.OperatorSpace.dmt_evolve! === dmt_evolve!
  @test MPSToolkit.OperatorSpace.dmt_step! === dmt_step!
  @test MPSToolkit.Chebyshev.chebyshev_moments === chebyshev_moments
  @test MPSToolkit.Chebyshev.spectral_function === spectral_function
end

@testset "pauli basis helpers" begin
  basis = pauli_basis()
  coeffs = pauli_components(ComplexF64[1 0; 0 0])

  @test length(basis) == 4
  @test basis[2].first == :X
  @test coeffs.I ≈ 0.5
  @test coeffs.Z ≈ 0.5
  @test coeffs.X ≈ 0.0
  @test coeffs.Y ≈ 0.0
end

@testset "dispatch points" begin
  psi = DummyState(1)
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1)
  trunc = BondDimTruncation(2)
  evolve!(psi, evo)
  project!(psi, trunc)
  @test psi.value == 4
  @test energy_density(psi, nothing) == 4
end

@testset "scarfinder step order" begin
  psi = TraceState(Symbol[])
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1)
  trunc = BondDimTruncation(2)
  target = EnergyTarget(0.0; operator=reshape(1:4, 2, 2))
  scarfinder_step!(psi, evo, trunc; target_energy=target)
  @test psi.trace == [:evolve, :project, :target]
end

@testset "scarfinder loop refinement" begin
  psi = TraceState(Symbol[])
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1)
  trunc = BondDimTruncation(2)
  scarfinder!(psi, evo, trunc, 2; refine=true, selector=EntropySelector())
  @test psi.trace == [:evolve, :project, :evolve, :project, :refine]
end

@testset "refinement preserves target energy" begin
  psi = RefineState(0)
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1)
  trunc = BondDimTruncation(2)
  target = EnergyTarget(0.0; operator=reshape(1:4, 2, 2))
  scarfinder!(psi, evo, trunc, 1; refine=true, selector=EntropySelector(), target_energy=target, refine_steps=1)
  @test psi.value == 306
end

@testset "fidelity selector uses selection context" begin
  psi = FidelityState(0)
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1)
  trunc = BondDimTruncation(2)
  selector = FidelitySelector()
  context = SelectionContext(; reference_state=3)

  scarfinder!(psi, evo, trunc, 1; refine=true, selector=selector, selector_context=context, refine_steps=3)

  @test psi.value == 3
end

@testset "scarfinder upgrades single-step evolutions" begin
  trunc = BondDimTruncation(2)

  local_state = StepCountState(0)
  local_evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1)
  @test_logs (:warn, r"ScarFinder") scarfinder_step!(local_state, local_evo, trunc)
  @test local_state.value == 10

  dmt_state = StepCountState(0)
  dmt_evo = DMTGateEvolution(reshape(1:4, 2, 2), 0.1; schedule=[1])
  @test_logs (:warn, r"ScarFinder") scarfinder_step!(dmt_state, dmt_evo, trunc)
  @test dmt_state.value == 10

  tdvp_state = StepCountState(0)
  tdvp_evo = TDVPEvolution(reshape(1:4, 2, 2), -0.1im; nsteps=1)
  @test_logs (:warn, r"ScarFinder") scarfinder_step!(tdvp_state, tdvp_evo, trunc)
  @test tdvp_state.value == 10
end

@testset "scarfinder preserves non-unit step counts" begin
  psi = StepCountState(0)
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1; nstep=3)
  trunc = BondDimTruncation(2)

  scarfinder_step!(psi, evo, trunc)

  @test psi.value == 3
end
