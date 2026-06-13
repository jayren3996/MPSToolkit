using ITensors
using ITensorMPS

@testset "configuration types" begin
  evo = LocalGateEvolution(rand(ComplexF64, 4, 4), 0.01; nstep=10)
  dmt_evo = DMTGateEvolution(rand(ComplexF64, 4, 4), 0.01; schedule=[1], maxdim=6, connector_buffer=6)
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

@testset "configuration validation" begin
  gate = ComplexF64[1 0; 0 1]

  @test_throws ArgumentError LocalGateEvolution(gate, 0.1; schedule=[1], nstep=-1)
  @test_throws ArgumentError LocalGateEvolution(gate, 0.1; schedule=[1], maxdim=-1)
  @test_throws ArgumentError DMTGateEvolution(gate, 0.1; schedule=[1], maxdim=0)
  @test_throws ArgumentError DMTGateEvolution(gate, 0.1; schedule=[1], connector_buffer=-1)
  @test_throws ArgumentError TDVPEvolution(gate, 0.1; nsteps=0)
  @test_throws ArgumentError TDVPEvolution(gate, 0.1; nsteps=1, solver_kwargs=(reverse_step=false,))
  @test_throws ArgumentError BondDimTruncation(0)
  @test_throws ArgumentError BondDimTruncation(2; cutoff=-1e-3)
  @test_throws ArgumentError EnergyTarget(0.0; maxstep=0)
  @test_throws ArgumentError EnergyTarget(0.0; tol=-1e-3)
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

mutable struct NaNScoreState
  value::Int
end

MPSToolkit.evolve!(psi::NaNScoreState, evo::LocalGateEvolution) = (psi.value += 1; psi)
MPSToolkit.project!(psi::NaNScoreState, trunc::BondDimTruncation) = psi
MPSToolkit.score(selector::EntropySelector, psi::NaNScoreState, context::SelectionContext) = NaN
MPSToolkit._assign_state!(psi::NaNScoreState, updated::NaNScoreState) = (psi.value = updated.value; psi)

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
  @test evolve!(psi, evo) === psi
  @test project!(psi, trunc) === psi
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

@testset "refinement warns on NaN selector scores" begin
  psi = NaNScoreState(0)
  # nstep=10 avoids the unrelated single-step normalization warning, isolating the NaN warning.
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1; nstep=10)
  trunc = BondDimTruncation(2)
  @test_logs (:warn, r"NaN") MPSToolkit.trajectory_refine!(psi, evo, trunc, EntropySelector(); refine_steps=2)
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

@testset "scarfinder single-step upgrade preserves evolution fields" begin
  # A DMTGateEvolution built to track a traceless operator (normalize=false) keeps the
  # default nstep=1, so ScarFinder's single-step upgrade rebuilds it. The rebuild must carry
  # the normalize field across; otherwise the upgraded evolution silently re-normalizes the
  # trajectory and destroys the unnormalized-trace diagnostic the field exists to preserve.
  dmt_evo = DMTGateEvolution(rand(ComplexF64, 4, 4), 0.01; schedule=[1], maxdim=6,
                             connector_buffer=6, normalize=false)
  dmt_upgraded = MPSToolkit._scarfinder_evolution(dmt_evo; warn=false)
  @test dmt_upgraded.nstep == 10
  @test dmt_upgraded.normalize == false

  # The TDVP rebuild already preserves normalize; pin it so the two drivers stay consistent.
  tdvp_evo = TDVPEvolution(reshape(1:4, 2, 2), -0.1im; nsteps=1, normalize=false)
  tdvp_upgraded = MPSToolkit._scarfinder_evolution(tdvp_evo; warn=false)
  @test tdvp_upgraded.normalize == false
end

@testset "scarfinder preserves non-unit step counts" begin
  psi = StepCountState(0)
  evo = LocalGateEvolution(reshape(1:4, 2, 2), 0.1; nstep=3)
  trunc = BondDimTruncation(2)

  scarfinder_step!(psi, evo, trunc)

  @test psi.value == 3
end

@testset "energy matching uses normalized energy" begin
  sites = siteinds("S=1/2", 1)
  psi = MPS(sites, n -> "Up")
  initial_norm = real(inner(psi, psi))
  h = ComplexF64[1 0; 0 0]
  evo = LocalGateEvolution(h, 0.1; schedule=[1], nstep=2, maxdim=4, cutoff=1e-12)
  trunc = BondDimTruncation(4; cutoff=1e-12)
  target = EnergyTarget(0.1; operator=h, alpha=1.0, maxstep=1, tol=1e-12)

  MPSToolkit.match_energy!(psi, evo, trunc, target)

  final_norm = real(inner(psi, psi))
  normalized_energy = energy_density(psi, h) / final_norm
  @test normalized_energy ≈ 1.0 atol = 1e-10
  @test final_norm ≈ initial_norm atol = 1e-10
end

@testset "MPO energy matching moves energy toward target" begin
  sites = siteinds("S=1/2", 4)
  psi = MPS(sites, n -> isodd(n) ? "Up" : "Dn")
  os = OpSum()
  for j in 1:3
    os += "Sz", j, "Sz", j + 1
  end
  for j in 1:4
    os += 0.9, "Sx", j
  end
  H = MPO(os, sites)
  evo = TDVPEvolution(H, -0.05im; time_step=-0.05im, nsteps=1,
                      solver_kwargs=(maxdim=16, cutoff=1e-12))
  trunc = BondDimTruncation(16; cutoff=1e-12)

  e_initial = energy_density(psi, H)
  target_value = e_initial - 0.1
  target = EnergyTarget(target_value; operator=H, alpha=1.0, maxstep=40, tol=1e-8)

  MPSToolkit.match_energy!(psi, evo, trunc, target)
  e_final = energy_density(psi, H)

  # The MPO correction must perform imaginary-time evolution that actually moves
  # the energy toward the target. A purely-imaginary tdvp time argument would do
  # unitary (real-time) evolution instead, conserving <H> and freezing the energy.
  @test abs(e_final - target_value) < abs(e_initial - target_value) - 1e-6
end

@testset "MPO energy matching recovers from overshoot" begin
  sites = siteinds("S=1/2", 4)
  psi = MPS(sites, n -> isodd(n) ? "Up" : "Dn")
  os = OpSum()
  for j in 1:3
    os += "Sz", j, "Sz", j + 1
  end
  for j in 1:4
    os += 0.9, "Sx", j
  end
  H = MPO(os, sites)
  evo = TDVPEvolution(H, -0.05im; time_step=-0.05im, nsteps=1,
                      solver_kwargs=(maxdim=16, cutoff=1e-12))
  trunc = BondDimTruncation(16; cutoff=1e-12)
  e_initial = energy_density(psi, H)
  target_value = e_initial - 0.02
  # alpha is large enough that the first cooling step overshoots the small gap; the
  # proportional rollback must land near the target rather than bouncing back to e_initial.
  target = EnergyTarget(target_value; operator=H, alpha=8.0, maxstep=12, tol=1e-9)
  MPSToolkit.match_energy!(psi, evo, trunc, target)
  e_final = energy_density(psi, H)
  @test abs(e_final - target_value) < 0.5 * abs(e_initial - target_value)
end

@testset "dense energy matching does not stall on small steps" begin
  X = ComplexF64[0 1; 1 0]
  Id = ComplexF64[1 0; 0 1]
  h = -(kron(X, Id) + kron(Id, X))
  sites = siteinds("S=1/2", 2)
  psi = MPS(sites, n -> "Up")
  evo = LocalGateEvolution(exp(0.0 * h), 0.1; schedule=[1], nstep=1, maxdim=16, cutoff=1e-12)
  trunc = BondDimTruncation(16; cutoff=1e-12)
  e_initial = energy_density(psi, h)
  target_value = e_initial - 0.3
  # alpha is small enough that one substep improves the error by less than tol; the loop
  # must accumulate progress over many steps instead of rolling back the first slow step
  # and returning the state unchanged.
  target = EnergyTarget(target_value; operator=h, alpha=0.05, maxstep=200, tol=0.05)
  MPSToolkit.match_energy!(psi, evo, trunc, target)
  e_final = energy_density(psi, h)
  @test abs(e_final - target_value) < 0.5 * abs(e_initial - target_value)
end

@testset "dense energy matching is schedule-duplication independent" begin
  X = ComplexF64[0 1; 1 0]
  Id = ComplexF64[1 0; 0 1]
  h = -(kron(X, Id) + kron(Id, X)) / 2
  sites = siteinds("S=1/2", 3)
  make_state() = MPS(sites, n -> "Up")
  trunc = BondDimTruncation(16; cutoff=1e-12)
  target_value = energy_density(make_state(), h) - 0.5
  mk_target() = EnergyTarget(target_value; operator=h, alpha=0.3, maxstep=1, tol=1e-12)
  gate = exp(0.0 * h)
  # One correction substep must move the energy by an amount independent of how many times
  # the schedule happens to list each bond (single pass vs. duplicated pass).
  evo_single = LocalGateEvolution(gate, 0.1; schedule=[1, 2], nstep=1, maxdim=16, cutoff=1e-12)
  evo_dup = LocalGateEvolution(gate, 0.1; schedule=[1, 2, 1, 2], nstep=1, maxdim=16, cutoff=1e-12)
  psi_single = make_state()
  MPSToolkit.match_energy!(psi_single, evo_single, trunc, mk_target())
  psi_dup = make_state()
  MPSToolkit.match_energy!(psi_dup, evo_dup, trunc, mk_target())
  @test energy_density(psi_single, h) ≈ energy_density(psi_dup, h) atol = 1e-9
end
