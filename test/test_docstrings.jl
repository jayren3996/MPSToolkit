using Base.Docs

function _has_binding_doc(mod::Module, name::Symbol)
  return Docs.hasdoc(mod, name)
end

@testset "docstring coverage" begin
  documented_bindings = [
    :Evolution,
    :Observables,
    :ScarFinder,
    :Bases,
    :OperatorSpace,
    :Models,
    :Chebyshev,
    :LocalGateEvolution,
    :DMTGateEvolution,
    :TDVPEvolution,
    :BondDimTruncation,
    :EnergyTarget,
    :SelectionContext,
    :EntropySelector,
    :FidelitySelector,
    Symbol("evolve!"),
    Symbol("project!"),
    :energy_density,
    :bond_entropy,
    :entanglement_spectrum,
    :fidelity_distance,
    Symbol("scarfinder_step!"),
    Symbol("scarfinder!"),
    Symbol("tebd_evolve!"),
    Symbol("dmt_step!"),
    Symbol("dmt_evolve!"),
    :local_gates_from_hamiltonians,
    :tebd_evolution_from_hamiltonians,
    :tebd_strang_schedule,
    :tebd_strang_evolution,
    Symbol("tdvp_evolve!"),
    Symbol("_assign_state!"),
    :_clone_state,
    Symbol("match_energy!"),
    Symbol("trajectory_refine!"),
    :_diagonal_values,
    :_entropy_from_values,
    :score,
    :_bond_start,
    :_dense_local_operator,
    :_operator_span,
    :_finite_local_wavefunction,
    :pauli_matrices,
    :pauli_basis,
    :pauli_components,
    :pauli_basis_state,
    :pauli_total_sz_state,
    :pauli_gate,
    :pauli_gate_from_hamiltonian,
    :pauli_lindblad_generator,
    :pauli_gate_from_lindbladian,
    :spinhalf_matrices,
    :spinhalf_xyz_bond_hamiltonian,
    :spinhalf_tfim_bond_hamiltonian,
    :ChebyshevRescaling,
    :chebyshev_rescaling,
    :rescale_hamiltonian,
    :SpectralFunction,
    :chebyshev_moments,
    Symbol("energy_cutoff!"),
    :jackson_damping,
    :jackson_kernel,
    :reconstruct_chebyshev,
    :spectral_function,
    :pauli_daoe_projector,
    :fdaoe_projector,
    :pauli_siteinds,
  ]

  for name in documented_bindings
    @test _has_binding_doc(MPSToolkit, name)
  end
end

@testset "docstring examples evaluate correctly" begin
  # Mirrors the `jldoctest` examples in the docstrings, verified here as robust value checks
  # (independent of REPL display formatting) so they are exercised by the normal test suite.
  # pauli_matrices docstring example
  @test pauli_matrices().Z == ComplexF64[1 0; 0 -1]
  @test keys(pauli_matrices(; include_identity=false)) == (:X, :Y, :Z)
  # pauli_basis_state docstring example
  sites = pauli_siteinds(3)
  rho = pauli_basis_state(sites, [:I, :Z, :I])
  @test length(rho) == 3
end
