using Base.Docs

@static if isdefined(Docs, :hasdoc)
  _has_binding_doc(mod::Module, name::Symbol) = Docs.hasdoc(mod, name)
else
  # `Docs.hasdoc` arrived in Julia 1.11, and Project.toml declares a 1.10 floor. Reproduce its
  # lookup: scan the doc metadata of every module that registered any, then follow the binding's
  # alias so a name re-exported from a submodule resolves to the definition carrying the
  # docstring. Verified against 1.10.11 and 1.12.6 to agree on documented bindings, undocumented
  # ones, the re-export path, and names that do not exist.
  _has_binding_doc(mod::Module, name::Symbol) = _has_binding_doc(Docs.Binding(mod, name))
  function _has_binding_doc(binding::Docs.Binding)
    for m in Docs.modules
      dict = Docs.meta(m; autoinit=false)
      if !isnothing(dict) && haskey(dict, binding)
        return true
      end
    end
    alias = Docs.aliasof(binding)
    return alias == binding ? false : _has_binding_doc(alias)
  end
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
    :operator_basis_matrices,
    :local_dimension,
    :pauli_basis_state,
    :pauli_total_sz_state,
    :pauli_domain_wall_state,
    :operator_siteinds,
    :operator_basis_state,
    :operator_product_state,
    :operator_local_sum_state,
    :pauli_gate,
    :pauli_gate_from_hamiltonian,
    :pauli_lindblad_generator,
    :pauli_gate_from_lindbladian,
    :operator_gate,
    :operator_gate_from_hamiltonian,
    :operator_gate_from_imaginary_time,
    :operator_lindblad_generator,
    :operator_gate_from_lindbladian,
    :spinhalf_matrices,
    :spinhalf_xyz_bond_hamiltonian,
    :spinhalf_tfim_bond_hamiltonian,
    :pxp_term_hamiltonian,
    :pxp_term_support,
    :pxp_constraint_mpo,
    :pauli_state_from_mpo,
    :pauli_superoperator_mpo,
    :operator_state_from_mpo,
    :operator_superoperator_mpo,
    :pauli_pxp_constraint_state,
    :pauli_pxp_constraint_projector,
    :pauli_trace,
    :pauli_expectation,
    :pauli_expectation_profile,
    :operator_trace,
    :operator_expectation,
    :operator_expectation_profile,
    :pauli_gate_from_imaginary_time,
    :pauli_gibbs_state,
    :operator_gibbs_state,
    Symbol("constrained_dmt_evolve!"),
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
    :pauli_fdaoe_projector,
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
