"""
MPSToolkit provides finite-MPS utilities for local-gate and TDVP evolution, observables,
operator-space workflows, ScarFinder refinement, and Chebyshev spectral reconstruction.
"""
module MPSToolkit

using LinearAlgebra
using Printf
using Serialization
using ITensors
using ITensorMPS

include("evolution/types.jl")
include("evolution/dispatch.jl")
include("evolution/tebd.jl")
include("evolution/tdvp.jl")
include("observables/dispatch.jl")
include("observables/entanglement.jl")
include("observables/energy.jl")
include("scarfinder/types.jl")
include("scarfinder/dispatch.jl")
include("scarfinder/selectors.jl")
include("scarfinder/algorithm.jl")
include("bases/pauli.jl")
include("bases/operator_basis.jl")
include("models/spinhalf.jl")
include("models/pxp.jl")
include("operator_space/states.jl")
include("operator_space/gates.jl")
include("operator_space/helpers.jl")
include("operator_space/vectorize.jl")
include("operator_space/pxp.jl")
include("operator_space/expectations.jl")
include("operator_space/thermal.jl")
include("operator_space/checkpoint.jl")
include("operator_space/dmt/lowrank.jl")
include("operator_space/dmt.jl")
# After dmt.jl: `_dmt_bond_truncate!` annotates its `cache` keyword with `_DMTEnvCache`, and a
# type annotation in a method signature is resolved when the method is defined.
include("operator_space/dmt/bond.jl")
include("operator_space/daoe.jl")
include("operator_space/constrained.jl")
include("chebyshev/types.jl")
include("chebyshev/reconstruction.jl")
include("chebyshev/moments.jl")

"""
Namespace for evolution configuration objects and TEBD/TDVP drivers.
"""
module Evolution
using ..MPSToolkit:
  LocalGateEvolution,
  DMTGateEvolution,
  TDVPEvolution,
  evolve!,
  tebd_evolve!,
  dmt_evolve!,
  tdvp_evolve!,
  local_gates_from_hamiltonians,
  tebd_evolution_from_hamiltonians,
  tebd_strang_schedule,
  tebd_strang_evolution
export LocalGateEvolution, DMTGateEvolution, TDVPEvolution, evolve!, tebd_evolve!, dmt_evolve!, tdvp_evolve!, local_gates_from_hamiltonians, tebd_evolution_from_hamiltonians, tebd_strang_schedule, tebd_strang_evolution
end

"""
Namespace for finite-MPS energy, entanglement, and fidelity observables.
"""
module Observables
using ..MPSToolkit: energy_density, bond_entropy, entanglement_spectrum, fidelity_distance
export energy_density, bond_entropy, entanglement_spectrum, fidelity_distance
end

"""
Namespace for ScarFinder projection, energy matching, and selector APIs.
"""
module ScarFinder
using ..MPSToolkit:
  BondDimTruncation,
  EnergyTarget,
  SelectionContext,
  EntropySelector,
  FidelitySelector,
  project!,
  match_energy!,
  trajectory_refine!,
  scarfinder_step!,
  scarfinder!
export BondDimTruncation,
  EnergyTarget,
  SelectionContext,
  EntropySelector,
  FidelitySelector,
  project!,
  match_energy!,
  trajectory_refine!,
  scarfinder_step!,
  scarfinder!
end

"""
Namespace for local basis helpers.
"""
module Bases
using ..MPSToolkit: pauli_matrices, pauli_basis, pauli_components, operator_basis_matrices, local_dimension
export pauli_matrices, pauli_basis, pauli_components, operator_basis_matrices, local_dimension
end

"""
Namespace for Pauli-basis operator-space states, gates, DMT, and DAOE projectors.
"""
module OperatorSpace
using ..MPSToolkit:
  pauli_siteinds,
  pauli_basis_state,
  pauli_total_sz_state,
  pauli_domain_wall_state,
  operator_siteinds,
  operator_basis_state,
  operator_product_state,
  operator_local_sum_state,
  pauli_gate,
  pauli_gate_from_hamiltonian,
  pauli_lindblad_generator,
  pauli_gate_from_lindbladian,
  operator_gate,
  operator_gate_from_hamiltonian,
  operator_gate_from_imaginary_time,
  operator_lindblad_generator,
  operator_gate_from_lindbladian,
  pauli_state_from_mpo,
  pauli_superoperator_mpo,
  operator_state_from_mpo,
  operator_superoperator_mpo,
  pauli_pxp_constraint_state,
  pauli_pxp_constraint_projector,
  pauli_trace,
  pauli_expectation,
  pauli_expectation_profile,
  operator_trace,
  operator_expectation,
  operator_expectation_profile,
  pauli_gate_from_imaginary_time,
  pauli_gibbs_state,
  operator_gibbs_state,
  constrained_dmt_evolve!,
  DMTOptions,
  dmt_step!,
  dmt_evolve!,
  DMTGateEvolution,
  pauli_daoe_projector,
  pauli_fdaoe_projector,
  fdaoe_projector
export pauli_siteinds,
  pauli_basis_state,
  pauli_total_sz_state,
  pauli_domain_wall_state,
  operator_siteinds,
  operator_basis_state,
  operator_product_state,
  operator_local_sum_state,
  pauli_gate,
  pauli_gate_from_hamiltonian,
  pauli_lindblad_generator,
  pauli_gate_from_lindbladian,
  operator_gate,
  operator_gate_from_hamiltonian,
  operator_gate_from_imaginary_time,
  operator_lindblad_generator,
  operator_gate_from_lindbladian,
  pauli_state_from_mpo,
  pauli_superoperator_mpo,
  operator_state_from_mpo,
  operator_superoperator_mpo,
  pauli_pxp_constraint_state,
  pauli_pxp_constraint_projector,
  pauli_trace,
  pauli_expectation,
  pauli_expectation_profile,
  operator_trace,
  operator_expectation,
  operator_expectation_profile,
  pauli_gate_from_imaginary_time,
  pauli_gibbs_state,
  operator_gibbs_state,
  constrained_dmt_evolve!,
  DMTOptions,
  dmt_step!,
  dmt_evolve!,
  DMTGateEvolution,
  pauli_daoe_projector,
  pauli_fdaoe_projector,
  fdaoe_projector
end

"""
Namespace for dense spin-half model-building helpers.
"""
module Models
using ..MPSToolkit:
  spinhalf_matrices,
  spinhalf_xyz_bond_hamiltonian,
  spinhalf_tfim_bond_hamiltonian,
  pxp_term_hamiltonian,
  pxp_term_support,
  pxp_constraint_mpo
export spinhalf_matrices,
  spinhalf_xyz_bond_hamiltonian,
  spinhalf_tfim_bond_hamiltonian,
  pxp_term_hamiltonian,
  pxp_term_support,
  pxp_constraint_mpo
end

"""
Namespace for Chebyshev moments, kernels, reconstruction, and energy-window cutoff tools.
"""
module Chebyshev
using ..MPSToolkit:
  ChebyshevRescaling,
  chebyshev_rescaling,
  rescale_hamiltonian,
  SpectralFunction,
  chebyshev_moments,
  energy_cutoff!,
  jackson_damping,
  jackson_kernel,
  reconstruct_chebyshev,
  spectral_function
export ChebyshevRescaling,
  chebyshev_rescaling,
  rescale_hamiltonian,
  SpectralFunction,
  chebyshev_moments,
  energy_cutoff!,
  jackson_damping,
  jackson_kernel,
  reconstruct_chebyshev,
  spectral_function
end

export Evolution, Observables, ScarFinder, Bases, OperatorSpace, Models, Chebyshev
export evolve!, project!, energy_density, bond_entropy, entanglement_spectrum, fidelity_distance
export scarfinder_step!, scarfinder!
export LocalGateEvolution, DMTGateEvolution, TDVPEvolution, BondDimTruncation, EnergyTarget, SelectionContext, EntropySelector, FidelitySelector
export tebd_evolve!, dmt_evolve!, tdvp_evolve!, local_gates_from_hamiltonians, tebd_evolution_from_hamiltonians, tebd_strang_schedule, tebd_strang_evolution
export pauli_matrices, pauli_basis, pauli_components, operator_basis_matrices, local_dimension
export spinhalf_matrices, spinhalf_xyz_bond_hamiltonian, spinhalf_tfim_bond_hamiltonian
export pxp_term_hamiltonian, pxp_term_support, pxp_constraint_mpo
export pauli_siteinds, pauli_basis_state, pauli_total_sz_state, pauli_domain_wall_state, operator_siteinds, operator_basis_state, operator_product_state, operator_local_sum_state, pauli_gate, pauli_gate_from_hamiltonian, pauli_lindblad_generator, pauli_gate_from_lindbladian, operator_gate, operator_gate_from_hamiltonian, operator_gate_from_imaginary_time, operator_lindblad_generator, operator_gate_from_lindbladian, DMTOptions, dmt_step!, dmt_evolve!, pauli_daoe_projector, pauli_fdaoe_projector, fdaoe_projector
export pauli_state_from_mpo, pauli_superoperator_mpo, operator_state_from_mpo, operator_superoperator_mpo, pauli_pxp_constraint_state, pauli_pxp_constraint_projector, pauli_trace, pauli_expectation, pauli_expectation_profile, operator_trace, operator_expectation, operator_expectation_profile
export pauli_gate_from_imaginary_time, pauli_gibbs_state, operator_gibbs_state, constrained_dmt_evolve!
export DMTCheckpoint, dmt_checkpoint_path, dmt_checkpoint_save, dmt_checkpoint_load, dmt_checkpoint_latest, dmt_checkpoint_resume
export ChebyshevRescaling, chebyshev_rescaling, rescale_hamiltonian, SpectralFunction, chebyshev_moments, energy_cutoff!, jackson_damping, jackson_kernel, reconstruct_chebyshev, spectral_function

end
