"""
MPSToolkit provides finite-MPS utilities for local-gate and TDVP evolution, observables,
operator-space workflows, ScarFinder refinement, and Chebyshev spectral reconstruction.
"""
module MPSToolkit

using LinearAlgebra
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
include("models/spinhalf.jl")
include("operator_space/helpers.jl")
include("operator_space/dmt.jl")
include("operator_space/daoe.jl")
include("chebyshev/reconstruction.jl")
include("chebyshev/types.jl")
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
using ..MPSToolkit: pauli_matrices, pauli_basis, pauli_components
export pauli_matrices, pauli_basis, pauli_components
end

"""
Namespace for Pauli-basis operator-space states, gates, DMT, and DAOE projectors.
"""
module OperatorSpace
using ..MPSToolkit:
  pauli_siteinds,
  pauli_basis_state,
  pauli_total_sz_state,
  pauli_gate,
  pauli_gate_from_hamiltonian,
  pauli_lindblad_generator,
  pauli_gate_from_lindbladian,
  DMTOptions,
  dmt_step!,
  dmt_evolve!,
  DMTGateEvolution,
  pauli_daoe_projector,
  fdaoe_projector
export pauli_siteinds,
  pauli_basis_state,
  pauli_total_sz_state,
  pauli_gate,
  pauli_gate_from_hamiltonian,
  pauli_lindblad_generator,
  pauli_gate_from_lindbladian,
  DMTOptions,
  dmt_step!,
  dmt_evolve!,
  DMTGateEvolution,
  pauli_daoe_projector,
  fdaoe_projector
end

"""
Namespace for dense spin-half model-building helpers.
"""
module Models
using ..MPSToolkit: spinhalf_matrices, spinhalf_xyz_bond_hamiltonian, spinhalf_tfim_bond_hamiltonian
export spinhalf_matrices, spinhalf_xyz_bond_hamiltonian, spinhalf_tfim_bond_hamiltonian
end

"""
Namespace for Chebyshev moments, kernels, reconstruction, and energy-window cutoff tools.
"""
module Chebyshev
using ..MPSToolkit:
  ChebyshevRescaling,
  SpectralFunction,
  chebyshev_moments,
  energy_cutoff!,
  jackson_damping,
  jackson_kernel,
  reconstruct_chebyshev,
  spectral_function
export ChebyshevRescaling,
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
export pauli_matrices, pauli_basis, pauli_components
export spinhalf_matrices, spinhalf_xyz_bond_hamiltonian, spinhalf_tfim_bond_hamiltonian
export pauli_siteinds, pauli_basis_state, pauli_total_sz_state, pauli_gate, pauli_gate_from_hamiltonian, pauli_lindblad_generator, pauli_gate_from_lindbladian, DMTOptions, dmt_step!, dmt_evolve!, pauli_daoe_projector, fdaoe_projector
export ChebyshevRescaling, SpectralFunction, chebyshev_moments, energy_cutoff!, jackson_damping, jackson_kernel, reconstruct_chebyshev, spectral_function

end
