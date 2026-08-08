"""
    _pxp_physical_sites(nsites)

Create the internal physical dimension-2 site indices used to build PXP constraint objects in
operator space.
"""
function _pxp_physical_sites(nsites::Integer)
  return [Index(2, "PXPSite,n=$(n)") for n in 1:Int(nsites)]
end

"""
    pauli_pxp_constraint_state(sites)

Return the vectorized PXP constraint projector `P_G = Π_j (1 - n_j n_{j+1})` as a Pauli-basis
operator-space `MPS` with bond dimension 2.

Up to normalization this is the infinite-temperature density matrix of the **constrained**
sector, `ρ_∞ ∝ P_G`, and is the natural `initial_state` for thermal preparation with
[`pauli_gibbs_state`](@ref) in PXP transport calculations: the unconstrained identity state
carries weight on blockade-violating configurations that the PXP Hamiltonian never couples
to physical dynamics.

# Arguments
- `sites`: Pauli-space site indices (dimension 4), typically from [`pauli_siteinds`](@ref).

# Returns
- An `MPS` over `sites` whose Pauli amplitudes are `tr(P_α† P_G)`.

# Notes
- Defined for spin-1/2 (local dimension 2) operator space only; throws `ArgumentError` for
  any other local dimension, since the PXP blockade constraint is a two-level construction.
"""
function pauli_pxp_constraint_state(sites)
  all(local_dimension(site) == 2 for site in sites) ||
    throw(ArgumentError("this helper is defined for spin-1/2 (local dimension 2) operator space only"))
  return pauli_state_from_mpo(pxp_constraint_mpo(_pxp_physical_sites(length(sites))), sites)
end

"""
    pauli_pxp_constraint_projector(sites)

Return the operator-space MPO implementing the two-sided PXP constraint projection
`ρ ↦ P_G ρ P_G` in the normalized Pauli basis, with bond dimension 4.

Every PXP term commutes with `P_G`, so exact evolution never leaves the constrained sector;
truncation (DMT or plain SVD) does leak weight out of it. Periodically applying this MPO
("checkpoints") removes the leaked weight — see [`constrained_dmt_evolve!`](@ref).

# Arguments
- `sites`: Pauli-space site indices (dimension 4).

# Returns
- An `MPO` over `sites`, idempotent up to truncation, that fixes any vectorized operator
  supported in the constrained sector and annihilates operators supported on
  blockade-violating configurations.

# Notes
- Defined for spin-1/2 (local dimension 2) operator space only; throws `ArgumentError` for
  any other local dimension, since the PXP blockade constraint is a two-level construction.
"""
function pauli_pxp_constraint_projector(sites)
  all(local_dimension(site) == 2 for site in sites) ||
    throw(ArgumentError("this helper is defined for spin-1/2 (local dimension 2) operator space only"))
  return pauli_superoperator_mpo(pxp_constraint_mpo(_pxp_physical_sites(length(sites))), sites)
end
