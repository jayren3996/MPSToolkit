"""
    pxp_term_support(nsites, site)

Return the physical-site range acted on by the PXP term assigned to `site` on an open chain.

# Arguments
- `nsites`: Number of sites in the open chain, at least 2.
- `site`: Center site of the term, satisfying `1 <= site <= nsites`.

# Returns
- The `UnitRange` `max(site-1, 1):min(site+1, nsites)`. Its first entry is the gate start
  position expected by the TEBD/DMT schedulers.

# Notes
- Bulk terms span three sites; the two edge terms (`site == 1` and `site == nsites`) span two.
"""
function pxp_term_support(nsites::Integer, site::Integer)
  nsites >= 2 || throw(ArgumentError("PXP requires at least two sites"))
  1 <= site <= nsites || throw(ArgumentError("PXP term site must lie in 1:nsites"))
  return max(Int(site) - 1, 1):min(Int(site) + 1, Int(nsites))
end

"""
    pxp_term_hamiltonian(nsites, site; omega=1.0)

Return the dense PXP Hamiltonian term assigned to `site` on an open chain.

The local basis is `(|0>, |1>)` with `|0>` the ground state, the projector is
`P = |0><0| = (I + Z) / 2`, and the open-chain Hamiltonian is

`H = Ω [ X_1 P_2 + Σ_{j=2}^{N-1} P_{j-1} X_j P_{j+1} + P_{N-1} X_N ]`,

one term per site. Summing the returned blocks embedded at
[`pxp_term_support`](@ref) over all `site` values reproduces `H` exactly.

# Arguments
- `nsites`: Number of sites in the open chain, at least 2.
- `site`: Center site of the term.

# Keyword Arguments
- `omega`: Overall Rabi-frequency prefactor `Ω`.

# Returns
- A dense `8 x 8` matrix in the bulk and a dense `4 x 4` matrix for the two edge terms, with
  kron factor order equal to chain order.
"""
function pxp_term_hamiltonian(nsites::Integer, site::Integer; omega::Real=1.0)
  support = pxp_term_support(nsites, site)
  paulis = pauli_matrices()
  ground = (paulis.I + paulis.Z) / 2
  factors = [s == site ? paulis.X : ground for s in support]
  return omega * foldl(kron, factors)
end

"""
    _pxp_constraint_transition(state, local_state)

Transition rule of the two-state PXP constraint automaton ("was the previous site excited").
Physical basis labels are `1 = |0>` and `2 = |1>`; adjacent excitations get coefficient zero.
"""
function _pxp_constraint_transition(state::Int, local_state::Int)
  local_state == 1 && return 1, 1.0
  return state == 1 ? (2, 1.0) : (2, 0.0)
end

"""
    pxp_constraint_mpo(sites)

Build the PXP (Rydberg-blockade) constraint projector `P_G = Π_j (1 - n_j n_{j+1})` as a
diagonal MPO with bond dimension 2 on physical dimension-2 sites.

# Arguments
- `sites`: Physical site indices with local dimension 2, ordered `(|0>, |1>)`.

# Returns
- A diagonal `MPO` that annihilates any configuration with two adjacent excitations and acts
  as the identity on the constrained subspace.

# Notes
- `P_G` commutes with every PXP term, so exact PXP evolution preserves the constrained
  sector; the projector is needed to remove truncation-induced leakage. See
  [`pauli_pxp_constraint_projector`](@ref) for the operator-space version.
"""
function pxp_constraint_mpo(sites)
  length(sites) >= 1 || throw(ArgumentError("PXP constraint MPO requires at least one site"))
  all(s -> dim(s) == 2, sites) || throw(ArgumentError("PXP constraint MPO requires dimension-2 sites"))
  return _diagonal_projector_mpo(sites, 2, _pxp_constraint_transition)
end
