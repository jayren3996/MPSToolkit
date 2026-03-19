"""
    spinhalf_matrices(; include_identity=true)

Return dense spin-1/2 matrices in the `(I, Sx, Sy, Sz)` convention.

# Keyword Arguments
- `include_identity`: If `true`, include the identity matrix as field `I`.

# Returns
- A named tuple of dense `2 x 2` matrices where `Sx = σx / 2`, `Sy = σy / 2`, and
  `Sz = σz / 2`.
"""
function spinhalf_matrices(; include_identity::Bool=true)
  paulis = pauli_matrices(; include_identity=include_identity)
  if include_identity
    return (I=paulis.I, Sx=paulis.X / 2, Sy=paulis.Y / 2, Sz=paulis.Z / 2)
  end
  return (Sx=paulis.X / 2, Sy=paulis.Y / 2, Sz=paulis.Z / 2)
end

"""
    spinhalf_xyz_bond_hamiltonian(; Jx=0.0, Jy=0.0, Jz=0.0)

Return the dense two-site spin-1/2 XYZ bond Hamiltonian
`Jx Sx⊗Sx + Jy Sy⊗Sy + Jz Sz⊗Sz`.

# Keyword Arguments
- `Jx`, `Jy`, `Jz`: Coupling strengths multiplying the corresponding spin-spin terms.

# Returns
- A dense `4 x 4` matrix acting on two spin-1/2 sites.
"""
function spinhalf_xyz_bond_hamiltonian(; Jx::Real=0.0, Jy::Real=0.0, Jz::Real=0.0)
  spins = spinhalf_matrices()
  return Jx * kron(spins.Sx, spins.Sx) + Jy * kron(spins.Sy, spins.Sy) + Jz * kron(spins.Sz, spins.Sz)
end

"""
    spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=1.0)

Return the dense two-site TFIM bond Hamiltonian with open-boundary field splitting:
`-J Sz⊗Sz - g (w_left Sx⊗I + w_right I⊗Sx)`.

# Arguments
- `nsites`: Number of sites in the open chain.
- `bond`: Bond index of the two-site term, satisfying `1 <= bond < nsites`.

# Keyword Arguments
- `J`: Ising coupling.
- `g`: Transverse-field strength.

# Returns
- A dense `4 x 4` matrix suitable for TEBD helper constructors.

# Notes
- The transverse field is split evenly between neighboring bonds in the bulk and assigned
  fully on the edges, which makes the sum over bond Hamiltonians reproduce the standard
  open-chain TFIM Hamiltonian.
"""
function spinhalf_tfim_bond_hamiltonian(nsites::Integer, bond::Integer; J::Real=1.0, g::Real=1.0)
  1 <= bond < nsites || throw(ArgumentError("bond index must satisfy 1 <= bond < nsites"))
  spins = spinhalf_matrices()
  left_weight = bond == 1 ? 1.0 : 0.5
  right_weight = bond == nsites - 1 ? 1.0 : 0.5
  return -J * kron(spins.Sz, spins.Sz) - g * (left_weight * kron(spins.Sx, spins.I) + right_weight * kron(spins.I, spins.Sx))
end
