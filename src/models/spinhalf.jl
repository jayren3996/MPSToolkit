"""
    spinhalf_matrices(; include_identity=true)

Return dense spin-1/2 matrices in the `(I, Sx, Sy, Sz)` convention.
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
"""
function spinhalf_xyz_bond_hamiltonian(; Jx::Real=0.0, Jy::Real=0.0, Jz::Real=0.0)
  spins = spinhalf_matrices()
  return Jx * kron(spins.Sx, spins.Sx) + Jy * kron(spins.Sy, spins.Sy) + Jz * kron(spins.Sz, spins.Sz)
end

"""
    spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=1.0)

Return the dense two-site TFIM bond Hamiltonian with open-boundary field splitting:
`-J Sz⊗Sz - g (w_left Sx⊗I + w_right I⊗Sx)`.
"""
function spinhalf_tfim_bond_hamiltonian(nsites::Integer, bond::Integer; J::Real=1.0, g::Real=1.0)
  1 <= bond < nsites || throw(ArgumentError("bond index must satisfy 1 <= bond < nsites"))
  spins = spinhalf_matrices()
  left_weight = bond == 1 ? 1.0 : 0.5
  right_weight = bond == nsites - 1 ? 1.0 : 0.5
  return -J * kron(spins.Sz, spins.Sz) - g * (left_weight * kron(spins.Sx, spins.I) + right_weight * kron(spins.I, spins.Sx))
end
