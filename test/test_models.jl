using ITensors
using ITensorMPS

@testset "spinhalf model helpers" begin
  mats = spinhalf_matrices()
  @test mats.Sx ≈ ComplexF64[0 0.5; 0.5 0]
  @test mats.Sz ≈ ComplexF64[0.5 0; 0 -0.5]

  xyz = spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=2.0, Jz=3.0)
  @test size(xyz) == (4, 4)
  sx = ComplexF64[0 0.5; 0.5 0]
  sy = ComplexF64[0 -0.5im; 0.5im 0]
  sz = ComplexF64[0.5 0; 0 -0.5]
  @test xyz ≈ kron(sx, sx) + 2 * kron(sy, sy) + 3 * kron(sz, sz)

  tfim = spinhalf_tfim_bond_hamiltonian(6, 1; J=1.0, g=0.5)
  @test size(tfim) == (4, 4)
  id2 = Matrix{ComplexF64}(I, 2, 2)
  @test tfim ≈ -kron(sz, sz) - 0.5 * (kron(sx, id2) + 0.5 * kron(id2, sx))

  bulk = spinhalf_tfim_bond_hamiltonian(6, 3; J=1.0, g=0.5)
  @test bulk ≈ -kron(sz, sz) - 0.5 * (0.5 * kron(sx, id2) + 0.5 * kron(id2, sx))

  right_edge = spinhalf_tfim_bond_hamiltonian(6, 5; J=1.0, g=0.5)
  @test right_edge ≈ -kron(sz, sz) - 0.5 * (0.5 * kron(sx, id2) + kron(id2, sx))
  @test_throws ArgumentError spinhalf_tfim_bond_hamiltonian(6, 0)
  @test_throws ArgumentError spinhalf_tfim_bond_hamiltonian(6, 6)

  paulis = pauli_matrices()
  mfi_left = spinhalf_mixed_field_ising_bond_hamiltonian(4, 1; J=0.7, gx=1.1, gz=-0.3)
  @test mfi_left ≈ 0.7 * kron(paulis.Z, paulis.Z) +
    1.1 * (kron(paulis.X, paulis.I) + 0.5 * kron(paulis.I, paulis.X)) -
    0.3 * (kron(paulis.Z, paulis.I) + 0.5 * kron(paulis.I, paulis.Z))
  mfi_bulk = spinhalf_mixed_field_ising_bond_hamiltonian(4, 2; J=0.7, gx=1.1, gz=-0.3)
  @test mfi_bulk ≈ 0.7 * kron(paulis.Z, paulis.Z) +
    0.55 * (kron(paulis.X, paulis.I) + kron(paulis.I, paulis.X)) -
    0.15 * (kron(paulis.Z, paulis.I) + kron(paulis.I, paulis.Z))
  @test_throws ArgumentError spinhalf_mixed_field_ising_bond_hamiltonian(4, 0)
end
