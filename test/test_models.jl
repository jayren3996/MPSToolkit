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
end
