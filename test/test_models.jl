using ITensors
using ITensorMPS

@testset "spinhalf model helpers" begin
  mats = spinhalf_matrices()
  @test mats.Sx ≈ ComplexF64[0 0.5; 0.5 0]
  @test mats.Sz ≈ ComplexF64[0.5 0; 0 -0.5]

  xyz = spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=2.0, Jz=3.0)
  @test size(xyz) == (4, 4)

  tfim = spinhalf_tfim_bond_hamiltonian(6, 1; J=1.0, g=0.5)
  @test size(tfim) == (4, 4)
end
