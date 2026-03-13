using ITensors
using ITensorMPS
using LinearAlgebra
using Test

function _exact_chebyshev_moments(energy::Real, order::Integer)
  moments = Vector{Float64}(undef, order)
  moments[1] = 1.0
  order == 1 && return moments

  moments[2] = energy
  for n in 3:order
    moments[n] = 2 * energy * moments[n - 1] - moments[n - 2]
  end
  return moments
end

@testset "Chebyshev moments match exact recursion" begin
  nsites = 4
  order = 8
  fields = [0.7, -0.2, 0.5, -0.4]

  sites = siteinds("S=1/2", nsites)
  psi = productMPS(sites, [isodd(j) ? "Up" : "Dn" for j in 1:nsites])

  os = OpSum()
  for j in 1:nsites
    os += fields[j], "Sz", j
  end
  h_mpo = MPO(os, sites)

  eigenvalue = sum(fields[j] * (isodd(j) ? 0.5 : -0.5) for j in 1:nsites)
  halfwidth = sum(abs, fields) / 2 + 0.1
  h_mpo = h_mpo / halfwidth

  exact = _exact_chebyshev_moments(eigenvalue / halfwidth, order)
  observed = chebyshev_moments(h_mpo, psi; order=order, maxdim=64, cutoff=1e-14)

  @test observed ≈ exact atol=1e-8
end

@testset "spectral reconstruction helpers" begin
  moments = [1.0, 0.2, -0.1, 0.05]
  kernel = jackson_kernel(length(moments))
  spectrum = spectral_function(moments; center=0.3, halfwidth=1.7)
  x = 0.25
  ω = 0.3 + 1.7 * x

  @test jackson_damping(0, length(moments)) ≈ kernel[1]
  @test reconstruct_chebyshev(x, moments; kernel=kernel) ≈ spectrum(ω) * 1.7 atol=1e-12
  @test spectrum(0.3 + 1.7 * 1.2) == 0.0
  @test MPSToolkit.Chebyshev.chebyshev_moments === chebyshev_moments
end

@testset "standalone Chebyshev energy cutoff" begin
  sites = siteinds("S=1/2", 2)
  os = OpSum()
  os += 0.5, "Sz", 1
  os += 0.5, "Sz", 2
  h_mpo = MPO(os, sites)

  in_window = productMPS(sites, ["Up", "Dn"])
  in_window_before = copy(in_window)
  energy_cutoff!(in_window, h_mpo; sweeps=2, krylovdim=4, window=0.3, tol=1e-12)
  @test abs(inner(in_window_before, in_window)) ≈ 1.0 atol=1e-8

  out_window = productMPS(sites, ["Up", "Up"])
  mixed = normalize(add(out_window, in_window; maxdim=4, cutoff=1e-14))
  energy_cutoff!(mixed, h_mpo; sweeps=4, krylovdim=4, window=0.3, tol=1e-12)
  normalize!(mixed)

  overlap_out = abs(inner(out_window, mixed))
  overlap_in = abs(inner(in_window, mixed))
  @test overlap_out < 1e-6
  @test overlap_in > 0.99

  fully_out_window = productMPS(sites, ["Up", "Up"])
  energy_cutoff!(fully_out_window, h_mpo; sweeps=4, krylovdim=4, window=0.3, tol=1e-12)
  @test isfinite(abs(inner(fully_out_window, fully_out_window)))
  @test abs(inner(fully_out_window, fully_out_window)) < 1e-12
  @test MPSToolkit.Chebyshev.energy_cutoff! === energy_cutoff!
end

@testset "Chebyshev moments energy cutoff integration" begin
  sites = siteinds("S=1/2", 2)
  os = OpSum()
  os += 0.5, "Sz", 1
  os += 0.5, "Sz", 2
  h_mpo = MPO(os, sites)

  in_window = productMPS(sites, ["Up", "Dn"])
  out_window = productMPS(sites, ["Up", "Up"])
  mixed = normalize(add(out_window, in_window; maxdim=4, cutoff=1e-14))

  default_moments = chebyshev_moments(h_mpo, mixed; order=4, maxdim=16, cutoff=1e-14)
  @test default_moments ≈ [1.0, 0.25, -0.75, -0.5] atol=1e-8

  filtered_moments = chebyshev_moments(
    h_mpo,
    mixed;
    order=4,
    maxdim=16,
    cutoff=1e-14,
    energy_cutoff=true,
    energy_cutoff_sweeps=4,
    krylovdim=4,
    window=0.3,
    energy_cutoff_tol=1e-12,
  )
  @test filtered_moments ≈ [1.0, 0.0, -0.5, 0.0] atol=1e-8
end
