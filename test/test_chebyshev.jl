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
