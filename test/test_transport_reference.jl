using Test

# Reproducibility guard (2026-06-13 transport-review follow-up): the PXP short-term
# superdiffusion claim must be reproducible FROM THE REPO, not from a volatile scratch
# directory. This re-fits the committed chi=48 production output (N=64, normalize=false;
# examples/operator_space/data/pxp_energy_correlator_chi48_N64.csv, re-fit by
# examples/operator_space/pxp_energy_correlator_reference.jl) so a regression in the committed
# data or the fit is caught in CI. No physics is run here — only the committed reference is re-fit.

@testset "committed chi=48 correlator reference pins the superdiffusive plateau" begin
  data = joinpath(dirname(@__DIR__), "examples", "operator_space", "data",
                  "pxp_energy_correlator_chi48_N64.csv")
  @test isfile(data)
  lines = readlines(data)
  t = Float64[]
  M2 = Float64[]
  total = Float64[]
  for line in lines[2:end]
    isempty(strip(line)) && continue
    f = split(line, ',')
    push!(t, parse(Float64, f[2]))
    push!(M2, parse(Float64, f[3]))
    push!(total, parse(Float64, f[4]))
  end
  @test length(t) >= 60
  @test issorted(M2)                              # M2 spreads monotonically
  # LSQ log-log slope over the converged short-term window t in [8,14]; the exponent is 2/z.
  idx = findall(i -> 8.0 <= t[i] <= 14.0, eachindex(t))
  lt = log.(t[idx])
  lm = log.(M2[idx])
  n = length(lt)
  sx = sum(lt)
  slope = (n * sum(lt .* lm) - sx * sum(lm)) / (n * sum(lt .^ 2) - sx^2)
  @test 1.40 <= slope <= 1.55                     # superdiffusive (>> diffusive 1.0), measured ~1.48
  # The normalize=false conserved-total drift IS the truncation error bar: small early, grows late.
  drift(tt) = (k = findfirst(≈(tt), t); abs(total[k] - total[1]) / abs(total[1]))
  @test drift(8.0) < 0.05
  @test drift(14.0) > 0.08
end
