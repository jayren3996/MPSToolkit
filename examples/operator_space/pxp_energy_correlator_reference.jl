# Committed chi=48 reference for PXP energy-correlator superdiffusion.
#
# This is a companion to pxp_energy_correlator.jl. The shipped demo there defaults to
# maxdim=32, which UNDERSHOOTS the KPZ exponent: at chi=32 truncation suppresses M2 growth, so
# the local slope lands just below the KPZ value 4/3 by cancellation, not by reaching the
# asymptote. The converged short-term signal only appears at chi >= 48, and the production run
# that establishes it (N=64, normalize=false, t up to 14) is expensive. To keep the
# superdiffusion claim reproducible FROM THE REPO rather than from a volatile scratch
# directory, the output of that run is committed here and re-fit by this script.
#
# Provenance — examples/operator_space/data/pxp_energy_correlator_chi48_N64.csv:
#   N=64, maxdim=48 (gate-level bond caps at 2*maxdim=96), normalize=false, project_every=1,
#   dt=0.1, t in 0.2 .. 14.0. Generated 2026-06-12 by the correlator protocol.
#   Columns: k, t, M2, total (conserved tr(P_G H O(t))), symerr, chi (achieved bond), secs.
#
# Reference results (see docs/superpowers/audits/2026-06-13-transport-review.md): the local M2
# slope plateaus at 1.46-1.50 over t=8-14 (LSQ 2/z = 1.480, z_eff ~ 1.35), robustly above the
# diffusive value 2/z = 1.0 and descending toward the KPZ value 4/3 = 1.333 FROM ABOVE. The
# conserved-total drift (the truncation error bar with normalize=false) grows 1.8% (t=8) ->
# 12.2% (t=14), so the plateau is trustworthy only to t ~ 10; the z=3/2 asymptote needs
# t >= 20-30 at chi > 48 and is NOT established here. test/test_pxp.jl pins this fit in CI.

using Printf

const DATA = joinpath(@__DIR__, "data", "pxp_energy_correlator_chi48_N64.csv")

function load_reference(path)
    lines = readlines(path)
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
    return t, M2, total
end

# LSQ slope of log(M2) vs log(t) over [tmin, tmax]; the exponent is 2/z (4/3 for KPZ z=3/2).
function fit_slope(t, M2, tmin, tmax)
    idx = findall(i -> tmin <= t[i] <= tmax && M2[i] > 0, eachindex(t))
    lt = log.(t[idx])
    lm = log.(M2[idx])
    n = length(lt)
    sx = sum(lt)
    slope = (n * sum(lt .* lm) - sx * sum(lm)) / (n * sum(lt .^ 2) - sx^2)
    return slope, n
end

# Local 2/z over a ~1.0 time-unit gap, matching the per-row slope pxp_energy_correlator.jl prints.
local_slope(t, M2, k, gap) = log(M2[k] / M2[k-gap]) / log(t[k] / t[k-gap])

t, M2, total = load_reference(DATA)
gap = 5  # report step is 2*dt = 0.2 t-unit, so gap*0.2 = 1.0 t-unit
drift(k) = abs(total[k] - total[1]) / abs(total[1])

println("PXP energy-correlator chi=48 reference (N=64, normalize=false)  [committed data]")
println(rpad("t", 8), rpad("M2", 12), rpad("local 2/z", 12), "drift (error bar)")
for k in eachindex(t)
    (abs(t[k] - round(t[k])) < 1e-9 && 6 <= t[k] <= 14) || continue
    ls = k > gap ? @sprintf("%.3f", local_slope(t, M2, k, gap)) : "--"
    @printf("%-8.1f%-12.4f%-12s%.2f%%\n", t[k], M2[k], ls, 100 * drift(k))
end

slope, n = fit_slope(t, M2, 8.0, 14.0)
println()
@printf("LSQ fit over t in [8,14] (%d pts):  2/z = %.3f  =>  z_eff = %.2f\n", n, slope, 2 / slope)
println("""
Diffusive 2/z = 1.0; KPZ (z = 3/2) 2/z = 4/3 = 1.333. The measured 2/z ~ 1.48 is robustly
superdiffusive and descends toward 4/3 FROM ABOVE; it has NOT reached the asymptote (which
needs t >= 20-30 at chi > 48). The drift column is the normalize=false truncation error bar:
trustworthy to ~6% only through t ~ 10. See docs/superpowers/audits/2026-06-13-transport-review.md.""")
