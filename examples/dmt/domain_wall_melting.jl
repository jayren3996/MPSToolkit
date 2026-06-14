# Charge transport from a melting infinitesimal domain wall in the spin-1/2 XXZ chain
#
#   H = sum_j [ Sx_j Sx_{j+1} + Sy_j Sy_{j+1} + Delta Sz_j Sz_{j+1} ]      (Sa = sigma^a / 2)
#
# The high-temperature spin transport of the XXZ chain has three regimes set by the anisotropy
# Delta, characterized by the dynamical exponent z (the magnetization profile spreads as t^(1/z)):
#
#   Delta < 1 (easy-plane)   ballistic            z = 1
#   Delta = 1 (Heisenberg)   superdiffusive (KPZ) z = 3/2
#   Delta > 1 (easy-axis)    diffusive            z = 2
#
# PROTOCOL -- infinitesimal domain-wall melting (linear response, infinite temperature).
# We prepare the traceless domain-wall operator  D = sum_j sign(j - kink) sigma^z_j  (the O(eps)
# part of the weakly polarized step state  rho ~ I + eps D) and evolve it under DMT with the
# Schroedinger map  D -> U D U^dagger,  U = prod_b exp(-i H_b dt). The wall melts. We read the
# magnetization profile
#       m(x,t) = tr(D(t) sigma^z_x) / anchor(t),
# normalized by the unmelted bulk so the far field reads +-1, and the CHARGE that has crossed the
# central bond,
#       T(t) = sum_{x > kink} [ 1 - m(x,t) ]   ~   t^(1/z).
# The running dynamical exponent is the logarithmic derivative
#       z(t) = 1 / ( d ln T / d ln t ),
# which starts near 1 (every Delta looks ballistic while the front is still local) and DRIFTS
# toward the regime target as the hydrodynamic tail develops.
#
# WHY THESE CHOICES.
#  * D is TRACELESS, so we evolve with normalize=false (per-sweep normalize! would rescale to unit
#    Frobenius norm, which drifts under truncation and corrupts absolute charges).
#  * anchor(t) = (1/2)[mean(p over rightmost K) - mean(p over leftmost K)] is a RATIO, hence
#    invariant to any overall rescaling (the unknown Pauli-basis prefactor, residual norm drift).
#  * sum_x m(x,t) is conserved at 0 (total S^z commutes with every XXZ bond); its drift away from
#    0 is the truncation error bar, printed as max|Sigma m|.
#  * B1 caveat (documented in docs/src/manual/dmt.md): DMT connector preservation is O(per-bond
#    truncation) at connector_buffer > 0, so we use connector_buffer = 4 and a converged maxdim,
#    read z(t) in the intermediate plateau, and treat it as an effective, resource-bound exponent.
#
# References:
#   - Ljubotina, Znidaric, Prosen, PRL 122, 210602 (2019)  (KPZ in the Heisenberg magnet).
#   - White, Zaletel, Mong, Refael, PRB 97, 035127 (2018)  (DMT).
#   - Yi-Thomas, Ware, Sau, White, arXiv:2310.06886         (DMT/DAOE hydrodynamics benchmark).

using ITensors
using ITensorMPS
using MPSToolkit
using Printf

# ---- parameters (default = cheap; raise NSITES / MAXDIM / T_MAX for a sharper "detailed" run) ----
const NSITES           = 80
const MAXDIM           = 32
const DT               = 0.1                 # per-sweep gate time; one evolve! advances t by 2*DT
const T_MAX            = 12.0
const NCALL            = round(Int, T_MAX / (2 * DT))   # forward+reverse DMT sweeps
const CONNECTOR_BUFFER = 4
const GATE_MAXDIM      = 4 * MAXDIM
const CUTOFF           = 1e-10
const K_ANCHOR         = 4                    # outermost sites (each side) defining the bulk anchor
const SLOPE_HALFWIDTH  = 3                    # +/- points in the running log-log slope window
const KINK             = NSITES ÷ 2           # domain wall on the central bond
const CASES            = ((0.5, "ballistic", 1.0), (1.0, "superdiffusive (KPZ)", 1.5),
                          (2.0, "diffusive", 2.0))
const CSV_PATH         = joinpath(@__DIR__, "domain_wall_melting.csv")
# Physical Pauli-Z (2x2) for the single-site magnetization measurement tr(D sigma^z_x).
const PAULI_Z          = pauli_matrices().Z
# -------------------------------------------------------------------------------------------------

# Running log-log slope d ln(y) / d ln(t) via a sliding-window least-squares fit. For y ~ t^p the
# slope -> p; the transferred charge obeys T ~ t^(1/z), so the caller inverts z = 1 / slope. Pure
# function -- its power-law oracle (`_assert_slope_oracle`) makes it unit-testable with no physics.
function running_loglog_slope(t::AbstractVector, y::AbstractVector; halfwidth::Int=SLOPE_HALFWIDTH)
    n = length(t)
    slope = fill(NaN, n)
    for k in 1:n
        lo, hi = k - halfwidth, k + halfwidth
        (lo >= 1 && hi <= n) || continue
        idx = [i for i in lo:hi if t[i] > 0 && y[i] > 0]
        length(idx) >= 2 || continue
        lt = log.(t[idx])
        ly = log.(y[idx])
        m = length(lt)
        sx = sum(lt)
        sy = sum(ly)
        denom = m * sum(abs2, lt) - sx^2
        denom == 0 && continue
        slope[k] = (m * sum(lt .* ly) - sx * sy) / denom
    end
    return slope
end

# Cheap self-check: an exact power law must recover its exponent to ~machine precision.
function _assert_slope_oracle()
    t = collect(1.0:0.5:12.0)
    s = running_loglog_slope(t, 3.0 .* t .^ (2 / 3); halfwidth=3)
    finite = filter(isfinite, s)
    @assert !isempty(finite) && maximum(abs.(finite .- 2 / 3)) < 1e-9 "running_loglog_slope power-law oracle failed"
    return nothing
end

# Bulk-anchored magnetization profile and the charge transferred across the central bond.
function profile_and_charge(state, kink, K)
    n = length(state)
    p = real.(pauli_expectation_profile(state, [(x, PAULI_Z) for x in 1:n]; normalize=false))
    anchor = (sum(@view p[(n - K + 1):n]) - sum(@view p[1:K])) / (2 * K)  # ~ +scale; right->+1, left->-1
    m = p ./ anchor
    transferred = sum(1.0 - m[x] for x in (kink + 1):n)                   # T(t) = sum_{x>kink}[1 - m]
    sigma_m = sum(m)                                                      # conservation oracle: ~ 0
    band = (n - 2K + 1):(n - K)                                           # K sites just inside the anchor
    front = maximum(abs(1.0 - m[x]) for x in band)                       # front-reached-the-edge guard
    return transferred, sigma_m, front
end

function _build_evolution(Delta)
    sites = pauli_siteinds(NSITES)
    D = pauli_domain_wall_state(sites; kink=KINK)
    hbond = spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=1.0, Jz=Delta)
    gate = pauli_gate_from_hamiltonian(hbond, DT)
    schedule = collect(1:(NSITES - 1))
    evo = DMTGateEvolution(gate, DT; schedule=schedule, reverse_schedule=reverse(schedule),
        maxdim=MAXDIM, cutoff=CUTOFF, gate_maxdim=GATE_MAXDIM,
        connector_buffer=CONNECTOR_BUFFER, normalize=false)  # traceless operator -> no rescaling
    return D, evo
end

# One regime: build D, DMT-evolve, record (t, T, |Sigma m|, front) at every sweep.
function melt(Delta)
    D, evo = _build_evolution(Delta)
    T0, s0, f0 = profile_and_charge(D, KINK, K_ANCHOR)
    times = [0.0]; charge = [T0]; drift = [abs(s0)]; front = [f0]
    for k in 1:NCALL
        evolve!(D, evo)
        T, sm, fr = profile_and_charge(D, KINK, K_ANCHOR)
        push!(times, 2 * DT * k); push!(charge, T); push!(drift, abs(sm)); push!(front, fr)
    end
    return times, charge, drift, front
end

# Extrapolate the full 3-regime cost from a few warmed-up sweeps before committing to the run.
function timing_probe(; calls=3)
    D, evo = _build_evolution(1.0)
    evolve!(D, evo)  # warm up / compile
    t0 = time()
    for _ in 1:calls
        evolve!(D, evo)
    end
    return (time() - t0) / calls
end

function main(; probe=true)
    _assert_slope_oracle()
    if probe
        per = timing_probe()
        @printf("timing probe: %.2f s / evolve!  ->  full run ~ %.1f min  (3 regimes x %d sweeps)\n",
                per, 3 * NCALL * per / 60, NCALL)
    end

    println("\nDomain-wall melting via operator-space DMT")
    println("  nsites=$NSITES  maxdim=$MAXDIM  connector_buffer=$CONNECTOR_BUFFER  dt=$DT  t_max=$T_MAX  kink=$KINK")
    println(rpad("Delta", 7), rpad("regime", 22), rpad("z(target)", 11), rpad("z(t_max)", 11),
            rpad("T(t_max)", 11), "max|Sigma m|")

    results = Dict{Float64,NamedTuple}()
    for (Delta, regime, zpred) in CASES
        t, T, drift, front = melt(Delta)
        @assert abs(T[1]) < 1e-9 "T(0) should vanish for an unmelted wall"   # structural oracle
        z = 1.0 ./ running_loglog_slope(t, T)
        results[Delta] = (; t, T, z, drift, front)
        zlast = let zf = filter(isfinite, z); isempty(zf) ? NaN : zf[end] end
        @printf("%-7.2f%-22s%-11.1f%-11.3f%-11.3f%.2e\n",
                Delta, regime, zpred, zlast, T[end], maximum(drift))
    end

    t = results[0.5].t
    open(CSV_PATH, "w") do io
        println(io, "t,T_d05,z_d05,drift_d05,T_d10,z_d10,drift_d10,T_d20,z_d20,drift_d20")
        for k in eachindex(t)
            @printf(io, "%.4f,%.6e,%.4f,%.3e,%.6e,%.4f,%.3e,%.6e,%.4f,%.3e\n", t[k],
                results[0.5].T[k], results[0.5].z[k], results[0.5].drift[k],
                results[1.0].T[k], results[1.0].z[k], results[1.0].drift[k],
                results[2.0].T[k], results[2.0].z[k], results[2.0].drift[k])
        end
    end
    println("\nwrote $CSV_PATH")
    println("z(t) is an effective, resource-bound exponent: it starts near 1 (ballistic front) and")
    println("DRIFTS toward the targets within reach of t_max. Raise NSITES / MAXDIM / T_MAX to sharpen.")
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
