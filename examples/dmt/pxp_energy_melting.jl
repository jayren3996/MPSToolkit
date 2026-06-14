# Energy transport from a melting infinitesimal energy domain wall in the PXP chain
#
#   H = Omega [ X_1 P_2 + sum_{j=2}^{N-1} P_{j-1} X_j P_{j+1} + P_{N-1} X_N ],   P = |0><0|,
#
# is the Rydberg-blockade PXP model: a site flips only when both neighbors are in the ground
# state, so the dynamics stays in the constrained sector (no two adjacent excitations) and never
# couples to blockade-violating configurations. Energy is the ONLY local conserved density, and it
# is generically DIFFUSIVE (z = 2) in the long-time limit; the KPZ value z = 3/2 is a long-lived
# INTERMEDIATE-time transient. Because bare PXP sits near a proximate integrable point, its energy
# transport stays NEAR-BALLISTIC (z ~ 1) out to t ~ 100, reaching the z = 3/2 crossover only at
# t ~ 100-300 (or sooner under a chemical-potential deformation). At the t <~ 50 reached here, an
# effective z just above 1 is the EXPECTED (literature-consistent) value, not a converged exponent:
#
#   Ljubotina, Desaules, Serbyn, Papic, "Superdiffusive energy transport in kinetically
#   constrained models", PRX 13, 011033 (2023);
#   McRoberts & Moessner, "Parametrically long lifetime of superdiffusion in non-integrable spin
#   chains", PRL 133, 256301 (2024)  (KPZ is a transient; diffusion is restored asymptotically).
#
# PROTOCOL -- infinitesimal energy-domain-wall melting (linear response). The energy analog of the
# XXZ magnetization example (examples/dmt/domain_wall_melting.jl). We prepare a weak energy domain
# wall
#       rho(0) ∝ exp(-beta H_L) ⊗ exp(+beta H_R)       (left half slightly cold, right slightly hot)
# at SMALL beta, so it is the linear-response / infinitesimal limit; terms straddling the wall get
# weight 0 so the two halves factorize exactly. The state is built by imaginary-time TEBD in
# operator space (`pauli_gibbs_state`) ON TOP OF the vectorized constraint projector P_G -- the
# infinite-temperature state of the constrained sector -- because the unconstrained identity
# carries unphysical weight on blockade-violating configurations. We evolve rho(t) = U rho U^dagger
# under DMT and read the energy profile
#       e_j(t) = tr(rho(t) h_j) / tr(rho),
# and the ENERGY transferred across the wall,
#       dE(t) = sum_{j > wall} [ e_j(t) - e_j(0) ]   ~   t^(1/z),
# whose running log-log slope inverts to the running exponent  z(t) = 1 / (d ln|dE| / d ln t).
#
# THE CONSTRAINT (the PXP-specific part). Exact PXP gates commute with P_G, but DMT TRUNCATION
# leaks Hilbert-Schmidt weight out of the constrained sector; `constrained_dmt_evolve!` re-applies
# the operator-space projector  rho -> P_G rho P_G  at every sweep. Two error bars are printed:
#   * `sector_leakage` = 1 - <rho, P_G rho P_G>/<rho, rho> -- the residual constraint violation;
#   * `energy drift`   = |sum_j e_j(t) - sum_j e_j(0)|     -- total energy is conserved, so its
#                                                            drift is the DMT truncation error bar.
#
# As in the XXZ example, z(t) is an EFFECTIVE, resource-bound exponent: PXP energy starts ballistic
# and bends to KPZ only slowly, so read z(t) in the intermediate plateau and raise
# NSITES / MAXDIM / T_MAX to sharpen.

using ITensors
using ITensorMPS
using MPSToolkit
using Printf

# ---- parameters (constrained single-case run; raise NSITES / MAXDIM / T_MAX to sharpen) ----------
const NSITES           = 60
const WALL             = NSITES ÷ 2           # energy domain wall between sites WALL and WALL+1
const BETA             = 0.1                   # weak imbalance -> linear-response (infinitesimal) wall
const OMEGA            = 1.0
const DT               = 0.1                   # per-sweep gate time; one sweep advances t by 2*DT
const T_MAX            = 16.0
const NCALL            = round(Int, T_MAX / (2 * DT))   # constrained forward+reverse sweeps
const MAXDIM           = 48
const GATE_MAXDIM      = 4 * MAXDIM
const CONNECTOR_BUFFER = 8
const CUTOFF           = 1e-12
const PROJECT_EVERY    = 1                     # constraint checkpoint after every sweep
const GIBBS_NSTEPS     = 6                     # imaginary-time Trotter steps for rho(0)
const SLOPE_HALFWIDTH  = 3                     # +/- points in the running log-log slope window
const FIT_WINDOW       = (6.0, 16.0)           # late-window fit for the (crossover) effective 1/z
const EDGE_TOL         = 0.05                  # drop fit times once the front reaches the edges
const Z_TARGET         = 1.5                   # KPZ superdiffusion
const CSV_PATH         = joinpath(@__DIR__, "pxp_energy_melting.csv")
# -------------------------------------------------------------------------------------------------

# Running log-log slope d ln(y)/d ln(t) via a sliding-window least-squares fit. For y ~ t^p the
# slope -> p; dE ~ t^(1/z), so the caller inverts z = 1 / slope. (Verbatim from domain_wall_melting.jl.)
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

# Late-window least-squares slope of |dE|(t); returns (1/z, npoints). Skips early transients,
# boundary-contaminated times, and (defensively) zeros.
function fit_inverse_z(times, transferred, edge_fraction, (tmin, tmax))
    idx = findall(
        i -> tmin <= times[i] <= tmax && abs(transferred[i]) > 0 && edge_fraction[i] < EDGE_TOL,
        eachindex(times),
    )
    length(idx) >= 2 || return (NaN, length(idx))
    lt = log.(times[idx])
    lm = log.(abs.(transferred[idx]))
    n = length(lt)
    sx = sum(lt)
    slope = (n * sum(lt .* lm) - sx * sum(lm)) / (n * sum(lt .^ 2) - sx^2)
    return slope, n
end

# Constrained PXP setup: Pauli sites, local energy terms, the weak energy domain wall rho(0)
# prepared in the constrained sector, the sector projector, and the DMT gate-evolution object.
function _build_state_and_evolution()
    sites = pauli_siteinds(NSITES)
    terms = [(first(pxp_term_support(NSITES, j)), pxp_term_hamiltonian(NSITES, j; omega=OMEGA)) for j in 1:NSITES]
    # Domain-wall weights: +beta for terms entirely in the left half, -beta entirely in the right,
    # 0 for the (two) wall-straddling terms -> rho(0) = e^{-beta H_L} ⊗ e^{+beta H_R}.
    function wall_weight(j)
        support = pxp_term_support(NSITES, j)
        last(support) <= WALL && return BETA
        first(support) > WALL && return -BETA
        return 0.0
    end
    weights = [wall_weight(j) for j in 1:NSITES]
    rho = pauli_gibbs_state(sites, terms, weights; nsteps=GIBBS_NSTEPS, maxdim=2 * MAXDIM,
        cutoff=CUTOFF, initial_state=pauli_pxp_constraint_state(sites))
    projector = pauli_pxp_constraint_projector(sites)
    gates = [pauli_gate_from_hamiltonian(h, DT) for (_, h) in terms]
    schedule = [start for (start, _) in terms]
    evo = DMTGateEvolution(gates, DT; schedule=schedule, reverse_schedule=reverse(schedule),
        maxdim=MAXDIM, cutoff=CUTOFF, gate_maxdim=GATE_MAXDIM, connector_buffer=CONNECTOR_BUFFER)
    return rho, evo, projector, terms
end

# Normalized energy-density profile e_j(t) = tr(rho h_j)/tr(rho), in one O(N) sweep.
energy_profile(state, terms) = real.(pauli_expectation_profile(state, terms))

# Hilbert-Schmidt weight outside the constrained sector: 1 - <rho, P_G rho P_G>/<rho, rho>.
function sector_leakage(state, projector)
    projected = apply(projector, state; maxdim=GATE_MAXDIM, cutoff=CUTOFF)
    return 1 - real(inner(state, projected)) / real(inner(state, state))
end

# Build, constrained-DMT-evolve, and record (t, dE, energy drift, leakage, edge) at every sweep.
function melt()
    rho, evo, projector, terms = _build_state_and_evolution()
    profile0 = energy_profile(rho, terms)
    e_right0 = sum(profile0[(WALL + 1):NSITES])
    total0 = sum(profile0)
    times = [0.0]; transferred = [0.0]; drift = [0.0]
    leak = [sector_leakage(rho, projector)]; edge = [0.0]
    for k in 1:NCALL
        constrained_dmt_evolve!(rho, evo, projector; project_every=PROJECT_EVERY)
        profile = energy_profile(rho, terms)
        delta = profile - profile0
        push!(times, 2 * DT * k)
        push!(transferred, sum(profile[(WALL + 1):NSITES]) - e_right0)        # dE(t) ~ t^(1/z)
        push!(drift, abs(sum(profile) - total0))                              # energy conservation bar
        push!(leak, sector_leakage(rho, projector))                          # constraint bar
        push!(edge, (abs(delta[1]) + abs(delta[2]) + abs(delta[end - 1]) + abs(delta[end])) / maximum(abs, delta))
    end
    return times, transferred, drift, leak, edge
end

# Warmed-up timing probe -- the cold probe under-estimates because bond dims grow as the wall melts.
function timing_probe(; warmup=12, calls=3)
    rho, evo, projector, _ = _build_state_and_evolution()
    for _ in 1:warmup
        constrained_dmt_evolve!(rho, evo, projector; project_every=PROJECT_EVERY)
    end
    t0 = time()
    for _ in 1:calls
        constrained_dmt_evolve!(rho, evo, projector; project_every=PROJECT_EVERY)
    end
    return (time() - t0) / calls
end

function main(; probe=true)
    _assert_slope_oracle()
    if probe
        per = timing_probe()
        @printf("timing probe (warmed): %.2f s / sweep  ->  full run ~ %.1f min  (%d sweeps + Gibbs prep)\n",
                per, NCALL * per / 60, NCALL)
    end

    println("\nPXP energy-domain-wall melting via constrained operator-space DMT")
    println("  nsites=$NSITES  maxdim=$MAXDIM  connector_buffer=$CONNECTOR_BUFFER  beta=$BETA  dt=$DT  t_max=$T_MAX  wall=$WALL")
    t, dE, drift, leak, edge = melt()
    z = 1.0 ./ running_loglog_slope(t, abs.(dE))

    println(rpad("t", 8), rpad("dE(t)", 14), rpad("z(t)", 9), rpad("energy drift", 14),
            rpad("leakage", 12), "front@edges")
    for k in eachindex(t)
        (k == 1 || (k - 1) % 5 == 0) || continue
        ztext = isfinite(z[k]) ? @sprintf("%.3f", z[k]) : "--"
        @printf("%-8.2f%-14.4e%-9s%-14.2e%-12.2e%.2e\n", t[k], dE[k], ztext, drift[k], leak[k], edge[k])
    end

    inv_z, npts = fit_inverse_z(t, dE, edge, FIT_WINDOW)
    zlast = let zf = filter(isfinite, z); isempty(zf) ? NaN : zf[end] end
    println()
    @printf("z(target)=%.1f   z(t_max)=%.3f   late-window fit z_eff=%.2f over t in %s (%d pts)\n",
            Z_TARGET, zlast, isnan(inv_z) ? NaN : 1 / inv_z, FIT_WINDOW, npts)

    open(CSV_PATH, "w") do io
        println(io, "t,dE,z,drift,leakage,edge")
        for k in eachindex(t)
            @printf(io, "%.4f,%.6e,%.4f,%.3e,%.3e,%.3e\n", t[k], dE[k], z[k], drift[k], leak[k], edge[k])
        end
    end
    println("wrote $CSV_PATH")
    println("z(t) is an effective, resource-bound exponent: PXP energy starts ballistic (z~1) and")
    println("DRIFTS UP toward the KPZ value 3/2. Raise NSITES / MAXDIM / T_MAX to sharpen; the")
    println("leakage and energy-drift columns are the constraint and truncation error bars.")
    return (; t, dE, z, drift, leak, edge)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
