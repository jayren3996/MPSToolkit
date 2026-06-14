# Energy diffusion in the mixed-field (tilted-field) Ising chain -- the canonical DMT benchmark
#
#   H = sum_i [ J Z_i Z_{i+1} + g_x X_i + g_z Z_i ]      (Pauli operators, eigenvalues +-1)
#
# with J = 1, g_x = 1.4, g_z = 0.9045 is the strongly chaotic, nonintegrable point on which the
# DMT method was developed and benchmarked. Energy is the ONLY conserved density (there is no U(1)),
# and it spreads DIFFUSIVELY: the dynamical exponent is z = 2 and the infinite-temperature energy
# diffusion constant is
#
#       D ~= 1.446
#
# This number is not from the original DMT paper (which reports method accuracy, not a transport
# coefficient) but from the multi-method comparison that pins it with DMT, DAOE/OST and TEBD all
# agreeing to ~1%:
#
#   White, Zaletel, Mong, Refael, "Quantum dynamics of thermalizing systems", PRB 97, 035127 (2018)
#       -- DMT, validated on THIS model (accuracy ~1e-3 vs semi-exact MPS; no D / z reported).
#   Yi-Thomas, Ware, Sau, White, arXiv:2310.06886 (PRB 2024)
#       -- the cross-method number D ~= 1.446 (z = 2), DMT converged to 0.23% at chi = 256.
#   Rakovszky, von Keyserlingk, Pollmann, PRB 105, 075131 (2022), arXiv:2004.05177
#       -- DAOE gives D ~= 1.40 at the same field point.
#
# This example verifies BOTH halves of the benchmark with the two protocols that match the two
# observables -- and they are deliberately routed through the two DIFFERENT operator-space drivers,
# because the trace structure of the operator decides which driver is valid:
#
#   (A) MELT -> z   [DMT].   rho(0) ~ exp(-beta H_L) (x) exp(+beta H_R), a weak energy domain wall,
#       is a near-infinite-temperature DENSITY operator (trace-ful). DMT protects the trace/identity
#       and the local Pauli data, so it is the correct driver. We read the energy transferred across
#       the wall, dE(t) ~ t^(1/z), and invert the running log-log slope to z(t) -> 2. This validates
#       that the DMT machinery recovers the diffusive exponent on a chaotic model.
#
#   (B) SPREAD -> D [TEBD + DAOE].  The central energy density h_c(t), evolved in the Heisenberg
#       picture, is a TRACELESS operator. DMT must NOT be used on it (it protects a trace a traceless
#       operator does not have); we use ordinary operator-space TEBD with no renormalization. We
#       measure the energy-density autocorrelation
#               C_b(t) = tr[ h_c(t) h_b ]                       (one O(N) sweep, normalize=false)
#       whose spatial variance grows as  Var(t) -> 2 D t,  giving D from a late-window linear fit.
#       (Var is the second-moment / mean-square-displacement form of the energy-current Green-Kubo
#       integral; sum_b C_b(t) = tr[h_c(t) H] is conserved -> a built-in truncation error bar.)
#       This chaotic model's operator entanglement defeats *pure* TEBD (it leaks O(5%) of the
#       conserved energy by t ~ 10), so we add DAOE dissipation -- damping Pauli strings longer than
#       l* by gamma -- which is precisely how Rakovszky et al. obtained D ~= 1.40. DAOE controls the
#       truncation (the M0 error bar drops to <1%) at the cost of a transport bias that vanishes as
#       gamma -> 0 / l* -> infinity, so the printed D is a resource-bound LOWER estimate that climbs
#       toward 1.446 as the dissipation is relaxed.
#
# This split is exactly the density-vs-traceless guardrail enforced by `dmt_evolve!` /
# `constrained_dmt_evolve!`: melt the density with DMT, spread the traceless density with TEBD(+DAOE).
#
# As with every finite-resource transport run, z(t) and D(t) are EFFECTIVE, resource-bound: the
# diffusive plateau emerges only after the fast (non-conserved) modes decay (t ~ a few / Omega) and
# before the front reaches the chain edges. Raise NSITES / MAXDIM / T_MAX (and, for D, lower
# MSD_DAOE_GAMMA / raise MSD_DAOE_LSTAR, extrapolating to the dissipationless limit) to sharpen; the
# energy-drift / M0-drift columns are the truncation error bars.

using ITensors
using ITensorMPS
using MPSToolkit
using Printf

# ---- model (Pauli convention; the field point with the published D) -----------------------------
const J  = 1.0
const GX = 1.4
const GZ = 0.9045
const D_TARGET = 1.446        # arXiv:2310.06886 (cross-method, z = 2)
const Z_TARGET = 2.0

# ---- shared discretization ----------------------------------------------------------------------
const NSITES           = 60                    # large enough that the front + tails stay off the edges
const CUTOFF           = 1e-12
const SLOPE_HALFWIDTH  = 3

# ---- (A) melt -> z  (DMT, trace-ful density) ----------------------------------------------------
const WALL             = NSITES ÷ 2           # energy domain wall on bond WALL (sites WALL, WALL+1)
const BETA             = 0.1                   # weak imbalance -> linear-response (infinitesimal) wall
const DT               = 0.1                   # per-gate time; one DMT evolve! advances t by 2*DT
const MELT_TMAX        = 16.0
const MELT_NCALL       = round(Int, MELT_TMAX / (2 * DT))
const MELT_MAXDIM      = 48
const MELT_GATE_MAXDIM = 4 * MELT_MAXDIM
const MELT_CONNECTOR_BUFFER = 8
const GIBBS_NSTEPS     = 6
const MELT_FIT_WINDOW  = (6.0, 16.0)
const EDGE_TOL         = 0.05

# ---- (B) spread -> D  (operator-space TEBD, with DAOE dissipation) ------------------------------
const CBOND            = NSITES ÷ 2           # central bond carrying the energy density h_c(0)
const MSD_DT           = 0.1                   # one Strang evolve! advances t by MSD_DT
const MSD_TMAX         = 16.0
const MSD_NSTEP        = round(Int, MSD_TMAX / MSD_DT)
const MSD_MAXDIM       = 96                    # DAOE caps operator entanglement, so a modest bond suffices
const MSD_FIT_WINDOW   = (8.0, 16.0)           # linear Var(t) ~ 2Dt window (after fast modes decay)
# DAOE dissipation: damp Pauli strings longer than LSTAR by GAMMA (per step; gamma_eff = GAMMA/MSD_DT
# per unit time). This tames the operator-entanglement growth that defeats *pure* TEBD on this chaotic
# model (which leaks O(5%) of the conserved energy by t~10). The extracted D is biased low at finite
# GAMMA / finite LSTAR; the unbiased value is the GAMMA -> 0, LSTAR -> infinity, t -> infinity limit.
# Set MSD_DAOE_GAMMA = 0.0 to recover pure TEBD (faithful but bond-starved -- see header).
const MSD_DAOE_LSTAR   = 3
const MSD_DAOE_GAMMA   = 0.04

const MELT_CSV = joinpath(@__DIR__, "mixed_field_ising_melt.csv")
const MSD_CSV  = joinpath(@__DIR__, "mixed_field_ising_diffusion.csv")
# -------------------------------------------------------------------------------------------------

# Mixed-field Ising bond Hamiltonian in the Pauli basis, with open-boundary field splitting so the
# sum over bonds reproduces H = J sum Z_iZ_{i+1} + g_x sum X_i + g_z sum Z_i exactly (bulk sites get
# half their field from each adjacent bond; edge sites get their full field from the single bond).
function mfi_bond_hamiltonian(nsites::Integer, bond::Integer; J::Real=J, gx::Real=GX, gz::Real=GZ)
    1 <= bond < nsites || throw(ArgumentError("bond index must satisfy 1 <= bond < nsites"))
    p = pauli_matrices()
    lw = bond == 1 ? 1.0 : 0.5
    rw = bond == nsites - 1 ? 1.0 : 0.5
    zz = J * kron(p.Z, p.Z)
    fx = gx * (lw * kron(p.X, p.I) + rw * kron(p.I, p.X))
    fz = gz * (lw * kron(p.Z, p.I) + rw * kron(p.I, p.Z))
    return zz + fx + fz
end

# Local energy densities as (start, dense) pairs -- the gate generators, the melt energy profile,
# and the spread autocorrelation targets all use this one consistent decomposition of H.
bond_terms(nsites::Integer; kwargs...) =
    [(b, mfi_bond_hamiltonian(nsites, b; kwargs...)) for b in 1:(nsites - 1)]

# Running log-log slope d ln(y)/d ln(t) via a sliding least-squares fit. For y ~ t^p the slope -> p;
# dE ~ t^(1/z), so the caller inverts z = 1/slope. (Verbatim from domain_wall_melting.jl.)
function running_loglog_slope(t::AbstractVector, y::AbstractVector; halfwidth::Int=SLOPE_HALFWIDTH)
    n = length(t)
    slope = fill(NaN, n)
    for k in 1:n
        lo, hi = k - halfwidth, k + halfwidth
        (lo >= 1 && hi <= n) || continue
        idx = [i for i in lo:hi if t[i] > 0 && y[i] > 0]
        length(idx) >= 2 || continue
        lt = log.(t[idx]); ly = log.(y[idx])
        m = length(lt); sx = sum(lt); sy = sum(ly)
        denom = m * sum(abs2, lt) - sx^2
        denom == 0 && continue
        slope[k] = (m * sum(lt .* ly) - sx * sy) / denom
    end
    return slope
end

# Running linear slope d y / d t via a sliding least-squares fit. For the spread variance
# Var ~ 2 D t, this returns 2 D(t); the caller halves it. (Linear sibling of the log-log slope.)
function running_linear_slope(t::AbstractVector, y::AbstractVector; halfwidth::Int=SLOPE_HALFWIDTH)
    n = length(t)
    slope = fill(NaN, n)
    for k in 1:n
        lo, hi = k - halfwidth, k + halfwidth
        (lo >= 1 && hi <= n) || continue
        idx = collect(lo:hi)
        tt = t[idx]; yy = y[idx]
        m = length(tt); sx = sum(tt); sy = sum(yy)
        denom = m * sum(abs2, tt) - sx^2
        denom == 0 && continue
        slope[k] = (m * sum(tt .* yy) - sx * sy) / denom
    end
    return slope
end

# Cheap self-checks: exact power law / exact line must recover their slope to ~machine precision.
function _assert_slope_oracles()
    t = collect(1.0:0.5:12.0)
    s = running_loglog_slope(t, 3.0 .* t .^ (1 / 2); halfwidth=3)
    sf = filter(isfinite, s)
    @assert !isempty(sf) && maximum(abs.(sf .- 1 / 2)) < 1e-9 "log-log slope oracle failed"
    l = running_linear_slope(t, 2.89 .* t .+ 1.7; halfwidth=3)
    lf = filter(isfinite, l)
    @assert !isempty(lf) && maximum(abs.(lf .- 2.89)) < 1e-9 "linear slope oracle failed"
    return nothing
end

# Least-squares slope over an explicit window (skips early transients and edge-contaminated times).
# For loglog we fit log|y| vs log t, so the window admits any nonzero y (dE is signed/negative).
function _window_slope(t, y, edge, (tmin, tmax); loglog::Bool)
    idx = findall(i -> tmin <= t[i] <= tmax && edge[i] < EDGE_TOL && (!loglog || (t[i] > 0 && y[i] != 0)),
                  eachindex(t))
    length(idx) >= 2 || return (NaN, length(idx))
    xs = loglog ? log.(t[idx]) : collect(float.(t[idx]))
    ys = loglog ? log.(abs.(y[idx])) : collect(float.(y[idx]))
    n = length(xs); sx = sum(xs)
    slope = (n * sum(xs .* ys) - sx * sum(ys)) / (n * sum(abs2, xs) - sx^2)
    return slope, n
end

# ===== (A) melt -> z : DMT on the trace-ful energy domain wall ====================================

function _build_melt(; nsites=NSITES, wall=WALL, beta=BETA, dt=DT, maxdim=MELT_MAXDIM,
                     gate_maxdim=MELT_GATE_MAXDIM, connector_buffer=MELT_CONNECTOR_BUFFER,
                     cutoff=CUTOFF, gibbs_nsteps=GIBBS_NSTEPS)
    sites = pauli_siteinds(nsites)
    terms = bond_terms(nsites)
    # +beta for bonds entirely left of the wall, -beta entirely right, 0 for the straddling bond.
    wall_weight(b) = b < wall ? beta : (b > wall ? -beta : 0.0)
    weights = [wall_weight(b) for (b, _) in terms]
    rho = pauli_gibbs_state(sites, terms, weights; nsteps=gibbs_nsteps, maxdim=2 * maxdim, cutoff=cutoff)
    gates = [pauli_gate_from_hamiltonian(h, dt) for (_, h) in terms]
    schedule = [b for (b, _) in terms]
    evo = DMTGateEvolution(gates, dt; schedule=schedule, reverse_schedule=reverse(schedule),
        maxdim=maxdim, cutoff=cutoff, gate_maxdim=gate_maxdim, connector_buffer=connector_buffer)
    return rho, evo, terms
end

function run_melt(; nsites=NSITES, wall=WALL, dt=DT, ncall=MELT_NCALL, kwargs...)
    rho, evo, terms = _build_melt(; nsites=nsites, wall=wall, dt=dt, kwargs...)
    profile0 = real.(pauli_expectation_profile(rho, terms))
    e_right0 = sum(profile0[(wall + 1):end]); total0 = sum(profile0)
    times = [0.0]; dE = [0.0]; drift = [0.0]; edge = [0.0]
    for k in 1:ncall
        evolve!(rho, evo)
        profile = real.(pauli_expectation_profile(rho, terms))
        delta = profile - profile0
        push!(times, 2 * dt * k)
        push!(dE, sum(profile[(wall + 1):end]) - e_right0)            # dE(t) ~ t^(1/z)
        push!(drift, abs(sum(profile) - total0))                      # energy conservation bar
        push!(edge, (abs(delta[1]) + abs(delta[2]) + abs(delta[end - 1]) + abs(delta[end])) / maximum(abs, delta))
    end
    return (; times, dE, drift, edge)
end

# ===== (B) spread -> D : TEBD on the traceless central energy density =============================

# Vectorize the central bond energy density h_c into an operator-space MPS, by exact two-site Pauli
# decomposition of the SAME dense bond term used everywhere else (so q_c is bit-consistent with the
# gates and the measurement targets). The overall (sqrt 2)^N normalization is irrelevant: every
# spread observable here is a ratio in which it cancels, so we omit it.
function _central_energy_state(sites, cbond; J::Real=J, gx::Real=GX, gz::Real=GZ)
    nsites = length(sites)
    h = mfi_bond_hamiltonian(nsites, cbond; J=J, gx=gx, gz=gz)
    paulis = collect(values(pauli_matrices()))          # [I, X, Y, Z]
    labels_of = (:I, :X, :Y, :Z)
    components = MPS[]
    for a in 1:4, b in 1:4
        coeff = tr(kron(paulis[a], paulis[b]) * h) / 4
        abs(coeff) < 1e-13 && continue
        labels = fill(:I, nsites)
        labels[cbond] = labels_of[a]; labels[cbond + 1] = labels_of[b]
        push!(components, pauli_basis_state(sites, labels; coefficient=coeff))
    end
    return reduce(+, components)
end

function _build_spread(; nsites=NSITES, dt=MSD_DT, maxdim=MSD_MAXDIM, cutoff=CUTOFF)
    evo = tebd_strang_evolution(nsites, dt;
        local_hamiltonian=(bond, weight) -> weight * mfi_bond_hamiltonian(nsites, bond),
        map_hamiltonian=pauli_gate_from_hamiltonian, maxdim=maxdim, cutoff=cutoff)
    return evo
end

# Spatial variance of the energy-density autocorrelation C_b(t) = tr[h_c(t) h_b] (bond position b).
function _spread_variance(qc, terms, cbond)
    C = real.(pauli_expectation_profile(qc, terms; normalize=false))   # traceless -> normalize=false
    bonds = 1:length(C)
    M0 = sum(C)                                                        # = tr[h_c(t) H], conserved
    M1 = sum(b * C[b] for b in bonds)
    M2 = sum(b^2 * C[b] for b in bonds)
    var = M2 / M0 - (M1 / M0)^2
    amax = maximum(abs, C)
    edge = (abs(C[1]) + abs(C[2]) + abs(C[end - 1]) + abs(C[end])) / amax  # front-at-edge guard
    return var, M0, edge
end

function run_spread(; nsites=NSITES, cbond=CBOND, dt=MSD_DT, nstep=MSD_NSTEP, maxdim=MSD_MAXDIM,
                    cutoff=CUTOFF, daoe_lstar=MSD_DAOE_LSTAR, daoe_gamma=MSD_DAOE_GAMMA)
    sites = pauli_siteinds(nsites)
    terms = bond_terms(nsites)
    qc = _central_energy_state(sites, cbond)
    evo = _build_spread(; nsites=nsites, dt=dt, maxdim=maxdim, cutoff=cutoff)
    dissipator = daoe_gamma > 0 ? pauli_daoe_projector(sites; lstar=daoe_lstar, gamma=daoe_gamma) : nothing
    var0, m00, e0 = _spread_variance(qc, terms, cbond)
    times = [0.0]; var = [var0]; m0 = [m00]; edge = [e0]
    for k in 1:nstep
        evolve!(qc, evo)
        if dissipator !== nothing
            qc[:] = apply(dissipator, qc; maxdim=maxdim, cutoff=cutoff)
        end
        v, mm, ee = _spread_variance(qc, terms, cbond)
        push!(times, dt * k); push!(var, v); push!(m0, mm); push!(edge, ee)
    end
    m0_drift = abs.(m0 .- m0[1]) ./ abs(m0[1])                          # conserved-trace error bar
    return (; times, var, m0_drift, edge)
end

# ===== timing probes (warmed up: bond dims grow as operators spread) ==============================

function _probe(build, step!; warmup=6, calls=3)
    state, ctx = build()
    for _ in 1:warmup; step!(state, ctx); end
    t0 = time()
    for _ in 1:calls; step!(state, ctx); end
    return (time() - t0) / calls
end

melt_probe(; kwargs...) = _probe(
    () -> (rho_evo = _build_melt(; kwargs...); (rho_evo[1], rho_evo[2])),
    (rho, evo) -> evolve!(rho, evo))

function spread_probe(; nsites=NSITES, cbond=CBOND, kwargs...)
    return _probe(
        () -> (sites = pauli_siteinds(nsites);
               (_central_energy_state(sites, cbond), _build_spread(; nsites=nsites, kwargs...))),
        (qc, evo) -> evolve!(qc, evo))
end

# ===== driver =====================================================================================

function main(; protocol::Symbol=:both, probe::Bool=true)
    _assert_slope_oracles()
    println("Mixed-field Ising energy diffusion  (J=$J, g_x=$GX, g_z=$GZ)   targets: z=$Z_TARGET, D~=$D_TARGET")
    println("  nsites=$NSITES  cutoff=$CUTOFF")

    if protocol in (:melt, :both)
        println("\n[A] MELT -> z   (DMT on the trace-ful energy domain wall; maxdim=$MELT_MAXDIM, beta=$BETA, t_max=$MELT_TMAX)")
        probe && @printf("    timing probe (warmed): %.2f s / DMT sweep  ->  ~%.1f min\n",
                         (p = melt_probe()), MELT_NCALL * p / 60)
        m = run_melt()
        z = 1.0 ./ running_loglog_slope(m.times, abs.(m.dE))
        println("    ", rpad("t", 8), rpad("dE(t)", 14), rpad("z(t)", 9), rpad("energy drift", 14), "front@edges")
        for k in eachindex(m.times)
            (k == 1 || (k - 1) % 5 == 0) || continue
            zt = isfinite(z[k]) ? @sprintf("%.3f", z[k]) : "--"
            @printf("    %-8.2f%-14.4e%-9s%-14.2e%.2e\n", m.times[k], m.dE[k], zt, m.drift[k], m.edge[k])
        end
        inv_z, npz = _window_slope(m.times, m.dE, m.edge, MELT_FIT_WINDOW; loglog=true)
        @printf("    => z_eff = %.2f over t in %s (%d pts)   [target z = %.1f]\n",
                isnan(inv_z) ? NaN : 1 / inv_z, MELT_FIT_WINDOW, npz, Z_TARGET)
        open(MELT_CSV, "w") do io
            println(io, "t,dE,z,drift,edge")
            for k in eachindex(m.times)
                @printf(io, "%.4f,%.6e,%.4f,%.3e,%.3e\n", m.times[k], m.dE[k], z[k], m.drift[k], m.edge[k])
            end
        end
        println("    wrote $MELT_CSV")
    end

    if protocol in (:spread, :both)
        daoe = MSD_DAOE_GAMMA > 0 ? "DAOE lstar=$MSD_DAOE_LSTAR gamma=$MSD_DAOE_GAMMA" : "pure TEBD"
        println("\n[B] SPREAD -> D  (TEBD on the traceless central energy density; maxdim=$MSD_MAXDIM, $daoe, t_max=$MSD_TMAX)")
        probe && @printf("    timing probe (warmed): %.2f s / Strang sweep  ->  ~%.1f min\n",
                         (p = spread_probe()), MSD_NSTEP * p / 60)
        s = run_spread()
        Dt = running_linear_slope(s.times, s.var) ./ 2
        println("    ", rpad("t", 8), rpad("Var(t)", 14), rpad("D(t)", 9), rpad("M0 drift", 14), "front@edges")
        for k in eachindex(s.times)
            (k == 1 || (k - 1) % 10 == 0) || continue
            dt_ = isfinite(Dt[k]) ? @sprintf("%.3f", Dt[k]) : "--"
            @printf("    %-8.2f%-14.4e%-9s%-14.2e%.2e\n", s.times[k], s.var[k], dt_, s.m0_drift[k], s.edge[k])
        end
        slope, npd = _window_slope(s.times, s.var, s.edge, MSD_FIT_WINDOW; loglog=false)
        @printf("    => D = %.3f over t in %s (%d pts)   [target D ~= %.3f]\n",
                isnan(slope) ? NaN : slope / 2, MSD_FIT_WINDOW, npd, D_TARGET)
        open(MSD_CSV, "w") do io
            println(io, "t,var,D,m0_drift,edge")
            for k in eachindex(s.times)
                @printf(io, "%.4f,%.6e,%.4f,%.3e,%.3e\n", s.times[k], s.var[k], Dt[k], s.m0_drift[k], s.edge[k])
            end
        end
        println("    wrote $MSD_CSV")
    end

    println("\nz(t) and D(t) are effective, resource-bound: the diffusive plateau lives between the fast-mode")
    println("decay (t ~ a few/Omega) and the front reaching the edges. Raise NSITES / MAXDIM / T_MAX to sharpen;")
    println("the energy-drift and M0-drift columns are the truncation error bars.")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
