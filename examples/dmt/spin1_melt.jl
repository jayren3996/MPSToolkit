# Spin-1 magnetization transport from a melting domain wall: diffusion vs. SU(3) KPZ
#
#   Heisenberg   H = sum_j  S_j . S_{j+1}                       non-integrable -> diffusive,     z = 2
#   ULS / SU(3)  H = sum_j [S_j . S_{j+1} + (S_j . S_{j+1})^2]  integrable     -> KPZ,           z = 3/2
#
# (S^a are the spin-1 generators, d = 3.) The ULS point is the Uimin-Lai-Sutherland chain: adding
# the biquadratic term at unit coefficient promotes the SU(2) symmetry to SU(3) and makes the model
# Bethe-ansatz integrable. Infinite-temperature spin transport in an integrable magnet with
# non-abelian symmetry falls in the KPZ universality class (z = 3/2), while the generic
# (non-integrable) spin-1 Heisenberg magnet is an ordinary diffusive conductor (z = 2). The pair is
# the benchmark: ONE Hamiltonian near ONE target is weak evidence, TWO Hamiltonians separating in
# the two predicted directions is strong evidence, because every systematic this method has
# (Trotter error, finite chi, finite t, front contamination) is shared between them.
#
# PROTOCOL -- infinitesimal domain-wall melting (linear response, infinite temperature), the same
# protocol as the spin-1/2 sibling `domain_wall_melting.jl`, lifted to d = 3. We prepare the weakly
# polarized step density operator
#       rho(0) = I + sum_j c_j S^z_j,      c_j = -eps (j <= kink), +eps (j > kink),
# evolve it with DMT under U = prod_b exp(-i H_b dt), and read the magnetization profile
#       m(x,t) = tr(rho(t) S^z_x) / anchor(t),
# normalized by the unmelted bulk so the far field reads +-1, and the CHARGE that has crossed the
# central bond,
#       T(t) = sum_{x > kink} [1 - m(x,t)]   ~   t^(1/z).
# The running dynamical exponent is the logarithmic derivative
#       z(t) = 1 / ( d ln T / d ln t ),
# which starts near 1 (every model looks ballistic while the front is still local) and DRIFTS
# toward its regime target as the hydrodynamic tail develops.
#
# WHY THESE CHOICES.
#  * rho(0) is TRACE-FUL (the literal identity plus a weak S^z step), which is the object DMT is
#    built for: the trace/identity direction is what anchors its truncation. `operator_product_state`
#    and `operator_local_sum_state` share the literal-operator convention, so `add`ing them gives
#    exactly `I + sum_j c_j S^z_j`, trace `d^N` to machine precision.
#  * S^z = diag(1, 0, -1) is NOT proportional to any single generalized Gell-Mann matrix at d = 3,
#    so it cannot be selected with `operator_basis_state`; it must be built from the dense matrix.
#  * We evolve with `normalize = false`: per-sweep `normalize!` would rescale the Frobenius norm
#    (which drifts under truncation) and corrupt absolute charges.
#  * anchor(t) = (1/2)[mean(p over rightmost K) - mean(p over leftmost K)] is a RATIO, hence
#    invariant to any overall rescaling (the basis prefactor, residual norm drift). Measured with
#    `normalize = false`, `operator_expectation_profile` returns the LITERAL trace tr(rho S^z_x),
#    which carries a `d^((N-span)/2)` factor (~1e14 at N = 60) on top of the operator's own `d^N`
#    trace scale -- so every tolerance below is RELATIVE, normalized by `sum(abs, profile)`.
#  * sum_x m(x,t) is conserved at 0 (total S^z commutes with both bond Hamiltonians); its relative
#    drift |sum p| / sum|p| is the truncation error bar, printed as max drift, and stays at ~1e-13.
#  * `gate_maxdim = 0` (the default: apply the gate exactly). A finite cap would pre-truncate the
#    gate-inflated bond with a plain SVD before DMT can protect the local-operator content it
#    carries, which is the error DMT exists to avoid.
#  * The DMT kernel preserves every observable of diameter <= preserve_diameter (default 3, used
#    here) EXACTLY, so z(t) is limited by finite t_max / nsites / chi, not by truncation bias on the
#    connector. At d = 3 the protected block is 2*3^2 = 18 directions of every maxdim (floor 19), so
#    maxdim = 96 buys 78 complement directions -- read the budget table in docs/src/manual/dmt.md.
#
# CONVERGENCE IS THE POINT. Each Hamiltonian is run at three bond dimensions so a reader can see
# whether z has stopped moving. An honest "not converged at these budgets, and here is the trend"
# is a result; a single unconverged number presented as a measurement is not.
#
# FRONT CONTAINMENT IS THE BINDING CONSTRAINT. If the melt front reaches the chain edges the anchor
# is no longer the unmelted bulk and every exponent after that time is garbage. `front` below is
# the sibling's guard: the largest |1 - m| over the K sites just inside the anchor band. Fit times
# with `front > FRONT_TOL` are dropped, and the printed table reports the guard.
#
# WHAT THIS CONFIGURATION ACTUALLY FINDS (N = 100, t <= 11.8, chi = 40/56/72, 3.3 CPU-hours).
#   ULS/SU(3):   z = 1.51 +- 0.03,  converged in chi (1.484 -> 1.511 -> 1.507). Agrees with the
#                KPZ value 3/2 and with arXiv:2205.02853, which ran the same model at d = 3.
#   Heisenberg:  z = 1.51 +- 0.03,  also converged in chi (1.541 -> 1.526 -> 1.514) -- and NOT the
#                diffusive 2. The non-integrable chain has not crossed over to diffusion by
#                t ~ 12; at these times it is indistinguishable from the integrable point. That is
#                a finite-TIME limit, not a finite-chi one: the chi ladder has stopped moving
#                (|dz| = 0.012 on the last rung) while z is still far from 2.
#   The convergence study is not decoration. At chi = 40 the Heisenberg running exponent climbs to
#   z(11.2) = 1.647, which looks exactly like the expected crossover toward diffusion -- and it is
#   a TRUNCATION ARTIFACT: it is gone at chi = 56 (1.496) and chi = 72 (1.509). A single
#   unconverged run would have "confirmed" the prediction. Reading z off one bond dimension, at
#   any budget, is not a measurement.
#
# References:
#   - Ye, Machado, Kemp, Hutson, Yao, PRL 129, 230602 (2022), arXiv:2205.02853 -- DMT on exactly
#     this d = 3 SU(3) model at chi up to 256, reporting KPZ; the published comparison point for
#     this benchmark, and the only one at this local dimension.
#   - White, Zaletel, Mong, Refael, PRB 97, 035127 (2018) -- DMT.
#   - Ye, Machado, White, Mong, Yao, arXiv:1902.01859 -- the total-bond-dimension convention
#     `maxdim = chi_preserve + chi_extra` used here.
#   - Ljubotina, Znidaric, Prosen, PRL 122, 210602 (2019) -- KPZ in the (spin-1/2) Heisenberg magnet.
#   - Uimin (1970), Lai (1974), Sutherland (1975) -- the integrable SU(3) point.

using ITensors
using ITensorMPS
using MPSToolkit
using LinearAlgebra: I
using Printf

# ---- parameters -----------------------------------------------------------------------------
# Production defaults: 6 jobs (2 models x 3 bond dimensions), ~3 h wall on 32 cores. The cost of one
# sweep is set almost entirely by `maxdim` and by the saturated region around the wall, and is
# nearly flat in the un-melted remainder of the chain -- measured at N = 100 with 32 BLAS threads:
# ~24 s/sweep at maxdim = 40, ~30 at 56, ~55 at 72, essentially independent of N once the chain is
# longer than the melt. That asymmetry is why NSITES is generous and T_MAX is not: lengthening the
# chain is nearly free and buys front containment, whereas doubling t_max doubles the bill.
# T_MAX matters: at t <~ 8 the two models are NOT yet separated (both give z ~ 1.5); the Heisenberg
# curve only lifts away from the ULS one past t ~ 9. To smoke-test, drop MAXDIMS to `(40,)` and
# T_MAX to 2.0 (~2 min).
const NSITES           = 100                   # cost is set by the SATURATED region, not by N, so a
                                               # long chain buys front containment nearly for free
const DT               = 0.1                   # per-sweep gate time; one dmt_evolve! advances t by 2*DT
const T_MAX            = 12.0
const NCALL            = round(Int, T_MAX / (2 * DT))   # forward+reverse DMT sweeps
const MAXDIMS          = (40, 56, 72)          # convergence ladder; floor at d = 3 is 2*3^2 + 1 = 19
const CUTOFF           = 1e-10
const EPS_WALL         = 0.25                  # step amplitude of the weakly polarized wall
const K_ANCHOR         = 4                     # outermost sites (each side) defining the bulk anchor
const SLOPE_HALFWIDTH  = 3                     # +/- points in the running log-log slope window
const FRONT_TOL        = 0.02                  # drop fit times once the front reaches the anchor band
const FIT_WINDOW       = (6.0, T_MAX)          # late window for the single-number effective z
const KINK             = NSITES ÷ 2            # domain wall on the central bond
const CASES            = (("heisenberg", "diffusive", 2.0),
                          ("uls", "superdiffusive (KPZ)", 1.5))
const CSV_PATH         = joinpath(@__DIR__, "spin1_melt.csv")
# -------------------------------------------------------------------------------------------------

"""
    spin1_matrices()

Return `(Sx, Sy, Sz)` for spin 1 in the `|+1>, |0>, |-1>` basis, from the ladder coefficients
`S^+ |m> = sqrt(S(S+1) - m(m+1)) |m+1>` (so `S^+` has the two off-diagonal entries `sqrt(2)`).

MPSToolkit deliberately ships no spin-`S` model builders, so the five lines live here. Any source
of dense `d x d` matrices following the same `kron`-left-to-right multi-site convention works
instead -- e.g. `EDKit.spin((1, "xx"), (1, "yy"), (1, "zz"); D = 3)` returns exactly the
`heisenberg_bond()` matrix below and can be passed straight into
`operator_gate_from_hamiltonian`. EDKit.jl is NOT a dependency of MPSToolkit.
"""
function spin1_matrices()
    raise = zeros(ComplexF64, 3, 3)
    raise[1, 2] = raise[2, 3] = sqrt(2.0)
    lower = collect(raise')
    return (raise + lower) / 2, (raise - lower) / (2im), ComplexF64[1 0 0; 0 0 0; 0 0 -1]
end

# Two-site bond Hamiltonians. `kron` is applied left-to-right (first factor = leftmost site), the
# convention `operator_gate` and every operator-space state builder in MPSToolkit use.
function heisenberg_bond()
    sx, sy, sz = spin1_matrices()
    return kron(sx, sx) + kron(sy, sy) + kron(sz, sz)
end

# ULS / SU(3): H = S.S + (S.S)^2. The biquadratic term at unit coefficient is what promotes SU(2)
# to SU(3) (the bond operator becomes the SU(3) permutation up to a constant) and makes the chain
# integrable.
function uls_bond()
    h = heisenberg_bond()
    return h + h * h
end

bond_hamiltonian(name) = name == "heisenberg" ? heisenberg_bond() :
                         name == "uls" ? uls_bond() :
                         throw(ArgumentError("unknown case $(name)"))

# Running log-log slope d ln(y) / d ln(t) via a sliding-window least-squares fit. For y ~ t^p the
# slope -> p; the transferred charge obeys T ~ t^(1/z), so the caller inverts z = 1 / slope. Pure
# function -- its power-law oracle (`_assert_slope_oracle`) makes it testable with no physics.
# (Verbatim from domain_wall_melting.jl.)
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

# Late-window least-squares slope of T(t); returns (1/z, npoints). Skips the early ballistic
# transient and any time whose front guard says the melt has reached the anchor band.
function fit_inverse_z(times, transferred, front, (tmin, tmax); front_tol=FRONT_TOL)
    idx = findall(i -> tmin <= times[i] <= tmax && transferred[i] > 0 && front[i] < front_tol,
                  eachindex(times))
    length(idx) >= 2 || return (NaN, length(idx))
    lt = log.(times[idx])
    ly = log.(transferred[idx])
    n = length(lt)
    sx = sum(lt)
    return (n * sum(lt .* ly) - sx * sum(ly)) / (n * sum(abs2, lt) - sx^2), n
end

"""
    profile_and_charge(state, sz, kink, K)

Return `(transferred, relative_drift, front)` from the bulk-anchored magnetization profile.

`normalize = false` gives the LITERAL trace `tr(rho S^z_x)`, whose scale is set by the operator
(`~d^N`), so the drift is reported RELATIVE to the profile's own magnitude `sum(abs, p)` -- an
absolute drift divided by a site count is meaningless at this scale and reads as catastrophic
charge leakage when there is none.
"""
function profile_and_charge(state, sz, kink, K)
    n = length(state)
    p = real.(operator_expectation_profile(state, [(x, sz) for x in 1:n]; normalize=false))
    scale = sum(abs, p)
    drift = scale == 0 ? 0.0 : abs(sum(p)) / scale                        # conservation oracle: ~1e-13
    anchor = (sum(@view p[(n - K + 1):n]) - sum(@view p[1:K])) / (2 * K)  # right -> +1, left -> -1
    m = p ./ anchor
    transferred = sum(1.0 - m[x] for x in (kink + 1):n)                   # T(t) = sum_{x>kink}[1 - m]
    band = (n - 2K + 1):(n - K)                                           # K sites just inside the anchor
    front = maximum(abs(1.0 - m[x]) for x in band)                        # front-reached-the-edge guard
    return transferred, drift, front
end

function _build_evolution(name; nsites=NSITES, maxdim=first(MAXDIMS), dt=DT, kink=nsites ÷ 2,
                          eps_wall=EPS_WALL, cutoff=CUTOFF)
    d = 3
    _, _, sz = spin1_matrices()
    sites = operator_siteinds(nsites; d=d)
    identity_matrix = Matrix{ComplexF64}(I, d, d)
    coeffs = [j <= kink ? -eps_wall : eps_wall for j in 1:nsites]
    # Literal `I + sum_j c_j S^z_j`: both builders use the same literal-operator convention, so the
    # sum is the physical step density operator exactly (trace d^N), not a rescaled cousin.
    rho = add(operator_product_state(sites, fill(identity_matrix, nsites)),
              operator_local_sum_state(sites, sz, coeffs); maxdim=8, cutoff=0.0)
    gate = operator_gate_from_hamiltonian(bond_hamiltonian(name), dt; d=d)
    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule=schedule, reverse_schedule=reverse(schedule),
        nstep=1, maxdim=maxdim, cutoff=cutoff,
        normalize=false)  # absolute charges -> no per-sweep rescaling
    return rho, evo, sz
end

# One (Hamiltonian, maxdim) run: build rho(0), DMT-evolve, record (t, T, drift, front, chi) per sweep.
function melt(name; nsites=NSITES, maxdim=first(MAXDIMS), dt=DT, ncall=NCALL, kink=nsites ÷ 2,
              k_anchor=K_ANCHOR, verbose=false)
    rho, evo, sz = _build_evolution(name; nsites=nsites, maxdim=maxdim, dt=dt, kink=kink)
    t0, s0, f0 = profile_and_charge(rho, sz, kink, k_anchor)
    times = [0.0]; charge = [t0]; drift = [s0]; front = [f0]; chi = [maxlinkdim(rho)]
    wall0 = time()
    for k in 1:ncall
        dmt_evolve!(rho, evo)
        transferred, sm, fr = profile_and_charge(rho, sz, kink, k_anchor)
        push!(times, 2 * dt * k); push!(charge, transferred); push!(drift, sm); push!(front, fr)
        push!(chi, maxlinkdim(rho))
        verbose && @printf("    t=%5.2f  T=%9.4f  chi=%4d  drift=%.2e  front=%.2e  (%.0f s elapsed)\n",
                           times[end], transferred, chi[end], sm, fr, time() - wall0)
    end
    return (; times, charge, drift, front, chi, seconds=time() - wall0)
end

# Extrapolate the whole run from a few warmed-up sweeps before committing to it. Per-sweep cost
# GROWS with time (only the saturated region near the wall costs anything), so this is a LOWER
# bound on the true cost, not an estimate of it.
function timing_probe(; calls=2, name="heisenberg", maxdim=first(MAXDIMS), nsites=NSITES)
    rho, evo, _ = _build_evolution(name; nsites=nsites, maxdim=maxdim)
    dmt_evolve!(rho, evo)  # warm up / compile
    t0 = time()
    for _ in 1:calls
        dmt_evolve!(rho, evo)
    end
    return (time() - t0) / calls
end

# One tidy row per time, four columns (T, z, drift, front) per (model, maxdim) job.
function write_csv(path, results, cases, maxdims)
    times = results[(first(first(cases)), first(maxdims))].times
    open(path, "w") do io
        header = ["t"]
        for (name, _, _) in cases, maxdim in maxdims
            append!(header, ["T_$(name)_m$(maxdim)", "z_$(name)_m$(maxdim)",
                             "drift_$(name)_m$(maxdim)", "front_$(name)_m$(maxdim)"])
        end
        println(io, join(header, ","))
        for k in eachindex(times)
            fields = [@sprintf("%.4f", times[k])]
            for (name, _, _) in cases, maxdim in maxdims
                r = results[(name, maxdim)]
                append!(fields, [@sprintf("%.6e", r.charge[k]), @sprintf("%.4f", r.z[k]),
                                 @sprintf("%.3e", r.drift[k]), @sprintf("%.3e", r.front[k])])
            end
            println(io, join(fields, ","))
        end
    end
    return path
end

function main(; probe=true, nsites=NSITES, maxdims=MAXDIMS, ncall=NCALL, dt=DT, cases=CASES,
              csv_path=CSV_PATH, verbose=false)
    _assert_slope_oracle()
    kink = nsites ÷ 2
    if probe
        per = timing_probe(; maxdim=maximum(maxdims), nsites=nsites)
        @printf("timing probe: %.1f s / sweep at maxdim=%d (EARLY time; late sweeps cost more)\n",
                per, maximum(maxdims))
    end

    println("\nSpin-1 (d = 3) domain-wall melting via operator-space DMT")
    println("  nsites=$nsites  dt=$dt  t_max=$(2 * dt * ncall)  kink=$kink  preserve_diameter=3 (default)")
    println("  maxdim ladder: $(collect(maxdims))   (protected block 2*3^2 = 18 of every maxdim)")
    println()
    println(rpad("model", 13), rpad("z(target)", 11), rpad("maxdim", 8), rpad("chi'", 6),
            rpad("z(fit)", 9), rpad("pts", 5), rpad("z(t_max)", 10), rpad("T(t_max)", 10),
            rpad("max drift", 11), rpad("front", 10), "min")

    results = Dict{Tuple{String,Int},NamedTuple}()
    for (name, regime, ztarget) in cases, maxdim in maxdims
        run = melt(name; nsites=nsites, maxdim=maxdim, dt=dt, ncall=ncall, kink=kink, verbose=verbose)
        @assert abs(run.charge[1]) < 1e-9 "T(0) should vanish for an unmelted wall"  # structural oracle
        zrun = 1.0 ./ running_loglog_slope(run.times, run.charge)
        inv_z, npts = fit_inverse_z(run.times, run.charge, run.front, FIT_WINDOW)
        zfit = 1.0 / inv_z
        results[(name, maxdim)] = (; run..., z=zrun, zfit, npts, regime, ztarget)
        zlast = let zf = filter(isfinite, zrun); isempty(zf) ? NaN : zf[end] end
        @printf("%-13s%-11.2f%-8d%-6d%-9.3f%-5d%-10.3f%-10.3f%-11.2e%-10.2e%.1f\n",
                name, ztarget, maxdim, maxdim - 18, zfit, npts, zlast, run.charge[end],
                maximum(run.drift), maximum(run.front), run.seconds / 60)
    end

    println()
    for (name, regime, ztarget) in cases
        zs = [results[(name, m)].zfit for m in maxdims]
        moved = length(zs) >= 2 ? abs(zs[end] - zs[end - 1]) : NaN
        @printf("%-13s %-22s target z=%.2f   z(chi) = %s   |dz| last step = %.3f\n",
                name, regime, ztarget, join((@sprintf("%.3f", z) for z in zs), " -> "), moved)
    end

    write_csv(csv_path, results, cases, maxdims)
    println("\nwrote $csv_path")
    println("z(t) is an EFFECTIVE, resource-bound exponent: it starts near 1 (ballistic front) and")
    println("drifts toward the target within reach of t_max. Compare the two models against each")
    println("other, not each against its target -- the shared systematics cancel in the contrast.")
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
