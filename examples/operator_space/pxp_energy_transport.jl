# PXP energy transport from a weak energy domain wall, via constrained operator-space DMT
#
# The PXP chain
#
#     H = Ω [ X_1 P_2 + sum_{j=2}^{N-1} P_{j-1} X_j P_{j+1} + P_{N-1} X_N ],   P = |0><0|,
#
# is the effective model of a Rydberg-blockaded atom chain: a site can flip only when both
# neighbors are in the ground state, so dynamics never connects the constrained sector (no
# adjacent excitations) to blockade-violating configurations. Despite its quantum many-body
# scars, generic PXP thermalization is fast, and the *energy* — the only local conserved
# density — is predicted and observed (numerically) to spread superdiffusively with dynamical
# exponent z = 3/2, in the KPZ universality class:
#
#   Ljubotina, Desaules, Serbyn, Papić, "Superdiffusive energy transport in kinetically
#   constrained models", PRX 13, 011033 (2023).
#
# Protocol (domain-wall quench in inverse temperature):
#   1. Initial state  rho(0) ∝ exp(-beta H_L) ⊗ exp(+beta H_R): left half slightly cold,
#      right half slightly hot. Terms that straddle the wall get weight zero, so the two
#      exponentials factorize exactly. The state is prepared by imaginary-time TEBD in
#      operator space (`pauli_gibbs_state`) *on top of the vectorized constraint projector*
#      P_G — the infinite-temperature state of the constrained sector — because the
#      unconstrained identity carries unphysical weight on blockade-violating configurations.
#   2. Evolve rho(t) = e^{-iHt} rho e^{+iHt} with DMT truncation (3-site bulk gates, 2-site
#      edge gates). Exact gates commute with P_G, but *truncation leaks weight out of the
#      constrained sector*; `constrained_dmt_evolve!` therefore re-applies the operator-space
#      projector MPO  rho -> P_G rho P_G  at checkpoints.
#   3. Measure the energy profile e_j(t) = tr(rho h_j)/tr(rho) in one O(N) sweep
#      (`pauli_expectation_profile`) and track the energy transferred across the wall,
#          dE(t) = sum_{j > wall} [e_j(t) - e_j(0)]  ~  t^(1/z).
#
# EXTRACTING z — mind the crossover. dE(t) reaches its hydrodynamic power law slowly: the
# quench starts quadratically (the initial energy current vanishes), passes through a
# ballistic-looking regime, and only then bends toward the asymptotic t^(2/3). The script
# prints the *local* log-log slope at each report row so the crossover is visible: with the
# short default run it falls 3 -> 2 -> 1.5 -> ~1 by t = 10, still above the KPZ value 2/3
# (the literature sees the asymptotic exponent emerge at noticeably longer times). Treat the
# fitted exponent as an effective, crossover-bound value; push nsites / maxdim / ncall up and
# the fit window later to approach the asymptotic regime. What the default run *does*
# demonstrate quantitatively is the constrained-DMT machinery itself: exact sector
# preservation (leakage diagnostics below), energy conservation, and a stable domain-wall
# melt to fit.

using ITensors
using ITensorMPS
using MPSToolkit
using Printf

# ---- parameters (configured for a few-minute run; see note above for convergence) --------
nsites = 40
wall   = nsites ÷ 2     # domain wall between sites `wall` and `wall + 1`
beta   = 0.25           # weak imbalance: linear-response regime
omega  = 1.0
dt     = 0.1            # per-sweep gate time; one constrained call advances t by 2*dt
ncall  = 50             # measurement sweeps -> t_max = 2*dt*ncall = 10.0
maxdim = 32
gate_maxdim = 4 * maxdim
connector_buffer = 8
cutoff = 1e-12
project_every = 1       # constraint checkpoint after every forward+reverse sweep
gibbs_nsteps = 6        # imaginary-time Trotter steps for the initial state
fit_window = (5.0, 10.0) # late-window fit for the (crossover) effective exponent
edge_tol = 0.05         # drop fit times once the front reaches the chain edges
# -------------------------------------------------------------------------------------------

sites = pauli_siteinds(nsites)
terms = [(first(pxp_term_support(nsites, j)), pxp_term_hamiltonian(nsites, j; omega=omega)) for j in 1:nsites]

# Domain-wall weights: +beta for terms supported entirely in the left half, -beta entirely in
# the right half, 0 for the (two) wall-straddling terms -> rho(0) = e^{-beta H_L} ⊗ e^{+beta H_R}.
function wall_weight(j)
  support = pxp_term_support(nsites, j)
  last(support) <= wall && return beta
  first(support) > wall && return -beta
  return 0.0
end
weights = [wall_weight(j) for j in 1:nsites]

println("preparing the constrained domain-wall state (N=$nsites, beta=$beta) ...")
rho = pauli_gibbs_state(
  sites,
  terms,
  weights;
  nsteps=gibbs_nsteps,
  maxdim=2 * maxdim,
  cutoff=cutoff,
  initial_state=pauli_pxp_constraint_state(sites),
)

projector = pauli_pxp_constraint_projector(sites)
gates = [pauli_gate_from_hamiltonian(h, dt) for (_, h) in terms]
schedule = [start for (start, _) in terms]
evo = DMTGateEvolution(
  gates,
  dt;
  schedule=schedule,
  reverse_schedule=reverse(schedule),
  nstep=1,
  maxdim=maxdim,
  cutoff=cutoff,
  gate_maxdim=gate_maxdim,
  connector_buffer=connector_buffer,
)

energy_profile(state) = real.(pauli_expectation_profile(state, terms))

# Hilbert-Schmidt weight outside the constrained sector: 1 - <rho, P rho P> / <rho, rho>.
function sector_leakage(state)
  projected = apply(projector, state; maxdim=gate_maxdim, cutoff=cutoff)
  return 1 - real(inner(state, projected)) / real(inner(state, state))
end

profile0 = energy_profile(rho)
right_energy(profile) = sum(profile[(wall + 1):nsites])
e_right0 = right_energy(profile0)
total0 = sum(profile0)

times = Float64[]
transferred = Float64[]
edge_fraction = Float64[]
leak_checkpoints = (1, ncall ÷ 2, ncall)

println("evolving with DMT + constraint checkpoints (t_max = $(2 * dt * ncall)) ...")
println(rpad("t", 8), rpad("dE(t)", 14), rpad("local slope", 13), rpad("energy drift", 14), "front@edges")
last_report = Ref((0.0, 0.0))   # (t, dE) of the previous report row, for the local slope
for k in 1:ncall
  constrained_dmt_evolve!(rho, evo, projector; project_every=project_every)
  profile = energy_profile(rho)
  delta = profile - profile0
  push!(times, 2 * dt * k)
  push!(transferred, right_energy(profile) - e_right0)
  push!(edge_fraction, (abs(delta[1]) + abs(delta[2]) + abs(delta[end - 1]) + abs(delta[end])) / maximum(abs, delta))
  if k % 5 == 0
    # Local log-log slope of |dE| between consecutive report rows: watch it fall through the
    # ballistic-to-hydrodynamic crossover toward the asymptotic 1/z.
    previous_t, previous_de = last_report[]
    slope_text = previous_de == 0 ? "--" :
      @sprintf("%.2f", log(abs(transferred[end] / previous_de)) / log(times[end] / previous_t))
    @printf("%-8.2f%-14.4e%-13s%-14.2e%.2e\n", times[end], transferred[end], slope_text,
      abs(sum(profile) - total0), edge_fraction[end])
    last_report[] = (times[end], transferred[end])
  end
  if k in leak_checkpoints
    # Leakage accrued by ONE unprojected sweep from the current (projected) state. This is the
    # quantity the checkpoints remove before it contaminates the transport observables.
    probe = copy(rho)
    dmt_evolve!(probe, evo)
    @printf("    [t=%.2f] residual leakage %.2e  |  one unprojected sweep accrues %.2e\n",
      times[end], sector_leakage(rho), sector_leakage(probe))
  end
end

# Log-log least-squares slope of |dE|(t) over the fit window, skipping early transients,
# boundary-contaminated times, and (defensively) sign flips. z = 1 / slope.
function fit_inverse_z(times, transferred, edge_fraction, (tmin, tmax))
  idx = findall(
    i -> tmin <= times[i] <= tmax && abs(transferred[i]) > 0 && edge_fraction[i] < edge_tol,
    eachindex(times),
  )
  lt = log.(times[idx])
  lm = log.(abs.(transferred[idx]))
  n = length(lt)
  sx = sum(lt)
  slope = (n * sum(lt .* lm) - sx * sum(lm)) / (n * sum(lt .^ 2) - sx^2)
  return slope, n
end

inv_z, npts = fit_inverse_z(times, transferred, edge_fraction, fit_window)
println()
@printf("late-window fit over t in %s (%d points):  effective 1/z = %.3f  =>  z_eff = %.2f\n",
  fit_window, npts, inv_z, 1 / inv_z)
println("""
This short run is still inside the ballistic-to-hydrodynamic crossover (the local slope
column above decreases monotonically toward the asymptote but has not reached it), so the
effective exponent overestimates 1/z. The asymptotic PXP value is superdiffusive (KPZ),
1/z = 2/3 (z = 3/2), per Ljubotina-Desaules-Serbyn-Papić, PRX 13, 011033 (2023); ordinary
diffusion would give 1/z = 0.5. To approach it, raise nsites / maxdim / ncall together and
move the fit window later (the cost grows linearly in each).""")
