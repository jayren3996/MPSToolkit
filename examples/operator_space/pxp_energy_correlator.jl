# PXP energy-transport exponent from the infinite-temperature energy correlator,
# via constrained operator-space DMT
#
# This is the COMPANION script to examples/dmt/pxp_energy_melting.jl, and the more efficient route to
# the dynamical exponent z. Instead of melting a domain wall, it Heisenberg-evolves a single
# local energy density in the constrained sector and watches it spread — the protocol of
#
#   Ljubotina, Desaules, Serbyn, Papić, "Superdiffusive energy transport in kinetically
#   constrained models", PRX 13, 011033 (2023),
#
# which finds KPZ-class superdiffusion, z = 3/2, for PXP energy transport.
#
# Protocol. The constrained infinite-temperature correlator
#
#     C(x, t) = tr( P_G h_x O(t) ),     O(0) ∝ P_G h_c P_G,   O(t) = e^{-iHt} O e^{+iHt},
#
# is the coefficient profile of the evolved (traceless) energy density. Its second moment
#
#     M2(t) = sum_x (x - c)^2 C(x, t) / sum_x C(x, t)  ~  t^(2/z)
#
# has asymptotic log-log slope 2/z = 4/3 for z = 3/2. The setup needs the same three
# constraint-aware ingredients as the domain-wall script: the vectorized constraint
# projector as the sector's reference state, a sector-projected initial operator, and
# constraint-projection checkpoints during DMT evolution (truncation leaks weight out of
# the sector even though every PXP term commutes with P_G).
#
# TWO THINGS THIS SCRIPT DOES DIFFERENTLY (both matter for traceless operators):
#   * `normalize=false` everywhere. DMT truncation sheds Hilbert-Schmidt norm faster than it
#     sheds the protected components, so the default per-sweep renormalization silently
#     inflates absolute traces over time (tens of percent over a long run) while leaving
#     ratio observables like M2 untouched. With `normalize=false` the conserved total
#     sum_x C(x, t) = tr(P_G H O(t)) becomes a real diagnostic: its drift IS the
#     truncation error bar on the run.
#   * Unnormalized profiles. O is traceless, so `pauli_expectation_profile` must be called
#     with `normalize=false` (the trace it would divide by vanishes).
#
# READING OFF z — crossover honesty. M2's local slope descends toward 4/3 from above
# through a slow crossover. Production-validated reference points (N=64, t<=14):
#     chi=48: local slope plateaus at ~1.46-1.49 over t = 8-14, still descending — the
#             converged crossover value at these times;
#     chi=32: plateaus lower, ~1.32-1.38 — truncation suppresses M2 growth, so the slope
#             UNDERSHOOTS the converged value (landing near 4/3 by partial cancellation,
#             not by having reached the asymptote).
# So: raise maxdim until the slope stops moving UP (truncation converged), then extend in
# time and watch it descend toward 4/3. The default parameters here are sized for a
# few-minute demonstration run; treat the printed exponent as an effective crossover value.

using ITensors
using ITensorMPS
using MPSToolkit
using Printf

# ---- parameters (sized for a short demonstration; see header for convergence) -----------
nsites = 48
center = nsites ÷ 2
maxdim = 32
dt     = 0.1            # per-sweep gate time; one constrained call advances t by 2*dt
ncall  = 60             # measurement sweeps -> t_max = 2*dt*ncall = 12.0
cutoff = 1e-12
gate_maxdim = 4 * maxdim
connector_buffer = 8
project_every = 1       # constraint checkpoint after every forward+reverse sweep
fit_window = (5.0, 12.0) # late-window fit for the (crossover) effective exponent
edge_tol = 1e-4          # drop fit times once the spreading front reaches the chain edges
# ------------------------------------------------------------------------------------------

sites = pauli_siteinds(nsites)
terms = [(first(pxp_term_support(nsites, j)), pxp_term_hamiltonian(nsites, j)) for j in 1:nsites]
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

# Initial operator: the PXP term at the center, vectorized and projected into the sector.
# "ProjUp" is the ITensor |0><0| projector ("Up" = our ground state |0>), so the OpSum below
# is exactly P X P centered on `center`.
println("building the sector-projected initial operator O(0) = P_G h_$(center) P_G ...")
phys = siteinds("S=1/2", nsites)
os = OpSum()
os += "ProjUp", center - 1, "X", center, "ProjUp", center + 1
O = pauli_state_from_mpo(MPO(os, phys), sites)
O = apply(projector, O; maxdim=gate_maxdim, cutoff=cutoff)
normalize!(O)   # one overall rescale at t=0 only; the evolution below never renormalizes

correlation_profile(state) = real.(pauli_expectation_profile(state, terms; normalize=false))

times = Float64[]
second_moment = Float64[]
edge_fraction = Float64[]
total_ref = Ref(0.0)        # conserved total at the first measurement
last_report = Ref((0.0, 0.0))

println("evolving with DMT + constraint checkpoints, normalize=false (t_max = $(2 * dt * ncall)) ...")
println(rpad("t", 8), rpad("M2(t)", 12), rpad("local slope", 13), rpad("cons. drift", 13), "front@edges")
for k in 1:ncall
  constrained_dmt_evolve!(O, evo, projector; project_every=project_every, normalize=false)
  p = correlation_profile(O)
  total = sum(p)
  k == 1 && (total_ref[] = total)
  push!(times, 2 * dt * k)
  push!(second_moment, sum((j - center)^2 * p[j] for j in 1:nsites) / total)
  push!(edge_fraction, (abs(p[1]) + abs(p[2]) + abs(p[end - 1]) + abs(p[end])) / maximum(abs, p))
  if k % 5 == 0
    previous_t, previous_m2 = last_report[]
    slope_text = previous_m2 == 0 ? "--" :
      @sprintf("%.2f", log(second_moment[end] / previous_m2) / log(times[end] / previous_t))
    @printf("%-8.2f%-12.4f%-13s%-13.2e%.2e\n", times[end], second_moment[end], slope_text,
      abs(total - total_ref[]) / abs(total_ref[]), edge_fraction[end])
    last_report[] = (times[end], second_moment[end])
  end
end

# Log-log least-squares slope of M2(t) over the fit window, skipping boundary-contaminated
# times. The exponent is 2/z, i.e. 4/3 for KPZ-class z = 3/2.
function fit_m2_slope(times, second_moment, edge_fraction, (tmin, tmax))
  idx = findall(
    i -> tmin <= times[i] <= tmax && second_moment[i] > 0 && edge_fraction[i] < edge_tol,
    eachindex(times),
  )
  lt = log.(times[idx])
  lm = log.(second_moment[idx])
  n = length(lt)
  sx = sum(lt)
  slope = (n * sum(lt .* lm) - sx * sum(lm)) / (n * sum(lt .^ 2) - sx^2)
  return slope, n
end

slope, npts = fit_m2_slope(times, second_moment, edge_fraction, fit_window)
println()
@printf("late-window fit over t in %s (%d points):  effective 2/z = %.3f  =>  z_eff = %.2f\n",
  fit_window, npts, slope, 2 / slope)
println("""
Reference: KPZ-class PXP energy transport has 2/z = 4/3 (z = 3/2), PRX 13, 011033 (2023).
At this maxdim the measured slope mixes two effects: the slow crossover (true slope still
ABOVE 4/3 at these times, descending) and truncation (which suppresses M2 growth and pulls
the slope DOWN — production runs measured ~1.46-1.49 at maxdim=48 vs ~1.32-1.38 at
maxdim=32 over t = 8-14). Raise maxdim until the slope stops moving up, then extend ncall;
the conservation-drift column is the run's truncation error bar.""")
