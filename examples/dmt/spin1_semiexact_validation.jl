# Does DMT actually beat plain SVD at d = 3? A semi-exact, equal-budget comparison.
#
# THIS IS NOT AN EXPONENT BENCHMARK. Its sibling `spin1_melt.jl` is, and that is exactly the
# problem: it extracted z = 1.507 (ULS/SU(3), KPZ target 1.5) and z = 1.514 (spin-1 Heisenberg,
# diffusive target 2.0) -- BOTH models landed on z ~ 1.51, so the run could not tell a
# superdiffusive model from a diffusive one, and the ULS "agreement" with KPZ could be a
# coincidence. A fitted slope over a finite-time window is the weakest evidence available, and
# it is what made the two models look alike.
#
# So test what DMT actually CLAIMS. The claim is not "it recovers the right exponent". It is:
#
#     at a fixed bond dimension, DMT preserves hydrodynamic content that plain SVD truncation
#     destroys.
#
# That is directly falsifiable against an EXACT reference at short times, and it is the
# validation White, Zaletel, Mong and Refael used for the method in the first place
# (PRB 97, 035127 (2018), Sec. IV.B): run the SAME Trotter circuit at a bond dimension large
# enough to be essentially exact, then compare truncation kernels at a small bond dimension where
# truncation is the ONLY difference between the arms.
#
# THREE ARMS, identical in every respect except the truncation kernel -- same initial operator,
# same gates, same bond order, same dt, same measurement, and (within any one comparison) the
# same cutoff:
#   1. REFERENCE  cutoff-only, no bond cap (maxdim = typemax(Int), cutoff = 1e-10).
#   2. DMT        `DMTGateEvolution` / `dmt_evolve!` at a small `maxdim`, cutoff 0.
#   3. SVD        `LocalGateEvolution` / `tebd_evolve!` at the SAME `maxdim`, cutoff 0, plain SVD.
# Arms 2 and 3 truncate by `maxdim` and nothing else (see ARM_CUTOFFS), so the single difference
# between them is WHICH `maxdim` directions the kernel keeps.
# Arms 1 and 3 are the same code path (`LocalGateEvolution`) and differ only in `maxdim`, so the
# reference is literally arm 3 with the cap removed. Arm 2 replays the same gate sequence: one
# `dmt_evolve!` call is a forward sweep over bonds 1..N-1 then a reverse sweep over N-1..1, which
# is precisely the `vcat(1:N-1, reverse(1:N-1))` schedule the other two arms run, so the circuits
# are identical gate for gate (asserted below by `_assert_arms_agree_uncapped`).
#
# PROTOCOL -- the same infinitesimal domain-wall melt as `spin1_melt.jl` and
# `domain_wall_melting.jl`, at N = 12 so an exact reference is reachable:
#       rho(0) = I + sum_j c_j S^z_j,   c_j = -eps (j <= kink), +eps (j > kink),
#       p_x(t) = tr(rho(t) S^z_x)                            [the LOCAL PROFILE, N numbers per t]
# and the metric is the profile error against the reference, site by site, on the eps-normalized
# profile q = p / eps (see EPS_LADDER: the exact melt is EXACTLY linear in eps, so q is one curve),
#       err_inf(t) = max_x |q_x^method(t) - q_x^ref(t)| / max_x |q_x^ref(t)|
#       err_l1(t)  = sum_x |q_x^method(t) - q_x^ref(t)| / sum_x |q_x^ref(t)|.
# `max_x |q^ref|` is the unmelted bulk step height, constant to 6 digits over the whole window, so
# err_inf reads directly as "worst-site error as a fraction of the full step".
#
# THE eps x CUTOFF GRID IS PART OF THE EXPERIMENT, NOT A ROBUSTNESS CHECK. The physics is linear
# in the wall amplitude eps but the truncation is not, and "plain SVD" is really two different
# competitors depending on whether the caller also sets a cutoff. The capped arms are therefore
# run on a 2 x 2 of (arm cutoff 0 or 1e-10) x (eps 0.25 or 0.01), at every rung of the maxdim
# ladder. The cutoff = 0 rows isolate the truncation KERNEL; the cutoff = 1e-10 rows are what a
# caller actually writes -- both sibling melts set exactly that -- and they are where a
# norm-relative criterion can quietly delete an infinitesimal signal. Reporting only one of the
# two would give a misleading verdict in whichever direction the author preferred. The whole grid
# is nearly free because ONE reference scores every cell (see EPS_LADDER).
#
# WHY A PROFILE AND NOT A SCALAR. A scalar summary of a melt can be blind by symmetry. The melt
# profile stays ANTISYMMETRIC about the wall, so for instance
# `sum(p[right half]) / sum(abs, p)` is identically 0.5 at every time, for every method, at every
# bond dimension -- it cannot move, and it will happily print `+0.500000` forever while the
# simulation falls apart. (Measured: it does.) The full profile carries N numbers per time step,
# is not protected by any symmetry, and is exactly the hydrodynamic content whose survival DMT
# claims as its distinguishing feature. It is also the quantity DMT preserves EXACTLY at each
# individual truncation (`preserve_diameter = 3` protects every single-site observable), so any
# error it accumulates arrives through the gates from discarded correlations -- which is the
# mechanism under test.
#
# WHAT THIS TESTS AND WHAT IT DOES NOT. Operator entanglement grows roughly exponentially in a
# melting spin-1 chain, so the uncapped reference is affordable only to t ~ 2 even at N = 12
# (measured chi: 17, 39, 72, 125, 217, 337, 540, 873 at t = 0.2 ... 1.6, and 1.6 s -> 193 s per
# sweep over the same range; the cost roughly triples every 0.2 in t). That is a SHORT window but
# a STRINGENT one: by t = 2 the exact state needs chi ~ 1500, so a maxdim = 20 arm is carrying
# ~1.3% of the bond space and there is ground truth to score it against. This measures the
# TRUNCATION MECHANISM, not the asymptotics. It cannot and does not say anything about z.
#
# THE d = 2 CONTROL is the single most informative addition beyond the core comparison. At d = 2
# operator entanglement grows far more slowly, the exact reference reaches much further in time,
# and DMT is already known to work (the spin-1/2 XXZ melt in `domain_wall_melting.jl` reproduces
# z = 1.02 / 1.62 / 1.79 across the ballistic / KPZ / diffusive regimes). Running the identical
# protocol at d = 2 answers "does d = 3 behave like the case we trust?" -- and, because its
# reference reaches 3x further in time, "does the d = 3 verdict survive past the window where the
# d = 3 reference dies?". The d = 2 ladder is matched in COMPLEMENT budget, not in raw `maxdim`:
# the DMT protected block is 2 d^2 directions of every `maxdim` (18 at d = 3, 8 at d = 2), so
# maxdim = 20 / 28 / 40 / 64 at d = 3 and 10 / 18 / 30 / 54 at d = 2 both buy
# chi' = 2 / 10 / 22 / 46 complement directions. That is the resource DMT actually spends.
#
# BUDGET FLOOR. At d = 3 with the default `preserve_diameter = 3` the protected block is
# 2 * 3^2 = 18 directions and the minimum legal `maxdim` is 19, so `maxdim = 20` (chi' = 2) is
# essentially the tightest budget the kernel admits -- a valid and deliberately extreme rung.
#
# ORACLES (all run before the physics, all cheap, all in this file):
#  * `_assert_step_oracle`  -- at t = 0 the profile must equal the ANALYTIC step
#    p_x = c_x tr((S^z)^2) d^(N-1) exactly (121.5 per site at d = 3, N = 6, eps = 1/4). This pins
#    the literal-operator convention shared by `operator_product_state`,
#    `operator_local_sum_state` and `operator_expectation_profile(...; normalize = false)`: the
#    last returns the LITERAL trace tr(rho O_x), not a normalized ratio.
#  * `_assert_arms_agree_uncapped` -- with NO truncation at all (cutoff 0, maxdim above the exact
#    bond dimension) on a tiny N = 4 chain, all three arms must produce the same profile to
#    machine precision. This is the load-bearing claim of the whole comparison: the arms differ
#    ONLY by truncation. If the schedules or gates disagreed, this fails.
#  * `_assert_eps_linearity` -- the untruncated profile must be exactly proportional to eps, which
#    is what lets one reference score every cell of the eps grid.
#  * `reference_floor` -- the reference is not exact either, it is cutoff-limited. Re-running it at
#    a 100x tighter cutoff over the early window measures its own error, which must sit far below
#    the errors being reported. Printed, not hidden.
#  * `front` -- the melt must not reach the chain edge, or the reference is measuring a reflection.
#    Reported as max_x-in-edge-band |p_x(t)/p_x(0) - 1| on the reference.
#
# WHAT THIS CONFIGURATION ACTUALLY FINDS (N = 12, dt = 0.1, ~55 min wall on 32 BLAS threads;
# numbers below are err_inf at each reference's LAST time, from spin1_semiexact_validation.csv).
#
#   REFERENCE REACH. Heisenberg d = 3 got to t = 1.8 (chi = 1433, 15.0 min) before the wall-clock
#   guard fired; ULS/SU(3) to t = 1.6 (chi = 1465, 10.4 min) -- ULS's operator entanglement grows
#   about twice as fast, so it buys 0.2 less time for the same money. The d = 2 control ran the
#   full t = 6.0 in 1.1 min (chi = 620). The reference's OWN error, measured against a 100x tighter
#   cutoff, is 3.0e-6 at t = 1.2 -- one to four orders of magnitude below every arm error quoted
#   below, except the maxdim = 64 DMT arm at early times, which is within ~10x of it.
#
#   THERE IS NO SINGLE VERDICT: THE ANSWER TURNS ON THE COMPLEMENT BUDGET chi'. At d = 3
#   Heisenberg, t = 1.8, cutoff 0, eps = 0.25 (DMT / SVD / ratio):
#       chi' =  2   4.15e-2 / 6.42e-3 / 0.16     DMT 6x WORSE
#       chi' = 10   3.02e-3 / 2.08e-3 / 0.69     DMT slightly worse
#       chi' = 22   2.81e-3 / 1.14e-3 / 0.41     DMT 2.5x worse
#       chi' = 46   6.33e-5 / 1.89e-4 / 2.98     DMT 3x BETTER, and the gap GROWS with time
#                                                (1.34 at t = 1.4, 2.08 at 1.6, 2.98 at 1.8)
#   The crossover is the structural overhead: the protected block is a FIXED 18 directions, so at
#   maxdim = 20 it is 90% of the budget and DMT is spending almost everything on a subspace plain
#   SVD would not have chosen, losing badly in Frobenius norm (where SVD is optimal by
#   Eckart-Young) and feeding that error back through the gates. At maxdim = 64 the overhead is
#   28%, the complement is resolved, and the protection starts paying. THE PRACTICAL READING: at
#   d = 3, `maxdim` must be several times `2 d^2` before DMT is worth using at all; near the floor
#   it is actively worse than the plain SVD it replaces.
#
#   THE d = 2 CONTROL SAYS d = 3 IS NOT THE PROBLEM. The identical pattern appears at d = 2 with
#   chi' matched: at t = 6.0, cutoff 0, eps = 0.25 the ratios run 1.23 (chi' = 2) / 0.48 (10) /
#   0.22 (22) / 2.14 (46), and at chi' = 46 the DMT advantage peaks at 3.5x around t = 4.2. So the
#   "DMT only wins at generous chi'" behaviour is a property of the METHOD, not of higher spin.
#   Note the d = 2 window runs well past its own front containment (see below), which is why its
#   absolute errors reach 10-20% at the tight rungs -- both methods are broken there, and a ratio
#   between two broken numbers is not worth much.
#
#   WHAT DMT UNAMBIGUOUSLY BUYS: INVARIANCE. DMT's error barely moves across the eps x cutoff grid
#   -- at d = 3, maxdim = 64, t = 1.8 it reads 6.326e-5 / 6.320e-5 / 6.326e-5 / 6.321e-5 over all
#   four cells, a 0.1% spread. Plain SVD's swings by 1.5x at the same rung and by far more earlier.
#   Two concrete consequences:
#     * TRACE. DMT holds tr(rho) to 6e-14..8e-13 in every cell; plain SVD drifts up to 1.3e-4
#       (d = 2, eps = 0.25). The trace is protected structurally, not because it wins a weight
#       contest.
#     * THE CUTOFF FLOOR. At eps = 0.01 with the ordinary cutoff = 1e-10 that both sibling melts
#       set, the SVD arm's error stops depending on maxdim: at t = 1.0 it is 1.02e-4 at chi' = 22
#       and 1.03e-4 at chi' = 46, and it is even NON-MONOTONE in time (2.26e-4 at t = 1.4 ->
#       1.41e-4 at t = 1.8). Buying bond dimension buys nothing, because a norm-relative cutoff
#       cannot see a signal that is 2.8% of the norm. DMT at chi' = 46 is 1.47e-6 there -- 70x
#       better -- and keeps improving with maxdim. THIS is the regime where the method's design
#       matters, and it is exactly the infinitesimal-wall linear-response regime the DMT
#       literature works in.
#
#   FRONT CONTAINMENT (measured, not assumed; edge drift over the outermost 2 sites per side).
#   Heisenberg d = 3 stays under 1e-3 through t = 1.6 (9.0e-4) and reaches 2.3e-3 only at the final
#   t = 1.8. ULS crosses earlier, at t = 1.4 (1.9e-3), reaching 5.4e-3 at t = 1.6. The d = 2
#   control is contained only to t = 2.6 (9.6e-4) and is at 1.4e-1 by t = 6.0 -- its late times are
#   a saturated box, not a melt. This does NOT bias the DMT-vs-SVD comparison, since all three arms
#   share the same finite chain and the reference is exact FOR THAT CHAIN, but the profiles past
#   those times are not bulk hydrodynamics and should not be read as such.
#
#   WHAT THIS DOES NOT SHOW. Nothing about z, and nothing asymptotic: the whole d = 3 window is
#   t <= 1.8. The tempting extrapolation -- "the chi' = 46 advantage grows with time, so it will be
#   large at t = 12" -- is exactly the kind of unfalsifiable claim this file exists to avoid.
#   Testing it needs an exact reference at t ~ 12, which is not reachable at any N.
#
# References:
#   - White, Zaletel, Mong, Refael, PRB 97, 035127 (2018) -- DMT, and the fixed-chi-against-exact
#     validation this script reproduces at d = 3 (their Sec. IV.B is the d = 2 version).
#   - Ye, Machado, Kemp, Hutson, Yao, PRL 129, 230602 (2022), arXiv:2205.02853 -- DMT at d = 3.
#   - Ye, Machado, White, Mong, Yao, arXiv:1902.01859 -- the maxdim = chi_preserve + chi_extra
#     budget convention used here.
#   - examples/dmt/spin1_melt.jl -- the exponent benchmark this script exists to backstop.

using ITensors
using ITensorMPS
using MPSToolkit
using LinearAlgebra: Diagonal, I
using Printf

# ---- parameters -----------------------------------------------------------------------------
# The REFERENCE runs are the entire cost; the capped arms are free by comparison (every capped
# sweep is well under a second at N = 12). Each reference therefore carries a wall-clock budget
# and a projection guard: before each sweep the driver extrapolates the next sweep's cost from the
# growth of the last two and stops if it would blow the budget. The reach is REPORTED, never
# assumed -- `t_reached` is a measured output of this script, not a parameter.
const NSITES        = 12          # 6 sites per side; front containment to t = 2 is verified, not assumed
const DT            = 0.1         # per-sweep gate time; one sweep (forward+reverse) advances t by 2*DT
const EPS_WALL      = 0.25        # reference wall amplitude; see EPS_LADDER for why this one
# The wall amplitude eps is NOT a cosmetic parameter here, and one reference covers all of them.
# The exact evolution is EXACTLY linear in eps: U (I + eps D) U^dagger = I + eps D(t), because
# U I U^dagger = I identically, so the measured profile is p^(eps)_x(t) = eps * P_x(t) with a
# single eps-INDEPENDENT P (tr(S^z_x) = 0 kills the identity's contribution). Every comparison
# below is therefore made on the eps-NORMALIZED profile q = p / eps, and the one reference run at
# EPS_WALL scores arms at every eps. `_assert_eps_linearity` checks this rather than assuming it.
#
# Truncation, however, is NOT linear in eps, and crossing eps with the arm cutoff separates two
# effects that are easy to confuse:
#  * A pure `maxdim` truncation should be nearly eps-INDEPENDENT, and it is worth seeing why. The
#    identity is a PRODUCT operator, so it occupies exactly ONE Schmidt direction at every bond no
#    matter how large it is; the remaining maxdim - 1 directions are ranked within the signal
#    sector eps D(t), which is uniformly scaled by eps, so the ordering -- and hence the kept
#    subspace -- does not move as eps changes. Measuring that rather than assuming it is what
#    makes the eps rows at cutoff = 0 worth printing: they are a null control, and if they are
#    NOT flat something in the comparison is wrong.
#  * A cutoff is a different animal, because it is measured against the FULL Hilbert-Schmidt norm.
#    The signal's share of that norm is 58% at eps = 0.25 but 2.8% at eps = 0.01 on this chain, so
#    the same numerical cutoff is ~20x coarser in signal terms at the smaller wall, and it can
#    discard hydrodynamic content that `maxdim` alone would have kept. DMT is structurally immune:
#    the protected block is carried through by construction, not because it won a weight contest.
# eps = 0.25 is the sibling `spin1_melt.jl` convention and the value at which the reference itself
# is most accurate (its own cutoff error is also norm-relative); eps = 0.01 is the genuinely
# infinitesimal, linear-response wall the DMT literature melts.
const EPS_LADDER    = (0.25, 0.01)
const CUTOFF        = 1e-10       # the REFERENCE's only truncation
# The capped arms are run at TWO cutoffs, because these are two different experiments and they do
# not have to give the same answer:
#   cutoff = 0     `maxdim` is the ONLY truncation, so the single difference between arms 2 and 3
#                  is WHICH maxdim directions the kernel keeps. This isolates the DMT rule itself
#                  and is the setting most favourable to plain SVD, whose choice is optimal in
#                  Frobenius norm by Eckart-Young.
#   cutoff = 1e-10 what a caller actually writes, and what both sibling melts use. ITensor's
#                  cutoff bounds DISCARDED WEIGHT RELATIVE TO THE FULL Hilbert-Schmidt norm --
#                  which the identity dominates once eps is small -- so it can discard signal that
#                  `maxdim` would have kept. DMT is structurally immune: its protected block
#                  survives the final refactorization by construction. Whether that immunity is
#                  worth anything in practice is a measurement, and it is made here.
const ARM_CUTOFFS   = (0.0, 1e-10)
const FLOOR_CUTOFF  = 1e-12       # tighter cutoff used to measure the reference's own error
const FLOOR_NCALL   = 6           # sweeps of the reference-floor probe (t <= 1.2, ~3 min)
const EDGE_BAND     = 2           # outermost sites per side used by the front-containment guard
const FRONT_TOL     = 1e-3        # a reference whose front exceeds this is measuring its own edge
const TINY_NSITES   = 6           # chain for the analytic t = 0 step oracle (no evolution, free)
const AGREE_NSITES  = 4           # chain for the uncapped three-arm agreement oracle: with NO
const AGREE_NCALL   = 3           # truncation the bond saturates at d^(2 * N/2), so keep N small
const CSV_PATH      = joinpath(@__DIR__, "spin1_semiexact_validation.csv")

# (label, model, d, t_max, maxdim ladder, reference wall-clock budget in seconds).
# Order matters: the core d = 3 Heisenberg comparison runs first, the d = 2 control second, and
# the more expensive ULS/SU(3) case last, so a truncated run still ships the core result and its
# control. Each case writes its own rows to the CSV as soon as it finishes.
# The two ladders are matched in COMPLEMENT budget chi' = maxdim - 2 d^2 = 2 / 10 / 22 / 46, not
# in raw maxdim -- chi' is the resource DMT spends, and the protected block is 18 at d = 3 vs 8 at
# d = 2. The ladder deliberately spans the point where the structural overhead stops dominating:
# at maxdim = 20 the protected block is 90% of the d = 3 budget, at maxdim = 64 only 28%, so if
# DMT's ranking against plain SVD depends on that fraction the ladder will show it.
const CASES = (
  (label = "heisenberg_d3", model = "heisenberg", d = 3, tmax = 2.0,
   maxdims = (20, 28, 40, 64), ref_budget = 1100.0),
  (label = "heisenberg_d2", model = "heisenberg", d = 2, tmax = 6.0,
   maxdims = (10, 18, 30, 54), ref_budget = 900.0),
  (label = "uls_d3", model = "uls", d = 3, tmax = 2.0,
   maxdims = (20, 28, 40, 64), ref_budget = 1100.0),
)
# -------------------------------------------------------------------------------------------------

"""
    spin_matrices(d)

Return `(Sx, Sy, Sz)` for spin `(d - 1) / 2` in the descending-`m` basis, from the ladder
coefficients `S^+ |m> = sqrt(S(S+1) - m(m+1)) |m+1>`.

`d = 2` gives the spin-1/2 matrices `sigma^a / 2` and `d = 3` the spin-1 generators (so
`S^z = diag(1, 0, -1)`). Written out here because MPSToolkit deliberately ships no spin-`S` model
builders; any source of dense `d x d` matrices using the same `kron`-left-to-right multi-site
convention works instead. EDKit.jl is NOT a dependency of MPSToolkit.
"""
function spin_matrices(d::Integer)
  d >= 2 || throw(ArgumentError("spin_matrices needs d >= 2"))
  s = (d - 1) / 2
  ms = [s - (k - 1) for k in 1:d]                      # m = s, s-1, ..., -s down the basis
  raise = zeros(ComplexF64, d, d)
  for k in 2:d                                          # S^+ |m> ~ |m+1>: column k -> row k-1
    m = ms[k]
    raise[k - 1, k] = sqrt(s * (s + 1) - m * (m + 1))
  end
  lower = collect(raise')
  return (raise + lower) / 2, (raise - lower) / (2im), Matrix{ComplexF64}(Diagonal(ms))
end

# Two-site bond Hamiltonians. `kron` is applied left-to-right (first factor = leftmost site), the
# convention `operator_gate` and every operator-space state builder in MPSToolkit use.
function heisenberg_bond(d::Integer)
  sx, sy, sz = spin_matrices(d)
  return kron(sx, sx) + kron(sy, sy) + kron(sz, sz)
end

# ULS / SU(3) at d = 3: H = S.S + (S.S)^2, the integrable Uimin-Lai-Sutherland point.
function uls_bond(d::Integer)
  d == 3 || throw(ArgumentError("the ULS point is defined here only at d = 3"))
  h = heisenberg_bond(d)
  return h + h * h
end

bond_hamiltonian(model, d) =
  model == "heisenberg" ? heisenberg_bond(d) :
  model == "uls" ? uls_bond(d) :
  throw(ArgumentError("unknown model $(model)"))

"""
    melt_state(d, nsites; eps_wall)

Build the literal step density operator `rho(0) = I + sum_j c_j S^z_j` with `c_j = -eps` left of
the central bond and `+eps` right of it, together with the `S^z` matrix used to read it back.

Both builders share the LITERAL-operator convention (an identity argument contributes amplitude
`sqrt(d)`), so `add`ing them gives exactly `I + sum_j c_j S^z_j`, trace `d^N` to machine
precision -- a traceful density operator, which is the object DMT is built for.
"""
function melt_state(d::Integer, nsites::Integer; eps_wall = EPS_WALL)
  _, _, sz = spin_matrices(d)
  sites = operator_siteinds(nsites; d = d)
  identity_matrix = Matrix{ComplexF64}(I, d, d)
  kink = nsites ÷ 2
  coeffs = [j <= kink ? -eps_wall : eps_wall for j in 1:nsites]
  rho = add(operator_product_state(sites, fill(identity_matrix, nsites)),
            operator_local_sum_state(sites, sz, coeffs); maxdim = 8, cutoff = 0.0)
  return rho, sz
end

# The local magnetization profile p_x = tr(rho S^z_x). `normalize = false` returns the LITERAL
# trace (it carries the d^((N - span) / 2) basis prefactor internally), so p_x(0) is the analytic
# c_x tr((S^z)^2) d^(N-1) exactly -- see `_assert_step_oracle`. Every error below is a RATIO
# against the reference's own maximum, hence blind to that prefactor.
profile(rho::MPS, sz) =
  real.(operator_expectation_profile(rho, [(x, sz) for x in 1:length(rho)]; normalize = false))

# The three arms. `svd_evolution` doubles as the reference: `maxdim = typemax(Int)` removes the
# cap and leaves the cutoff as the only truncation, and `maxdim = 20` makes it the plain-SVD arm.
# Nothing else differs.
function svd_evolution(gate, dt, nsites; maxdim, cutoff = CUTOFF)
  forward = collect(1:(nsites - 1))
  return LocalGateEvolution(gate, 2 * dt; schedule = vcat(forward, reverse(forward)),
    nstep = 1, maxdim = maxdim, cutoff = cutoff)
end

function dmt_evolution(gate, dt, nsites; maxdim, cutoff = CUTOFF)
  forward = collect(1:(nsites - 1))
  return DMTGateEvolution(gate, dt; schedule = forward, reverse_schedule = reverse(forward),
    nstep = 1, maxdim = maxdim, cutoff = cutoff,
    normalize = false)  # absolute profiles -> no per-sweep rescaling
end

"""
    _assert_step_oracle()

At `t = 0` the measured profile must equal the analytic step `c_x tr((S^z)^2) d^(N-1)` exactly.

This is a convention oracle with no physics in it: it fails if `operator_product_state`,
`operator_local_sum_state` or the unnormalized `operator_expectation_profile` prefactor ever stop
agreeing on what "the literal operator `I + sum_j c_j S^z_j`" means.
"""
function _assert_step_oracle()
  for d in (2, 3)
    rho, sz = melt_state(d, TINY_NSITES)
    p = profile(rho, sz)
    kink = TINY_NSITES ÷ 2
    trace_sz2 = real(sum(abs2, sz))                                   # tr((S^z)^2)
    want = [(j <= kink ? -EPS_WALL : EPS_WALL) * trace_sz2 * Float64(d)^(TINY_NSITES - 1)
        for j in 1:TINY_NSITES]
    @assert maximum(abs, p .- want) <= 1e-9 * maximum(abs, want) "t = 0 step oracle failed at d = $d"
  end
  return nothing
end

"""
    _assert_arms_agree_uncapped()

With truncation switched off entirely, the three arms must be the SAME circuit.

`maxdim` above the exact bond dimension (`d^(2 * (N/2))` at the center) and `cutoff = 0` means no
arm discards anything, so any disagreement is a difference in the gate sequence or the bond order,
not in the truncation kernel -- and the whole benchmark is meaningless if the arms differ that
way. Run on `TINY_NSITES` sites, where the exact bond dimension is small enough to carry.
"""
function _assert_arms_agree_uncapped(; d = 3, nsites = AGREE_NSITES, ncall = AGREE_NCALL)
  exact_dim = d^(2 * (nsites ÷ 2))
  gate = operator_gate_from_hamiltonian(bond_hamiltonian("heisenberg", d), DT; d = d)
  rho0, sz = melt_state(d, nsites)
  states = (deepcopy(rho0), deepcopy(rho0), deepcopy(rho0))
  evos = (svd_evolution(gate, DT, nsites; maxdim = typemax(Int), cutoff = 0.0),
      svd_evolution(gate, DT, nsites; maxdim = exact_dim, cutoff = 0.0),
      dmt_evolution(gate, DT, nsites; maxdim = exact_dim, cutoff = 0.0))
  for _ in 1:ncall
    evolve!(states[1], evos[1])
    evolve!(states[2], evos[2])
    dmt_evolve!(states[3], evos[3])
  end
  reference = profile(states[1], sz)
  scale = maximum(abs, reference)
  worst = maximum(maximum(abs, profile(s, sz) .- reference) / scale for s in states[2:3])
  @assert worst <= 1e-10 "uncapped arms disagree by $(worst); they are not the same circuit"
  return worst
end

"""
    _assert_eps_linearity()

The untruncated profile must be EXACTLY proportional to the wall amplitude `eps`.

This is what lets one reference run score arms at every `eps` in `EPS_LADDER`:
`U (I + eps D) U^dagger = I + eps D(t)` because `U I U^dagger = I`, and `tr(S^z_x) = 0`, so
`p^(eps)(t) = eps P(t)` with an `eps`-independent `P`. Checked with truncation switched off, since
truncation is precisely the thing that breaks the proportionality -- which is the effect under
study, not a bug.
"""
function _assert_eps_linearity(; d = 3, nsites = AGREE_NSITES, ncall = AGREE_NCALL)
  exact_dim = d^(2 * (nsites ÷ 2))
  gate = operator_gate_from_hamiltonian(bond_hamiltonian("heisenberg", d), DT; d = d)
  normalized = map(EPS_LADDER) do eps
    rho, sz = melt_state(d, nsites; eps_wall = eps)
    evo = svd_evolution(gate, DT, nsites; maxdim = exact_dim, cutoff = 0.0)
    for _ in 1:ncall
      evolve!(rho, evo)
    end
    return profile(rho, sz) ./ eps
  end
  scale = maximum(abs, first(normalized))
  worst = maximum(maximum(abs, q .- first(normalized)) / scale for q in normalized)
  @assert worst <= 1e-10 "eps-normalized profiles disagree by $(worst); the melt is not linear in eps"
  return worst
end

"""
    run_arm(kind, model, d, nsites, dt, ncall; maxdim, cutoff, budget, verbose)

Evolve one arm and return its per-sweep profiles.

`kind` is `:svd` (plain SVD truncation via `LocalGateEvolution` / `tebd_evolve!`; also the
reference when `maxdim = typemax(Int)`) or `:dmt` (`DMTGateEvolution` / `dmt_evolve!`).

`budget` is a wall-clock guard for the uncapped reference. Before each sweep the driver
extrapolates the next sweep's cost from the growth ratio of the last two measured sweeps and
stops early rather than overrunning; the achieved `ncall` is returned so the capped arms can be
scored on exactly the times the reference reached.
"""
function run_arm(kind::Symbol, model, d, nsites, dt, ncall;
                 maxdim, eps_wall = EPS_WALL, cutoff = CUTOFF, budget = Inf, verbose = false)
  gate = operator_gate_from_hamiltonian(bond_hamiltonian(model, d), dt; d = d)
  rho, sz = melt_state(d, nsites; eps_wall = eps_wall)
  evo = kind === :dmt ? dmt_evolution(gate, dt, nsites; maxdim = maxdim, cutoff = cutoff) :
                        svd_evolution(gate, dt, nsites; maxdim = maxdim, cutoff = cutoff)
  times = [0.0]
  profiles = [profile(rho, sz)]
  chis = [maxlinkdim(rho)]
  traces = [real(operator_trace(rho))]
  elapsed = 0.0
  per_sweep = Float64[]
  stopped = "completed"
  for k in 1:ncall
    if length(per_sweep) >= 2
      projected = per_sweep[end] * max(1.0, per_sweep[end] / per_sweep[end - 1])
      if elapsed + projected > budget
        stopped = @sprintf("stopped on budget: next sweep projected %.0f s, %.0f s of %.0f s used",
                  projected, elapsed, budget)
        break
      end
    end
    wall = time()
    kind === :dmt ? dmt_evolve!(rho, evo) : evolve!(rho, evo)
    push!(per_sweep, time() - wall)
    elapsed += per_sweep[end]
    push!(times, 2 * dt * k)
    push!(profiles, profile(rho, sz))
    push!(chis, maxlinkdim(rho))
    push!(traces, real(operator_trace(rho)))
    if verbose
      @printf("      %-4s maxdim=%-8s t=%4.1f  chi=%6d  %8.2f s  (cum %7.1f s)\n",
          kind, maxdim == typemax(Int) ? "inf" : string(maxdim), times[end], chis[end],
          per_sweep[end], elapsed)
      flush(stdout)
    end
  end
  return (; times, profiles, chis, traces, seconds = elapsed, stopped)
end

# Profile errors against the reference, both normalized by the reference's OWN magnitude so they
# are blind to the d^((N - span) / 2) basis prefactor and to any overall norm drift in the
# reference. `err_inf` is the worst site as a fraction of the unmelted step height.
relative_inf_error(p, reference) = maximum(abs, p .- reference) / maximum(abs, reference)
relative_l1_error(p, reference) = sum(abs, p .- reference) / sum(abs, reference)

"""
    front_curve(profiles, band)

Fractional change of the outermost `band` sites per side, one entry per recorded time.

The melt must stay away from the chain edge: once it arrives, the reference is measuring a
reflection off its own boundary. That does NOT invalidate the DMT-vs-SVD comparison -- every arm
runs on the same finite chain and is scored against the reference for THAT chain -- but it does
mean the profile past the containment time is no longer a bulk hydrodynamic profile, so the whole
curve is printed and the crossing time reported rather than a single pass/fail.
"""
function front_curve(profiles, band::Integer)
  p0 = first(profiles)
  n = length(p0)
  edge = vcat(1:band, (n - band + 1):n)
  return [maximum(abs(p[x] / p0[x] - 1) for x in edge) for p in profiles]
end

# Last time at which the edge drift is still below `tol`; `nothing` if it never crosses.
function front_contained_until(times, fronts, tol)
  breach = findfirst(>=(tol), fronts)
  isnothing(breach) && return last(times)
  return breach == 1 ? times[1] : times[breach - 1]
end

"""
    reference_floor(model, d, nsites, dt, ncall)

Measure the reference's own truncation error by re-running it at a 100x tighter cutoff.

The "exact" arm is cutoff-limited, not exact. ITensor's `cutoff` bounds the DISCARDED WEIGHT, so a
`1e-10` cutoff costs ~`1e-5` in amplitude per truncation and the errors accumulate over the
`2 (N - 1)` gates of every sweep. This returns the largest profile disagreement between the
shipped reference and a tighter one over the early window; it must sit far below the DMT and SVD
errors being reported, and it is printed so a reader can check that rather than take it on trust.
"""
function reference_floor(model, d, nsites, dt, ncall)
  loose = run_arm(:svd, model, d, nsites, dt, ncall; maxdim = typemax(Int), cutoff = CUTOFF)
  tight = run_arm(:svd, model, d, nsites, dt, ncall; maxdim = typemax(Int), cutoff = FLOOR_CUTOFF)
  n = min(length(loose.profiles), length(tight.profiles))
  worst = maximum(relative_inf_error(loose.profiles[k], tight.profiles[k]) for k in 1:n)
  return (; worst, t = 2 * dt * (n - 1), chi_loose = loose.chis[n], chi_tight = tight.chis[n],
      seconds = loose.seconds + tight.seconds)
end

# One case: the reference, then the DMT and SVD arms at every rung of the ladder, scored on the
# times the reference actually reached.
function run_case(case; nsites = NSITES, dt = DT, verbose = true)
  ncall = round(Int, case.tmax / (2 * dt))
  protected = 2 * case.d^2
  @printf("\n[%s]  d=%d  model=%s  N=%d  dt=%.2f  t_max=%.1f  protected block 2*%d^2 = %d\n",
      case.label, case.d, case.model, nsites, dt, case.tmax, case.d, protected)
  println("  reference: cutoff-only, maxdim = typemax(Int), cutoff = $(CUTOFF)")
  flush(stdout)

  reference = run_arm(:svd, case.model, case.d, nsites, dt, ncall;
            maxdim = typemax(Int), budget = case.ref_budget, verbose = verbose)
  reached = length(reference.times) - 1
  fronts = front_curve(reference.profiles, EDGE_BAND)
  front = maximum(fronts)
  contained = front_contained_until(reference.times, fronts, FRONT_TOL)
  @printf("  reference reached t = %.1f (%d/%d sweeps, chi = %d, %.1f min) -- %s\n",
      reference.times[end], reached, ncall, reference.chis[end], reference.seconds / 60,
      reference.stopped)
  @printf("  front containment: contained (edge drift < %.0e) through t = %.1f; max %.2e at t = %.1f\n",
      FRONT_TOL, contained, front, reference.times[end])
  println("  front curve: ",
      join((@sprintf("%.1f:%.1e", reference.times[k], fronts[k]) for k in 2:(reached + 1)), "  "))
  flush(stdout)

  # The reference is scored in eps-NORMALIZED form q = p / eps, which by `_assert_eps_linearity`
  # is the same curve for every wall amplitude -- so this single (expensive) run is the ground
  # truth for every rung of the eps ladder, and the eps sweep costs only capped arms.
  qref = [p ./ EPS_WALL for p in reference.profiles]

  arms = Dict{Tuple{Float64,Float64,Symbol,Int},NamedTuple}()
  arm_seconds = 0.0
  for cut in ARM_CUTOFFS, eps in EPS_LADDER, maxdim in case.maxdims, kind in (:dmt, :svd)
    arm = run_arm(kind, case.model, case.d, nsites, dt, reached;
                  maxdim = maxdim, eps_wall = eps, cutoff = cut)
    q = [p ./ eps for p in arm.profiles]
    errs_inf = [relative_inf_error(q[k], qref[k]) for k in 1:(reached + 1)]
    errs_l1 = [relative_l1_error(q[k], qref[k]) for k in 1:(reached + 1)]
    arms[(cut, eps, kind, maxdim)] = (; arm..., errs_inf, errs_l1)
    arm_seconds += arm.seconds
  end
  @printf("  capped arms: %d runs (%d cutoffs x %d eps x %d rungs x 2 kinds) in %.1f min total\n",
      2 * length(case.maxdims) * length(EPS_LADDER) * length(ARM_CUTOFFS),
      length(ARM_CUTOFFS), length(EPS_LADDER), length(case.maxdims), arm_seconds / 60)

  for cut in ARM_CUTOFFS, eps in EPS_LADDER
    @printf("\n  arm cutoff = %.0e   eps = %.3f\n", cut, eps)
    println("  ", rpad("maxdim", 8), rpad("chi'", 6), rpad("t", 6), rpad("err_inf DMT", 14),
        rpad("err_inf SVD", 14), rpad("SVD/DMT", 10), rpad("err_l1 DMT", 13),
        rpad("err_l1 SVD", 13), "chi_ref")
    for maxdim in case.maxdims
      dmt, svd = arms[(cut, eps, :dmt, maxdim)], arms[(cut, eps, :svd, maxdim)]
      for k in 2:(reached + 1)
        ratio = dmt.errs_inf[k] > 0 ? svd.errs_inf[k] / dmt.errs_inf[k] : NaN
        @printf("  %-8d%-6d%-6.1f%-14.3e%-14.3e%-10.3g%-13.3e%-13.3e%d\n",
            maxdim, maxdim - protected, reference.times[k], dmt.errs_inf[k],
            svd.errs_inf[k], ratio, dmt.errs_l1[k], svd.errs_l1[k], reference.chis[k])
      end
    end
  end
  flush(stdout)
  return (; case, reference, arms, reached, front, protected)
end

function append_csv(io, result)
  case, reference, arms = result.case, result.reference, result.arms
  for k in 1:(result.reached + 1)
    @printf(io, "%s,%d,%s,reference,%.1e,%.4f,0,0,%.4f,%d,%.6e,%.6e,%.6e,%.6e\n",
        case.label, case.d, case.model, CUTOFF, EPS_WALL, reference.times[k],
        reference.chis[k], 0.0, 0.0,
        abs(reference.traces[k] / reference.traces[1] - 1), result.front)
  end
  for cut in ARM_CUTOFFS, eps in EPS_LADDER, maxdim in case.maxdims, kind in (:dmt, :svd)
    arm = arms[(cut, eps, kind, maxdim)]
    for k in 1:(result.reached + 1)
      @printf(io, "%s,%d,%s,%s,%.1e,%.4f,%d,%d,%.4f,%d,%.6e,%.6e,%.6e,%.6e\n",
          case.label, case.d, case.model, kind, cut, eps, maxdim,
          maxdim - result.protected, reference.times[k], arm.chis[k],
          arm.errs_inf[k], arm.errs_l1[k],
          abs(arm.traces[k] / arm.traces[1] - 1), result.front)
    end
  end
  flush(io)
  return io
end

function main(; cases = CASES, nsites = NSITES, dt = DT, csv_path = CSV_PATH, verbose = true,
              floor_probe = true)
  println("Semi-exact validation of the DMT truncation kernel at d = 3")
  println("  three arms, identical circuits, differing ONLY in truncation:")
  println("    reference  LocalGateEvolution, maxdim = typemax(Int), cutoff = $(CUTOFF)")
  println("    dmt        DMTGateEvolution   at a small maxdim (preserve_diameter = 3)")
  println("    svd        LocalGateEvolution at the SAME maxdim, plain SVD truncation")
  println("  metric: max_x |q_method - q_ref| / max_x |q_ref| on the eps-normalized S^z profile")
  @printf("  capped arms run on cutoff = %s x eps = %s x the maxdim ladder\n",
      string(ARM_CUTOFFS), string(EPS_LADDER))

  _assert_step_oracle()
  worst = _assert_arms_agree_uncapped()
  linearity = _assert_eps_linearity()
  @printf("\noracles: t=0 analytic step OK (d = 2 and 3); uncapped three-arm agreement %.2e; ",
      worst)
  @printf("eps-linearity %.2e over eps = %s\n", linearity, string(EPS_LADDER))
  flush(stdout)

  if floor_probe
    fl = reference_floor("heisenberg", 3, nsites, dt, FLOOR_NCALL)
    @printf("reference floor: cutoff %.0e vs %.0e agree to %.2e up to t = %.1f (chi %d vs %d, %.1f min)\n",
        CUTOFF, FLOOR_CUTOFF, fl.worst, fl.t, fl.chi_loose, fl.chi_tight, fl.seconds / 60)
    flush(stdout)
  end

  results = NamedTuple[]
  open(csv_path, "w") do io
    println(io, "case,d,model,arm,cutoff,eps,maxdim,chi_prime,t,chi,err_inf,err_l1,trace_rel_err,front_ref")
    for case in cases
      result = run_case(case; nsites = nsites, dt = dt, verbose = verbose)
      append_csv(io, result)
      push!(results, result)
    end
  end

  println("\n", "="^108)
  println("VERDICT TABLE -- profile error at the LAST time each reference reached")
  println(rpad("case", 16), rpad("t_ref", 7), rpad("cutoff", 9), rpad("eps", 7), rpad("maxdim", 8),
      rpad("chi'", 6), rpad("err_inf DMT", 14), rpad("err_inf SVD", 14), rpad("SVD/DMT", 10),
      "front")
  for result in results
    for cut in ARM_CUTOFFS, eps in EPS_LADDER, maxdim in result.case.maxdims
      dmt = result.arms[(cut, eps, :dmt, maxdim)]
      svd = result.arms[(cut, eps, :svd, maxdim)]
      de, se = dmt.errs_inf[end], svd.errs_inf[end]
      @printf("%-16s%-7.1f%-9.0e%-7.3f%-8d%-6d%-14.3e%-14.3e%-10.3g%.1e\n",
          result.case.label, result.reference.times[end], cut, eps, maxdim,
          maxdim - result.protected, de, se, de > 0 ? se / de : NaN, result.front)
    end
  end
  println("\nwrote $csv_path")
  println("READ THIS BEFORE THE TABLE. 'SVD/DMT' is the ratio of plain-SVD to DMT profile error")
  println("at the SAME bond dimension, same circuit, same times: > 1 means DMT is buying")
  println("something, < 1 means plain SVD is strictly better at that budget. It is NOT constant")
  println("down the maxdim ladder: the DMT protected block is a fixed 2 d^2 directions, so a")
  println("budget close to that floor is spent almost entirely on structure and DMT loses, while")
  println("a budget several times larger resolves the complement and DMT wins. Read the cutoff")
  println("rows too: at cutoff = 0 the eps rows are a null control and should be flat (they are),")
  println("whereas at a nonzero cutoff -- what a caller actually writes -- the SVD arm can hit an")
  println("error floor that no extra maxdim removes, because the cutoff is relative to the FULL")
  println("norm and the identity dominates it at small eps. DMT is structurally immune to that.")
  println("The window is SHORT by construction (the exact reference dies around t = 2 at d = 3):")
  println("this measures the truncation MECHANISM at severe budgets, and says nothing about")
  println("dynamical exponents or asymptotics.")
  return results
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
