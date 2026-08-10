# Does DMT beat plain SVD at `d = 3`? A semi-exact validation

**Answer: near the `d = 3` budget floor, no — and that negative result is what this benchmark
establishes.** At `maxdim = 20` (the tightest budget the kernel admits at `d = 3`) DMT's profile
error is **at least 38.7x worse** than plain SVD's at the same bond dimension, after subtracting
the reference's own error; a `chi'`-matched `d = 2` control reproduces the loss, so it is a
property of the method and not of higher spin. No budget-only rule predicts where this reverses:
at the same `d`, the same complement budget and the same time, spin-1 Heisenberg is a resolved
loss while ULS is a resolved win. The narrower positive claim that survives is that DMT is immune
to a norm-relative cutoff which pins plain SVD at an error floor no extra bond dimension removes.

> *The numbered sections below were written before the reference's own error was measured, and
> overstate the positive half. See both Addenda; the figures immediately above are the corrected
> ones.*

> **Read the Addendum first.** Sections 1-6 are the original write-up. A review found the
> negative result sound but *under*-sold, and the positive half (the `chi' = 46` "win" and the
> "0.1% invariance") unsupported — both were largely measurements of the reference's own error.
> The Addendum at the end supersedes those specific claims, with per-case floors now measured and
> subtracted. The negative result is unchanged except that it is stronger.

Everything below is measured on branch `feat/higher-spin-dmt`, script
`examples/dmt/spin1_semiexact_validation.jl`, data `examples/dmt/spin1_semiexact_validation.csv`
(untracked, matching its siblings). Total compute ~55 min on 32 BLAS threads.

---

## 1. What was run

The premise: `spin1_melt.jl` extracted `z = 1.507` (ULS/SU(3), KPZ target 1.5) and `z = 1.514`
(spin-1 Heisenberg, diffusive target 2.0) — **both ≈ 1.51**, so it could not discriminate a
superdiffusive model from a diffusive one and the ULS "agreement" could be coincidence. Rather than
chase a longer exponent run, this benchmark tests DMT's actual claim: *at fixed bond dimension it
preserves hydrodynamic content that plain SVD destroys*. That is falsifiable against an exact
reference at accessible times — the validation of White/Zaletel/Mong/Refael (PRB 97 035127,
Sec. IV.B), lifted to `d = 3`.

**Three arms, identical circuits, differing only in truncation.** Same initial operator, same
gates, same bond order, same `dt`, same measurement:

| arm | driver | truncation |
|---|---|---|
| reference | `LocalGateEvolution` / `evolve!` | `maxdim = typemax(Int)`, `cutoff = 1e-10` |
| dmt | `DMTGateEvolution` / `dmt_evolve!` | `maxdim` from the ladder, `preserve_diameter = 3` |
| svd | `LocalGateEvolution` / `tebd_evolve!` | the *same* `maxdim`, plain SVD |

One `dmt_evolve!` call is a forward sweep over bonds `1..N-1` then a reverse sweep over `N-1..1`;
the other two arms run exactly that as `vcat(1:N-1, reverse(1:N-1))`, so the circuits are identical
gate for gate. The reference is literally the SVD arm with the cap removed.

**Protocol.** `rho(0) = I + sum_j c_j S^z_j`, `c_j = ∓eps` about the central bond, `N = 12`,
`dt = 0.1`, spin matrices built inline from ladder coefficients (no EDKit dependency). Measured
quantity is the full local profile `p_x(t) = tr(rho(t) S^z_x)`, all `N` sites at every time.

**Metric.** `err_inf(t) = max_x |q_x^method − q_x^ref| / max_x |q_x^ref|` on the eps-normalized
profile `q = p / eps`; `err_l1` likewise with sums. `max_x |q^ref|` is the unmelted step height
(constant to 6 digits), so `err_inf` reads as "worst-site error as a fraction of the full step".
This is scale-free, so the `d^((N-span)/2)` prefactor of
`operator_expectation_profile(...; normalize = false)` cancels.

*On the scalar trap you flagged:* `sum(p[right half]) / sum(abs, p)` is identically `0.5` by the
antisymmetry of the melt and cannot move. The profile comparison has no such blind spot — it
carries `N` numbers per time and is exactly the quantity DMT claims to protect.

**Grid.** Four `maxdim` rungs × two arm cutoffs × two wall amplitudes, per case:

| | `d = 3` | `d = 2` | `chi' = maxdim − 2d²` |
|---|---|---|---|
| ladder | 20 / 28 / 40 / 64 | 10 / 18 / 30 / 54 | 2 / 10 / 22 / 46 |

The ladders are matched in **complement budget** `chi'`, not raw `maxdim`, because `chi'` is the
resource DMT actually spends (protected block 18 at `d = 3`, 8 at `d = 2`). Arm cutoffs `0` and
`1e-10`; wall amplitudes `eps = 0.25` (the sibling's convention) and `eps = 0.01` (linear response).

**One reference scores every cell.** The exact evolution is exactly linear in `eps`
(`U I U† = I`, and `tr(S^z_x) = 0`), so `p^(eps)(t) = eps · P(t)` with an `eps`-independent `P`.
Asserted in-file to `1.2e-12`, not assumed. So the `eps` sweep costs only cheap capped arms.

**Oracles (all in-file, all passed).**

| oracle | result |
|---|---|
| `t = 0` profile equals analytic `c_x tr((S^z)²) d^(N−1)`, `d = 2` and `3` | exact |
| three arms with truncation off produce the same profile | `0.00e+00` |
| eps-linearity of the untruncated profile | `1.18e-12` |
| reference's own error (cutoff `1e-10` vs `1e-12`, to `t = 1.2`) | **`3.03e-06`** |

That last one is the honest bound on the whole experiment: the "exact" arm is cutoff-limited, and
its error sits 1–4 orders of magnitude below every arm error quoted below.

## 2. Reference reach — where it had to stop

Operator entanglement grows roughly exponentially, so each reference carries a wall-clock budget
with a projection guard; the reach is a measured output, not a parameter.

| case | reached | `chi` there | cost | why it stopped |
|---|---|---|---|---|
| Heisenberg `d = 3` | **t = 1.8** | 1433 | 15.0 min | next sweep projected 1990 s, budget 1100 s |
| ULS/SU(3) `d = 3` | **t = 1.6** | 1465 | 10.4 min | next sweep projected 912 s |
| Heisenberg `d = 2` | **t = 6.0** | 620 | 1.1 min | completed |

Reference `chi(t)` at `d = 3` Heisenberg: 17, 39, 72, 125, 217, 337, 540, 873, **1433** at
`t = 0.2 … 1.8` — sweep cost roughly triples every `0.2` in `t`. ULS grows about twice as fast
(`chi = 585` at `t = 1.0` against Heisenberg's 217), which is why it buys 0.2 less time for the
same money. At `t = 1.8` a `maxdim = 20` arm is carrying **1.4%** of the bond space with ground
truth to score it against — a short window, but a stringent one.

## 3. Front containment (measured, not assumed)

Max fractional drift of the outermost 2 sites per side, threshold `1e-3`:

| case | contained through | value at last time |
|---|---|---|
| Heisenberg `d = 3` | **t = 1.6** (9.0e-4) | 2.3e-3 at t = 1.8 |
| ULS `d = 3` | **t = 1.2** (5.5e-4) | 5.4e-3 at t = 1.6 |
| Heisenberg `d = 2` | **t = 2.6** (9.6e-4) | 1.4e-1 at t = 6.0 |

The `d = 3` breaches are mild and only at the final point(s). The `d = 2` control runs well past
containment — by `t = 6` it is a saturated box, not a melt. **This does not bias the DMT-vs-SVD
comparison**: all three arms share the same finite chain and the reference is exact *for that
chain*. It does mean profiles past those times are not bulk hydrodynamics.

## 4. The error tables

`err_inf` at each reference's **last** time. Ratio > 1 means DMT wins.

### Heisenberg `d = 3`, t = 1.8 (χ_ref = 1433)

| chi' (maxdim) | cut 0, eps .25 | cut 0, eps .01 | cut 1e-10, eps .25 | cut 1e-10, eps .01 |
|---|---|---|---|---|
| **2** (20) DMT | 4.15e-2 | 1.11e-2 | 4.09e-2 | 1.08e-2 |
| **2** (20) SVD | 6.42e-3 | 6.63e-3 | 6.84e-3 | 6.50e-3 |
| **2** SVD/DMT | **0.16** | 0.60 | 0.17 | 0.60 |
| **10** (28) DMT | 3.02e-3 | 3.87e-3 | 3.02e-3 | 3.87e-3 |
| **10** (28) SVD | 2.08e-3 | 2.21e-3 | 2.19e-3 | 2.08e-3 |
| **10** SVD/DMT | 0.69 | 0.57 | 0.72 | 0.54 |
| **22** (40) DMT | 2.81e-3 | 2.89e-3 | 2.81e-3 | 2.89e-3 |
| **22** (40) SVD | 1.14e-3 | 1.14e-3 | 1.14e-3 | 9.23e-4 |
| **22** SVD/DMT | 0.41 | 0.39 | 0.41 | 0.32 |
| **46** (64) DMT | 6.33e-5 | 6.32e-5 | 6.33e-5 | 6.32e-5 |
| **46** (64) SVD | 1.89e-4 | 1.88e-4 | 1.27e-4 | 1.41e-4 |
| **46** SVD/DMT | **2.98** | **2.98** | **2.01** | **2.23** |

Time dependence at `chi' = 46`, cutoff 0, eps 0.25 — the advantage *grows*:

| t | 0.8 | 1.0 | 1.2 | 1.4 | 1.6 | 1.8 |
|---|---|---|---|---|---|---|
| SVD/DMT | 0.96 | 0.87 | 0.97 | 1.34 | 2.08 | **2.98** |

### ULS/SU(3) `d = 3`, t = 1.6 (χ_ref = 1465)

| chi' (maxdim) | DMT | SVD | SVD/DMT |
|---|---|---|---|
| 2 (20) | 1.90e-1 | 1.07e-2 | **0.056** |
| 10 (28) | 1.41e-2 | 1.33e-2 | 0.94 |
| 22 (40) | 8.42e-3 | 1.41e-2 | **1.67** |
| 46 (64) | 3.47e-3 | 3.69e-3 | 1.06 |

(cutoff 0, eps 0.25; the other three cells agree to within 3%.) Same shape — catastrophic at the
floor, DMT ahead once the complement is resolved.

### `d = 2` control, t = 6.0 (χ_ref = 620), chi'-matched

| chi' (maxdim) | DMT | SVD | SVD/DMT | peak SVD/DMT over the run |
|---|---|---|---|---|
| 2 (10) | 1.61e-1 | 1.98e-1 | 1.23 | 1.23 at t = 6.0 |
| 10 (18) | 1.52e-1 | 7.27e-2 | 0.48 | 1.13 at t = 0.6 |
| 22 (30) | 2.69e-2 | 5.81e-3 | 0.22 | 4.03 at t = 3.4 |
| 46 (54) | 1.05e-3 | 2.24e-3 | **2.14** | 3.47 at t = 4.2 (cutoff 1e-10) |

**The control's verdict: `d = 3` is not anomalous.** The same "loses near the floor, wins at
generous `chi'`" pattern appears at `d = 2`, where DMT is independently trusted (this branch
reproduces XXZ exponents 1.02/1.62/1.79). Note the tight-rung `d = 2` numbers are 10–20% errors on
both arms — a ratio between two broken numbers, past front containment, and not worth much.

## 5. What DMT unambiguously buys: invariance

> **Superseded by Addendum 1 §Critical 1 and Addendum 2 §6.** The "0.1% spread" figure below is
> the reference's own error reported four times. The trace and cutoff-floor findings stand.


This is the cleanest signal in the whole run, and it is not visible in a single-configuration test.

**Error invariance.** At `d = 3`, `chi' = 46`, `t = 1.8`, DMT reads `6.326e-5 / 6.320e-5 /
6.326e-5 / 6.321e-5` across all four (cutoff, eps) cells — a **0.1% spread**. Plain SVD swings by
1.5x at the same rung and far more at earlier times.

**Trace.** `max |tr(rho)/tr(rho_0) − 1|` over each run:

| | DMT | SVD |
|---|---|---|
| Heisenberg `d = 3` | 6.4e-14 … 8.9e-14 | up to **3.2e-6** |
| Heisenberg `d = 2` | 6.7e-13 … 7.5e-13 | up to **1.3e-4** |
| ULS `d = 3` | 6.9e-14 … 7.7e-14 | up to 3.6e-7 |

The trace survives because it is protected structurally, not because it wins a weight contest.

**The cutoff floor — the practically important one.** At `eps = 0.01` with the ordinary
`cutoff = 1e-10` *that both sibling melts set*, the SVD arm's error stops depending on `maxdim`
(Heisenberg `d = 3`):

| t | 0.6 | 0.8 | 1.0 | 1.2 | 1.4 | 1.6 | 1.8 |
|---|---|---|---|---|---|---|---|
| SVD `chi' = 22` | 1.73e-5 | 3.44e-5 | 1.02e-4 | 1.73e-4 | 1.40e-4 | 1.51e-4 | 9.23e-4 |
| SVD `chi' = 46` | 1.73e-5 | 3.44e-5 | 1.03e-4 | 1.86e-4 | 2.26e-4 | 2.00e-4 | 1.41e-4 |
| DMT `chi' = 46` | 1.52e-7 | 4.68e-7 | 1.47e-6 | 3.28e-6 | 9.04e-6 | 2.60e-5 | 6.32e-5 |
| SVD/DMT at 46 | **114** | **74** | **70** | **57** | 25 | 7.7 | 2.2 |

Doubling plain SVD's bond dimension changes nothing — its error is even non-monotone in time — because
a cutoff measured against the **full** Hilbert-Schmidt norm cannot see a signal that is 2.8% of it.
DMT keeps improving. This is exactly the infinitesimal-wall linear-response regime the DMT
literature works in.

**A null control that behaved.** At `cutoff = 0` the `eps` rows are flat to <1% at every rung
`chi' ≥ 10` (the identity is a *product* operator occupying one Schmidt direction, so a pure
`maxdim` ranking is `eps`-blind). They are *not* flat at `chi' = 2`, where DMT's error moves 4x
between `eps = 0.25` and `0.01` — with two complement directions the kernel's choice is genuinely
state-dependent.

## 6. Verdict

> **Superseded by the Revised bottom line at the end of Addendum 2.** The `chi' = 46` "3x
> better" below is not resolved by this data, and the mechanism attributed to overhead fraction is
> falsified in Addendum 2 §1.


**Does DMT measurably beat plain SVD truncation at equal bond dimension at `d = 3`?**

**Not at tight budgets — it is measurably worse.** At `maxdim = 20` (`chi' = 2`) DMT's profile
error at `t = 1.8` is `4.15e-2` against plain SVD's `6.42e-3`: DMT is **6x worse**. At
`maxdim = 40` (`chi' = 22`) it is still 2.5x worse (`2.81e-3` vs `1.14e-3`). This is not noise —
it holds at every time from `t = 0.4` onward, in all four (cutoff, eps) cells, and in the ULS model
(0.056 at `chi' = 2`) and the `d = 2` control (0.22 at `chi' = 22`).

**Yes at generous budgets, and increasingly so with time.** At `maxdim = 64` (`chi' = 46`) DMT's
error is `6.33e-5` against SVD's `1.89e-4` — **3x better** — with the ratio climbing 0.97 → 1.34 →
2.08 → 2.98 over `t = 1.2 … 1.8`.

**The mechanism is the fixed structural overhead.** The protected block is `2 d² = 18` directions
regardless of budget: 90% of `maxdim = 20`, but 28% of `maxdim = 64`. Near the floor DMT spends
almost everything on a subspace plain SVD would not have chosen, loses in Frobenius norm (where SVD
is optimal by Eckart-Young), and that error feeds back through the gates. **Practical rule: at
`d = 3`, `maxdim` must be several times `2 d²` before DMT is worth using at all.** The `d = 2`
control confirms this is the method's behaviour, not higher spin's.

**And DMT buys robustness that no `maxdim` buys plain SVD.** Error invariant to 0.1% across wall
amplitude and cutoff; trace exact to `1e-13`; and at an infinitesimal wall with a default cutoff,
plain SVD hits a floor that more bond dimension does not remove while DMT is 70x better.

**What this does not show.** Nothing about `z`, and nothing asymptotic — the whole `d = 3` window
is `t ≤ 1.8`. The tempting extrapolation ("the `chi' = 46` advantage grows with time, so it will be
large at `t = 12`") is precisely the unfalsifiable move this benchmark exists to avoid; testing it
would need an exact reference at `t ≈ 12`, which is not reachable at any `N`.

**Bearing on `spin1_melt.jl` (extrapolation, explicitly labelled).** That run used
`maxdim = 40/56/72`, i.e. `chi' = 22/38/54`, so its ladder straddles the `chi'` range where this
benchmark's verdict changes sign. Two reasons to treat the inference as an extrapolation and not a
result: it transfers a crossover measured at `N = 12, t ≤ 1.8` to `N = 100, t ≈ 12`, which is the
move ruled out immediately above; and `spin1_melt.jl`'s ULS arm is the model in which `chi' = 22`
is on the *winning* side here (1.67x at `t = 1.6`), so the concern applies mainly to its Heisenberg
arm. The recommendation stands — raise the lowest rung — but as a precaution to verify, not a
demonstrated defect.

---

# Addendum: corrections after review

Review found the negative result sound and **under**-sold, and the positive half not supported by
the data. Both are addressed below. **Nothing in the negative result was softened** — with the
reference's own error now measured per case and subtracted, it is 3-6x stronger than first
reported.

## Critical 1 — the `chi' = 46` "win" was largely a measurement of the reference

**Option taken: both.** The floor probe was extended (per case, per time — see Critical 2) *and*
the rung is now restated one-sidedly. The extension is what showed the one-sided restatement is
not merely prudent but mandatory.

Measured `d = 3` Heisenberg floor (reference at `cutoff = 1e-10` vs `1e-12`):

| t | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.2 |
|---|---|---|---|---|---|---|
| floor | 1.57e-9 | 4.36e-8 | 1.01e-7 | 6.11e-7 | 1.71e-6 | 3.03e-6 |
| DMT `chi' = 46` / floor | 1.00 | 1.59 | 1.51 | **0.77** | **0.86** | 1.08 |
| SVD `chi' = 46` / floor | 1.00 | 1.59 | 1.51 | **0.73** | **0.74** | 1.05 |

Both arms are at or *below* the reference's own error at every time the floor is known. The
reviewer's reading is confirmed and is worse than stated: it is not only the first two points of
`0.97 → 1.34 → 2.08 → 2.98` that are floor-pinned — **the entire cutoff-0 `chi' = 46` rung is
unresolved**, and since the probe stops at `t = 1.2` (the `1e-12` reference costs ~2.4x per sweep
and its next sweep projected at 514 s) the `t = 1.6` and `t = 1.8` values `2.60e-5` and `6.33e-5`
cannot be shown to be above it either. On the `d = 2` evidence, floors keep climbing into the
`1e-5` range.

**Withdrawn:** "3x better and the gap grows with time" as a measured kernel-only result; and the
flagship `6.326e-5 / 6.320e-5 / 6.326e-5 / 6.321e-5` "0.1% spread", which is the reference's error
reported four times.

**Replaced by** a guaranteed interval. With `m = |arm − reference|` and floor `f`, the triangle
inequality gives `|arm − exact| ∈ [m − f, m + f]`, so the true ratio lies in
`[(m_S − f)/(m_D + f), (m_S + f)/(m_D − f)]`. The script now prints the lower end per cell and
marks intervals containing 1 as `UNRESOLVED`. `ratio_bounds` implements it.

The positive claim that **does** survive this treatment is the cutoff pathology, where the SVD arm
is 35-170x the floor and solidly resolved while DMT sits at it:

| case | rung | t | guaranteed DMT advantage |
|---|---|---|---|
| Heisenberg `d = 3` | `chi' = 22` | 1.2 | **≥ 2.2x** |
| Heisenberg `d = 3` | `chi' = 46` | 1.2 | **≥ 29x** |
| Heisenberg `d = 2` | `chi' = 22` | 2.6 | **≥ 19x** |
| Heisenberg `d = 2` | `chi' = 46` | 2.6 | **≥ 35x** |

(`eps = 0.01`, `cutoff = 1e-10`.) This is now the positive headline: it is corroborated at both
local dimensions, at two rungs each, with the floor subtracted.

## Critical 2 — floors now measured for all three cases

`reference_floor` took a hardcoded `("heisenberg", 3, …)` and a single `FLOOR_NCALL`. It now takes
a per-case `floor_budget`, returns the **curve** rather than a max, and `run_case` runs it for
every case.

- **`d = 2` (full window, `t = 6.0`, completed):** floor rises to `2.2e-5` at `t = 2.2`, then
  *saturates and drifts down* to `3.0e-5` by `t = 6.0`. At `t = 2.6` it is `1.66e-5` — and the
  `chi' = 46` arms read `1.6154e-5` (DMT) and `1.6184e-5` (SVD). The reviewer's diagnosis is exact:
  **that cell is the reference's error reported twice**, and its `1.002` ratio is that error
  divided by itself.
- **ULS floor (to `t = 1.2`):** 1.5e-9, 2.1e-8, 3.6e-7, 1.8e-6, 4.4e-6, 6.4e-6.

**So the `d = 2` control corroborates the loss at `chi' = 2` and `chi' = 10` and nothing more.**
At `chi' = 22` its guaranteed interval is `[0.21, 1.73]` — straddling 1, no winner. At `chi' = 46`
it is pure floor. It cannot corroborate any crossover, and the report no longer claims it does.

## Important 3 — mechanism corrected

The overhead-*fraction* story is falsified by our own control: `d = 2, maxdim = 30` has 26.7%
overhead, essentially the same as `d = 3, maxdim = 64`'s 28.1%, yet DMT loses there. Rewritten.

The replacement is weaker than "absolute `chi'` decides", because the middle of the ladder is
non-monotonic: `d = 3, chi' = 10` reads 1.37 at `t = 1.6` but 0.69 at `t = 1.8`, and flips between
metrics (`err_inf` 0.689 vs `err_l1` 1.108 at `t = 1.8`). What the data supports: **only the ends
of the ladder are ordered** — `chi' = 2` loses everywhere, `chi' = 46` is never resolved as a loss
— and no single-parameter rule fits the middle.

## Important 4 — `d = 2` requoted at `t = 2.6`

`t = 6.0` is 2.3x past that case's containment, and at `chi' = 2` it reads 1.23 (DMT winning),
contradicting the headline. At the last contained time `t = 2.6` (cutoff 0, `eps = 0.25`) the
control reads **0.076 / 0.023 / 0.673 / 1.002** — far more supportive. Both the file and
`docs/src/manual/dmt.md` now quote `t = 2.6`.

## Important 5 — docs de-generalized

`docs/src/manual/dmt.md` claimed "at `maxdim = 40` (`chi' = 22`) still 2.5x worse" for all of
`d = 3`; that is Heisenberg-only, and **ULS is 1.67x better at that exact rung**. It also asserted
unconditional flatness across wall amplitude. The admonition is rewritten: negative result stated
as guaranteed bounds, no kernel-only win claimed, cutoff-immunity separated out as the resolved
positive claim, and per-rung model dependence stated explicitly.

## Important 6 — the two invariance claims split

"the eps rows should be flat (they are)" was false at the tight rungs. Measured DMT eps-spread at
`cutoff = 0`, `d = 3`, `t = 1.6`: `chi' = 2` → **272%**, `chi' = 10` → 3.7% (28% at `t = 1.8`),
`chi' = 22` → 0.6%, `chi' = 46` → 0.1%.

The clean half now leads: **cutoff-invariance**, quoted at `chi' = 22` where DMT's error
(`1.09e-3`) is three orders above the floor. DMT moves 0.019% / 0.001% / 0.019% between the two
cutoffs at `chi' = 10 / 22 / 46`; plain SVD moves 34% / 180% / 268%.

## Important 7 — `spin1_melt.jl` bearing labelled

The arithmetic stands but it transfers an `N = 12, t ≤ 1.8` crossover to `N = 100, t ~ 12`, which
§6 forbids, and it ignores that `spin1_melt.jl` runs ULS, where `chi' = 22` is on the *winning*
side (1.67x at `t = 1.6`). Kept as a recommendation, labelled extrapolation, with the ULS caveat.

## Important 8 — headline quoted as a range

At the last front-contained time the loss at `chi' = 2` is **2.9x-11.2x** (`t = 1.6`, `eps = 0.01`
to `0.25`); at `t = 1.8` it is 1.7x-6.5x; at `t = 1.2`, with the floor measured and subtracted, it
is **≥ 39x**. The sign is identical in all four cells at every time past `t = 0.2`; only the
magnitude moves. Note the loss *shrinks* with time as the SVD arm's own error catches up.

## The corrected user-facing rule

Replaces "`maxdim` must be several times `2 d^2`", which is a ratio rule and would license
`maxdim ≈ 24-32` at `d = 2` — where this run shows DMT still behind:

> Budget the **absolute complement** `chi' = maxdim − 2 d^(2n)`, not its ratio to the protected
> block. On this observable DMT was behind at `chi' = 2` in every model tested and at `chi' = 10`
> in the `d = 2` control, and was never resolved as a loss at `chi' = 46`. Take `chi' ≈ 40-50`
> (`maxdim ≳ 2 d^2 + 45`) as a **starting point with real uncertainty**: the crossover moved by
> ~2x in `chi'` between two models at the *same* `d` (ULS is ahead already at `chi' = 22`, spin-1
> Heisenberg is not). Treat the equal-`maxdim` DMT-vs-SVD comparison as a convergence check to run
> for your own model, exactly as this script does — not as a constant to look up.

## Minors

- The reference's cutoff (`1e-10`) is **looser** than the arms' (`0`), which is the root cause of
  Critical 1. Now stated at the `FLOOR_CUTOFF` definition and in the "how to read any number" note.
- The front curve is printed per time (with the containment crossing time) rather than only as a
  per-run maximum. It is left out of the CSV to keep the shipped schema reproducible against the
  55-minute dataset; the per-time values are in the run log.
- `chi' = 10` flips verdict between metrics at `t = 1.8` (`err_inf` 0.689, `err_l1` 1.108). Noted
  as part of the non-monotonic middle of the ladder.

## Revised bottom line

The benchmark's result is **a negative one, and it is strong**: at the tightest budget the DMT
kernel admits, DMT is 13x-44x worse than the plain SVD it replaces (guaranteed bounds, three
models, two local dimensions, last front-contained time), and the `d = 2` control shows this is
the method's behaviour and not higher spin's. The positive counterpart is narrower than first
claimed but real and resolved: DMT is immune to a norm-relative cutoff that pins plain SVD at an
error floor no extra bond dimension removes — ≥ 29x at `d = 3` and ≥ 35x at `d = 2`. The
kernel-only crossover at large `chi'` is **not** established by this data and is no longer claimed.

---

# Addendum 2: corrections after re-review

Re-review verified the interval arithmetic (1.6e7 adversarial draws, zero containment violations,
both endpoints attained) and confirmed all eight prior items addressed. Seven items remained. All
are fixed below except one, which is downgraded to "not established" because the data cannot
support it either way.

## 1. `d = 2, chi' = 22` — I reintroduced the error I set out to remove

Correct, and it was load-bearing twice. That cell is **UNRESOLVED at `cutoff = 0`**, interval
`[0.207, 1.729]`, with the arms agreeing to five digits — the same pure-floor signature the first
addendum correctly diagnosed at `chi' = 46`. The only cell in which it resolves is
`cutoff = 1e-10`, `eps = 0.25`, at a marginal `≥ 1.02x`. Asserting it as a loss in five places
(script x3, docs, report x2) was exactly the mistake the floor treatment exists to prevent. All
five are removed.

**Both arguments are re-grounded on resolved cells, and the conclusion gets stronger.** The
decisive pair is at the *same* `d = 3`, *same* `chi' = 22`, *same* `t = 1.2`:

| model | DMT | SVD | floor | interval | verdict |
|---|---|---|---|---|---|
| spin-1 Heisenberg | 7.47e-5 | 3.76e-5 | 3.03e-6 | [0.445, 0.567] | **LOSS ≥ 1.7x** |
| ULS/SU(3) | 9.52e-4 | 1.69e-3 | 6.35e-6 | [1.756, 1.794] | **WIN ≥ 1.7x** |

Same local dimension, same complement budget, same time, both solidly resolved, opposite verdicts.
**The model decides, not the budget.**

- *Overhead-fraction mechanism — falsified on resolved cells.* ULS at `chi' = 10`
  (overhead 64.3% of `maxdim = 28`) is a resolved **win**, `≥ 5.0x` at `t = 0.8` and `≥ 1.7x` at
  `t = 1.2`. `d = 2` Heisenberg at the same `chi' = 10` (overhead only 44.4% of `maxdim = 18`) is a
  resolved **loss** of `≥ 28.9x`. The *higher*-overhead cell is the one that wins. No appeal to the
  marginal `d = 2, chi' = 22` cell is needed.
- *Ratio-rule rejection — **not established**, and now stated as such.* Every cell a ratio rule
  would license is unresolved in this data (`d = 2` at `chi' = 22` and `46`, `d = 3` at
  `chi' = 46`), so there is no resolved counterexample to it, and at `d = 3` the ratio and absolute
  readings coincide at `maxdim = 64` anyway. The rejection of a lookup constant no longer rests on
  falsifying one particular constant; it rests on the pair above, which rules out any budget-only
  rule of *any* functional form.

## 2. `>= 38.7x` relabelled — two clocks, not one

It is the last **floor-covered** time (`t = 1.2`), not the last front-**contained** time. They
coincide for ULS (`t = 1.2`) and `d = 2` (`t = 2.6`) but not for the headline model: Heisenberg
`d = 3` is contained to `t = 1.6`, where the raw loss is 11.2x and the floor is unmeasured. Both
are now quoted with their own labels in the script, the docs and below, and the script's findings
block carries an explicit "two clocks, do not conflate them" note.

## 3. VERDICT TABLE brought under the floor treatment

It printed raw ratios at `times[end]` with no floor — including the withdrawn `2.98` and the
`d = 2` rows at `t = 6.0` where `chi' = 2` reads 1.23 with DMT *winning*. It now prints every cell
at each case's **last floor-covered time**, with the interval verdict, and states which time range
is left unsigned:

```
case            t     cutoff   eps  maxdim chi'  err_inf DMT  err_inf SVD  raw    floor-corrected verdict
heisenberg_d3   1.2   0e+00    0.250  20     2   1.748e-02    4.476e-04    0.0256 LOSS >= 38.7x
...             (reference reached t = 1.8; floor known only to t = 1.2, so t in (1.2, 1.8] is unsigned)
```

## 4. Floor persisted to the CSV

`append_csv` now writes a `floor` column for every row (NaN where the probe did not reach), and the
CSV was **regenerated by a full rerun** rather than patched. Every bound in this report is now
reproducible from the artifact alone: `lo = (svd − floor)/(dmt + floor)`,
`hi = (svd + floor)/(dmt − floor)`.

## 5. Bounds truncated, never rounded; "guaranteed" qualified

Round-half-up on a lower bound is invalid. `_trunc_down` truncates toward 1 at one decimal, so the
quoted figures are now `≥ 38.7` (exact 38.78), `≥ 43.7` (43.74), `≥ 28.9` (28.96), `≥ 2.1` (2.182),
`≥ 34.7` (34.79/34.80). The `ratio_bounds` docstring now also states that `floor = |loose − tight|`
bounds the reference's error only insofar as the tighter run is itself converged — an estimate, not
a certificate. Inflating every floor by 50% (beyond the reviewer's 10%) moves no verdict in the
shipped data, so the signs are robust; the second digit of any bound is not.

## 6. Both interval ends emitted

`lo` alone is vacuous for a resolved loss — the `chi' = 2` headline cell printed `>= 0.0254`, a
bound in the direction the cell does not support. `verdict_label` now emits whichever end is
informative: `WIN >= lo` for a win, `LOSS >= 1/hi` for a loss, `UNRESOLVED` when the interval
contains 1, `no floor` when the reference's error is unknown at that time. The `ratio_bounds`
docstring, which said only `lo` may be quoted, is corrected.

## 7. Model spread widened

Was "~2x in `chi'`". The resolved cells give **at least 2.2x** — ULS is a resolved win already at
`chi' = 10` while spin-1 Heisenberg is a resolved loss at `chi' = 22` — and the true spread is
plausibly larger, since neither model's crossover is bracketed from both sides by resolved cells.
Stated that way rather than as a single figure.

## Cheap items

- **Counterexample to the positive headline, now stated.** ULS at `chi' = 46`, `cutoff = 1e-10`,
  `eps = 0.01` is a resolved **win of ≥ 3.5x at `t = 1.0`** and a resolved **loss of ≥ 1.5x at
  `t = 1.2`**. DMT's cutoff *invariance* holds there (error moves ≤ 0.02% between cutoffs); what
  fails is the inference that invariance must translate into an advantage over SVD at every rung
  and time. In the script and the docs.
- **`err_l1` given the floor treatment** for the metric-flip argument: at `d = 3`, `chi' = 10`,
  `cutoff = 0`, `eps = 0.25`, both metrics are resolved losses at the floor-covered times
  (`t = 1.2`: `err_inf` [0.603, 0.621], `err_l1` [0.751, 0.844]). The flip to `err_l1 = 1.108`
  happens at `t = 1.8`, which has no floor — so the metric flip is *not* a resolved phenomenon and
  is no longer offered as evidence, only as a reason to distrust the unsigned late columns.
- **ULS `chi' = 22` requoted at `t = 1.2`** (guaranteed `≥ 1.7x`) instead of `t = 1.6`, which is
  past both containment and the floor probe.
- **Report lede rewritten** so it no longer states the two retracted claims in bold as "the
  answer", and §5/§6 carry inline superseded markers.

## Revised bottom line (supersedes §6 and Addendum 1)

**The negative result is the finding, and it is strong.** At the tightest budget the DMT kernel
admits, DMT is **≥ 38.7x** (spin-1 Heisenberg), **≥ 43.7x** (ULS) and **≥ 13.0x** (`d = 2` control)
worse than the plain SVD it replaces — guaranteed bounds at each case's last floor-covered time,
same sign in every wall-amplitude and cutoff cell, corroborated at a second local dimension.

**No budget-only rule explains where that reverses.** At the same `d`, `chi'` and `t`, one model is
a resolved loss and another a resolved win. Protected-block fraction is falsified outright.

**The kernel-only advantage at large `chi'` is not established by this data** — both arms sit at
the reference's resolution there.

**The resolved positive claim is narrow but real:** DMT is immune to a norm-relative cutoff that
pins plain SVD at an error floor no extra bond dimension removes (≥ 29.0x at `d = 3`, ≥ 34.7x at
`d = 2`), it holds the trace to `1e-13`, and its error is cutoff-invariant to ≤ 0.02%. Even this is
not universal across rungs and times — ULS at `chi' = 46` reverses within one time step.

## Postscript: a third clock, found by the item-3 fix itself

Rerunning with the corrected verdict table surfaced a bug in that very fix. Summarizing at "the
last **floor-covered** time" is right for `d = 3`, where the floor probe stops first — but at
`d = 2` the floor probe covers the *whole* window, so the rule selected `t = 6.0`, which is 2.3x
past containment, and the table printed `chi' = 2 … WIN >= 1.2x`. That is the identical failure the
previous round's item 4 flagged, reintroduced from the other direction.

There are **three clocks**, and they order differently per case:

| case | reference reaches | floor known to | front contained to | **summarized at** |
|---|---|---|---|---|
| Heisenberg `d = 3` | 1.8 | 1.2 | 1.6 | **1.2** |
| ULS `d = 3` | 1.6 | 1.2 | 1.2 | **1.2** |
| Heisenberg `d = 2` | 6.0 | 6.0 | 2.6 | **2.6** |

The table now summarizes at the **slowest** of the three and prints all four columns so the reader
can see which constraint bound. A ratio without a floor cannot be signed; a ratio past containment
is signed but describes a saturated box rather than a melt.

The regenerated CSV is unaffected — it carries every time, both metrics and the floor, so any
table is derivable from it. The bounds quoted throughout the addenda were already at `t = 1.2`
(`d = 3`) and `t = 2.6` (`d = 2`) and do not move.

**Reproducibility of the rerun.** The regeneration reproduced the first run exactly where it
matters: reference reach `t = 1.8` at `chi = 1433` (Heisenberg `d = 3`), `t = 1.6` at `chi = 1465`
(ULS), `t = 6.0` at `chi = 620` (`d = 2`), and identical floor curves to the digits printed. The
wall-clock guards fired at the same sweep.

## What I could not fix

**The ratio-rule rejection is not established, and I have stopped claiming it is.** Every cell that
a rule of the form "`maxdim` a few times `2 d^2`" would license is unresolved in this data
(`d = 2` at `chi' = 22` and `46`, `d = 3` at `chi' = 46`), so there is no resolved counterexample to
it; and at `d = 3` the ratio and absolute readings coincide at `maxdim = 64`, so this benchmark
cannot separate the two functional forms at all. Establishing it would need a reference tighter than
`1e-12` at `d = 2, maxdim = 30`, which costs more than the whole run.

What replaces it is stronger and does not depend on the functional form: at the same `d`, the same
`chi'` and the same time, two models give **opposite resolved verdicts**, so no budget-only rule of
any shape can be correct. The recommendation is therefore the convergence check, not a constant.

**The `d = 3` floor beyond `t = 1.2`** remains unmeasured; the `1e-12` reference's next sweep was
projected at 514 s (Heisenberg) and 836 s (ULS). Everything in `t ∈ (1.2, 1.8]` is reported
unsigned rather than estimated.
