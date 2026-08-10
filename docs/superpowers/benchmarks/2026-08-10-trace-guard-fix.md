# The `normalize=true` trace guard: a units error, not a conditioning test

**Answer: the guard was comparing the wrong two numbers, and putting them in the same units fixes
it outright — with a provable no-false-positive guarantee that the cancellation test the fix was
specified as cannot have.** `expectations.jl` compared the all-identity *coefficient*
`c_I = tr(ρ)/d^(N/2)` against `sqrt(eps) ||ρ||_HS`. Restoring the `d^(N/2)` makes the comparison
`|tr(ρ)| > sqrt(eps) ||ρ||_HS` in physical units, where `tr(ρ) = Σ_i λ_i ≥ sqrt(Σ_i λ_i²) =
||ρ||_HS` holds for **every** positive operator, so no thermal state can ever be rejected, at any
temperature, chain length or local dimension. Both sides are now evaluated in log space, which
also removes the second failure mode (`norm(rho)` overflowing to `NaN` and silently disabling the
guard). A cancellation-based test was built and measured first, as specified; it is reported in
§4 and was rejected on evidence — it misses the traceless case it exists to catch.

Branch `feat/higher-spin-dmt`. Everything below is measured; the probe scripts are throwaway and
their contents are reproduced inline where the construction matters.

---

## 1. The defect, reproduced

`ρ = ⊗_j e^{-β Z_j}` is positive, bond dimension 1, exactly representable, and has trace exactly
`(2 cosh β)^L`. Its single-site expectation is `<Z_x> = -tanh β` at every site and every `L`.

At `HEAD = bcc46ba`, `operator_expectation(rho, Z, x; normalize=true)`:

| β | L=120 | L=240 | L=400 |
|---|---|---|---|
| 0.35 | ok | ok | **throws** |
| 0.50 | ok | **throws** | **throws** |
| 1.00 | **throws** | **throws** | ok — *and that is worse* |

with `tr(ρ)/‖ρ‖_HS` measured at 1.4e6 … 7.9e50 across the grid. Nothing is remotely
ill-conditioned. The shipped condition was

```julia
abs(denominator) <= sqrt(eps(Float64)) * norm(rho)      # denominator = c_I = tr(ρ)/d^(N/2)
```

i.e. `tr(ρ)/‖ρ‖_HS ≤ sqrt(eps) · d^(N/2)`. The threshold therefore *grows* as `2^(L/2)`
(2.4e52 at `L = 400`) while the Gibbs ratio grows only as `(cosh β / sqrt(cosh 2β))^{-L}`. Which
side wins is a function of chain length, which is why the failures form a staircase.

The `(β, L) = (1.0, 400)` cell is the second failure mode. `‖ρ‖_HS = 1e175`, so `‖ρ‖²  = 1e350`
overflows; `norm(rho)` returns `NaN`; `x <= sqrt(eps) * NaN` is `false`; the guard does not fire.
**Simultaneously too strict at `L = 120` and unenforced at `L = 400`.**

## 2. Why the naive repairs fail

- **Delete it.** Restores the hazard it exists for: a traceless operator whose post-truncation
  trace is an `O(eps)` cancellation residue, where dividing amplifies every entry by `~1/eps`.
- **Loosen the tolerance.** The mismatch is `d^(N/2)`, not a constant. Any fixed threshold is
  wrong at some `L`.
- **Compare against `prod_j ‖A_j‖`.** Measured in §4 as "candidate A". It decays exponentially in
  `L` for legitimate positive states with `chi > 1` — `10^-33` for an `L = 80` XXZ Gibbs state.

## 3. The fix

`|tr(ρ)| > sqrt(eps) ‖ρ‖_HS`, evaluated as

```
(N/2) log d + log|c_I|  >  lognorm(ρ) + log(sqrt(eps))
```

Two new documented helpers in `src/operator_space/expectations.jl`: `_log_trace_resolution`
returns the pair `(log|tr ρ|, log ‖ρ‖_HS)`, and `_reject_unresolvable_trace` compares them. The
`normalize=false` path is untouched (the check is still behind `normalize &&`, so the norm is
never computed there), and `pauli_*` remain plain delegations.

**Why this is the right discriminator.** For positive semidefinite `ρ` with eigenvalues `λ_i ≥ 0`,
`tr(ρ) = Σ λ_i ≥ sqrt(Σ λ_i²) = ‖ρ‖_HS`, with equality only at rank 1. So a positive operator
clears the threshold by the full `1/sqrt(eps) = 6.7e7`, *independent of `N`, `d` and `β`. That is a
theorem, not a tuned margin.* What decays exponentially in `L` for a cold Gibbs state is the
identity **component** `c_I`; its **trace** does not, and the fixed check now reads the trace.

**Why log space is load-bearing, not cosmetic.** The `(β, L) = (1.0, 400)` state has
`‖ρ‖_HS = 1e175` and `|tr ρ| = 5.9e195`. `lognorm` never forms `‖ρ‖²`. A non-finite trace or norm
is raised as an `ArgumentError` in its own right rather than compared — `NaN <= x` is `false`, and
a `false` here is a silent pass.

### Measured: positive operators (must not throw)

`log10(|tr ρ| / ‖ρ‖_HS)`; the rejection line is `-7.83`.

| state | ratio (log10) |
|---|---|
| `⊗ e^{-0.35 Z}`, L = 120 / 240 / 400 | 15.27 / 30.54 / 50.90 |
| `⊗ e^{-0.5 Z}`, L = 120 / 240 / 400 | 13.02 / 26.04 / 43.39 |
| `⊗ e^{-1.0 Z}`, L = 120 / 240 / 400 | 6.14 / 12.28 / **20.47** |
| `⊗ e^{-β Z}` rescaled, β = 5 / 20 / 100, L = 400 | 0.01 / 0.00 / **0.00** |
| XXZ Gibbs (`Jz = 0.8`), β = 0.5 / 2 / 6, L = 60, chi ≤ 64 | 8.44 / 4.14 / 1.27 |
| `d = 3`, `⊗ e^{-3 h}`, L = 200 | 4.32 |
| `d = 4`, `⊗ e^{-3 h}`, L = 200 | 11.80 |
| identity + 5% random noise, L = 80, chi = 32 | 12.04 |

The `β = 100, L = 400` row is the theorem's equality case (rank-1 local blocks) and lands at
exactly `1.00`, 7.83 decades clear of the line. Every one of the nine reviewer grid cells now
returns `-tanh β` to `1e-12` against the exact oracle.

### Measured: traceless operators (must throw)

| state | ratio (log10) | verdict |
|---|---|---|
| bare `pauli_domain_wall_state`, L = 12 / 24 / 64 | `-Inf` | throw |
| DMT-melted wall, L = 12, sweeps 1 / 6 / 12 | -13.17 / -11.29 / -10.71 | throw |
| DMT-melted wall, L = 20, sweeps 1 / 6 / 12 | -11.81 / -10.52 / -10.40 | throw |
| DMT-melted wall, L = 40, sweeps 1 / 6 / 12 | -9.18 / -8.50 / **-7.48** | throw / throw / **pass** |
| TEBD-evolved single `Z`, L = 16, steps 1 / 12 | -15.97 / -11.41 | throw |
| TEBD-evolved single `Z`, L = 32, steps 1 / 12 | -13.56 / -9.12 | throw |

**The one acceptance is honest and is documented in the docstring.** At `L = 40, maxdim = 24`,
twelve sweeps of DMT have grown the residue to `3.3e-8` relative — genuinely above half
precision. The criterion is a statement about how many digits the trace still has, not about
tracelessness, and it cannot be tightened past `sqrt(eps)` without eating into the positive
operators' guaranteed floor of `1.0`. Traceless transport operators should carry
`normalize=false` by construction; the guard is a backstop, not a classifier.

### Measured: the threshold itself

An indefinite operator has no positivity floor, so it is what actually exercises the line.
`ρ = ⊗_j diag(1, -a)` has the closed forms `|tr ρ|/‖ρ‖_HS = ((1-a)/sqrt(1+a²))^L` and
`tr(ρ A_1)/tr(ρ) = tr(A²)/tr(A)`. At `a = 0.44` the per-site factor is `0.51258`:

| L | 6 | 12 | 20 | 26 | 28 | 34 |
|---|---|---|---|---|---|---|
| `\|tr ρ\|/‖ρ‖_HS` | 1.8e-2 | 3.3e-4 | 1.6e-6 | 2.8e-8 | 7.5e-9 | 1.4e-10 |
| new form | pass | pass | pass | pass | throw | throw |
| old form (÷ `2^(L/2)`) | pass | pass | **throw** | **throw** | throw | throw |

`L = 20` and `L = 26` are the discriminating pair: two to four decades of resolvable trace that
the old comparison threw away on a units error. The returned value at `L = 26` matches
`tr(A²)/tr(A) = 2.13142857…` to 1e-15, so a small trace-to-norm ratio is not by itself an accuracy
problem for an exactly representable state — which is the whole point.

## 4. The cancellation test that was specified, built, and rejected

The fix was specified as: accumulate running magnitudes during the right-environment sweep and
flag when `|c_I|` falls far below the no-cancellation scale. Two forms were implemented in log
space with per-step rescaling (so nothing overflows) and measured against the same battery.

**Candidate A — `|c_I| / prod_k ‖M_k‖_F`**, where `M_k = <I| ρ[k]` is the site tensor's identity
slice. Exactly `1` for a bond-dimension-1 product state at every `β` and `L`. But
`‖M_1 ⋯ M_N‖ ≪ prod ‖M_k‖` generically once `chi > 1` — the submultiplicative bound is loose by a
factor per site — so it collapses on legitimate positive states: `10^-11.5` for XXZ at `L = 40,
chi = 16`, `10^-33.1` at `L = 80, chi = 64`, with no sign cancellation anywhere. Same
exponential-in-`L` defect as the guard being replaced. Rejected.

**Candidate B — `|c_I| / A_1` with `A_k = |M_k| A_{k+1}`** (elementwise absolute values: the value
the contraction would take if nothing cancelled at any level). This is the tight, componentwise
form, and it is *excellent* on positive operators — `10^0.00` for the product Gibbs states at
every `(β, L)`, `10^-0.09` for XXZ at `β = 2, L = 80`, worst case `10^-2.28` on identity-plus-noise
at `chi = 32`. **But it does not detect the case it exists for:**

| traceless state | candidate B ratio (log10) | verdict at `sqrt(eps)` |
|---|---|---|
| DMT-melted wall, L = 12, sweep 1 | **-0.65** | **pass** |
| DMT-melted wall, L = 12, sweep 3 | **-1.20** | **pass** |
| DMT-melted wall, L = 20, sweep 1 | **-0.88** | **pass** |
| TEBD-evolved `Z`, L = 16, step 1 | **-0.22** | **pass** |
| DMT-melted wall, L = 12, sweep 6 | -9.79 | throw |

The reason is structural, and it is why no amount of tuning rescues the approach: for a freshly
evolved traceless operator there **is no cancellation in the environment sweep**. The identity
slices `M_k` are each individually small (`‖M_k‖ ≈ 0.07` per site at `L = 12`), and their product
is small because it is a product of small numbers — the magnitudes multiply down to `1e-13.8` with
no sign structure at all. The cancellation that makes the trace a residue happened earlier, inside
the gate application and the truncation, in the operator's own construction history, which
`operator_expectation_profile` cannot see. Confirming this with a forward-error estimate: the
componentwise relative error of `c_I` for that state is `8e-15` — the residue is computed to
fourteen digits; it is simply a residue.

So: **a cancellation-based test cannot be made reliable here, and the fallback the specification
named — remove the guard, document the hazard, add an opt-in check — was not needed either.** The
units fix is a third option that is strictly better than both: it is the direct remedy for the
diagnosis (an implicit `d^(N/2)`), it keeps the guard automatic, and it is the only one of the
three with a proof that positive operators are never rejected.

## 5. The `_dmt_connector` sibling: **not** affected

`src/operator_space/dmt/lowrank.jl:49` tests `|q0' S r0| > sqrt(eps) * norm(bond_matrix)`. The
same *shape*, but the scales are matched, for one specific reason: `q0` and `r0` are passed
through `_unit_direction` (`bond.jl:258-259`) before the comparison. Writing `L_I` and `R_I` for
the unnormalized identity directions of the two protected blocks,

```
|q0' S r0| = |c_I| / (‖L_I‖ ‖R_I‖)          ‖S‖_F = ‖ρ‖_HS   (canonical gauge at the bond)
```

and `‖L_I‖`, `‖R_I‖` decay with the number of sites they span at exactly the rate `c_I` does. The
normalization divides the exponential out. The `expectations.jl` guard had no such normalization:
it compared the raw `c_I` (which *is* `L_I' S R_I` with unnormalized environments) against a
whole-operator norm.

Measured, replicating lines 193-260 of `bond.jl` and taking the worst bond of each chain:

| state | worst `\|q0' S r0\|/‖S‖_F` |
|---|---|
| XXZ Gibbs β = 0.5, L = 20 / 60 / 120 | 0.9335357449853706 / …707 / …707 |
| XXZ Gibbs β = 2.0, L = 20 / 60 / 120 | 0.42369526816915 / 0.42369554568706 / 0.42369554568706 |
| XXZ Gibbs β = 6.0, L = 20 / 60 / 120 | 0.29166709 / 0.29346445 / 0.29346447 |
| near-product Gibbs β = 2.0, L = 20 / 60 / 120 / 240 | 0.517 / 0.402 / 0.629 / 0.581 |
| DMT-melted wall, L = 12 / 20 / 40 (largest bond) | 7.9e-14 / 2.9e-14 / 3.1e-15 |

The XXZ rows are the decisive signature: they agree to 16 significant figures (β = 0.5) or 8
(β = 6.0, where the state has not yet saturated its correlation length at `L = 20`) across a 6x
change in `L`, and to 14 figures between `L = 60` and `L = 120`. The ratio has **no `L`
dependence** once the state is converged, so there is no exponential to be swamped by, and the
traceless rows still land at `1e-14`. No change made, none needed. (The weakest
statement that survives as a *bound* rather than a measurement is `|q0' S r0| ≥ |c_I| ≥
d^{-N/2}‖S‖_F`, since `‖L_I‖, ‖R_I‖ ≤ 1`; that is only the worst case and the measurements are 13
to 60 decades better than it. If this is ever revisited, the thing to check is whether some
`ρ` can drive `‖L_I‖‖R_I‖` far above `|c_I|/‖ρ‖`, not the `L` scaling.)

## 6. Verification

- `julia --project=. -e 'using Pkg; Pkg.test()'` — **passed**, 9m04s wall, foreground, no
  new failures and no skips.
- 75 new assertions in `test/test_operator_space.jl`, testset
  *"normalized expectations reject an unresolvable trace and nothing else"* (+4.7s):
  the nine-cell reviewer grid against the `-tanh β` oracle; the `L = 400` overflow cell
  (`2·lognorm(ρ) > log(floatmax)`, i.e. `norm(rho)` is `NaN`, yet the call succeeds and is
  correct); `_log_trace_resolution` pinned against `operator_trace`/`norm` at `d = 2` and `d = 3`
  (this is the regression test for the `(√d)^N`); the six-point indefinite straddle of §3 with
  its closed-form verdict *and* closed-form value; the bare and TEBD-melted domain wall; two
  non-finite operators that must raise rather than pass; and `pauli_* == operator_*` on both a
  passing and a throwing input.
- The pre-existing traceless rejection in `test_pxp.jl:545` (PXP energy correlator) still throws.
- `normalize=false` is unchanged everywhere, including the DMT preservation and transport-reference
  suites, which measure exclusively on that path.

## 7. Residual concerns

1. **The guard is a conditioning test, not a tracelessness test.** §3 documents the one measured
   acceptance (`L = 40`, twelve DMT sweeps, residue grown to `3.3e-8` relative). Callers evolving
   traceless operators should pass `normalize=false` rather than rely on the backstop. The
   docstring says this.
2. **Positivity is assumed, never checked.** The `tr ρ ≥ ‖ρ‖_HS` guarantee covers thermal states
   and anything else positive. An *indefinite* operator with a legitimately tiny trace-to-norm
   ratio (below `sqrt(eps)`) is still rejected — correctly, since normalizing it really does lose
   more than half the digits, but it is a rejection with no positivity escape hatch. §3's
   straddle table is exactly this case, deliberately.
3. **`lognorm` cost.** Identical to the `norm(rho)` it replaces (one `loginner`), and only on the
   `normalize=true` path. Not a regression, but it remains `O(N chi³ d²)` per call, which matters
   for the diameter sweeps that call the profile thousands of times — those all use
   `normalize=false` and never reach it.
4. **`lognorm` on an already-orthogonalized MPS** takes ITensorMPS' `log(norm(M[oc]))` shortcut,
   which overflows once `‖ρ‖_HS > 1e308` (measured: `Inf` for the `β = 1, L = 800` state after
   `orthogonalize!`). At that point the centre tensor itself holds `Inf` entries and the operator
   is not representable at all, so the resulting `ArgumentError` is the correct outcome — but it
   is an error raised on the *norm*, not on the trace, and the message says so.

---

## Addendum (2026-08-10, post-review): the guarantee is about conditioning, not about range

**The fix is endorsed and unchanged.** An independent review reproved `tr(ρ) ≥ ‖ρ‖_HS` over 5000
random PSD trials with zero counterexamples and equality iff rank ≤ 1, cross-checked the log-space
quantities against its own dense oracle, forced both the `NaN` and `Inf` sub-modes and confirmed
the explicit branch catches each, verified `normalize=false` and the `pauli_*` delegations are
untouched, and reproduced the `_dmt_connector` verdict on XXZ Gibbs states at `L = 20/60/120` to
~9 significant figures. It also reimplemented the cancellation test from a different construction
— plain TEBD rather than the DMT melt of §4 — and got **0.134 / 0.133 / 0.285 / 0.044** over
sweeps 1-4 (no detection at all) and then **2.8e-15 / 1.6e-14** at sweeps 5-6. Erratic on exactly
the case it exists for, which is a sharper statement of §4's finding than §4 makes: the test does
not merely miss early, it has no monotone relationship to the residue.

**One thing in the write-up was wrong, and it was prose, not code.** Three sentences added by this
change claimed the guard never rejects a thermal state, full stop:

- `docs/src/manual/operator-space.md` — "no thermal state is ever rejected, at any temperature or
  chain length"
- `expectations.jl` `_reject_unresolvable_trace` — "independent of `N`, `d` and temperature"
- `expectations.jl` `operator_expectation_profile` — "always satisfies `tr(ρ) ≥ ‖ρ‖_HS` and is
  never rejected"

That is true of the *comparison* and false of the *computation*. `identity_coefficient`
(`denominator`, the scalar of the identity-environment sweep) is an ordinary linear-space
contraction with no progressive rescaling, unlike the norm side which goes through `lognorm`. It
can therefore overflow before the guard compares anything. For `ρ = ⊗_j e^{-β Z_j}` the identity
coefficient is exactly `(√d cosh β)^N`, and bisecting at `d = 2`:

| N | last finite β | `c_I` there | first non-finite β | `c_I` there |
|---|---|---|---|---|
| 400 | 2.10 | 1.538e307 | 2.11 | `Inf` |
| 200 | 3.89 | 6.522e307 | 3.90 | `Inf` |

(It becomes `Inf` at the crossover and `NaN` a little further past it, once an `Inf × 0` appears in
the contraction; the review saw the `NaN` at `β = 2.13, N = 400`, closed-form `log10 c_I = 312.2`.
Both are caught by the same `isnan(…) || … == Inf` branch.)

This is **a different mechanism from residual concern 4** in the original write-up: no prior
`orthogonalize!` is involved, it hits a fresh ungauged MPS, and it is on the trace side rather
than the norm side. It is **pre-existing and non-regressing** — the pre-fix code hit the identical
overflow and then returned `NaN` silently, through exactly the `NaN <= x` mechanism this change
closes — so no code change was made. It fails loudly with a correctly-labelled representability
error.

**Scoping, and the escape hatch.** The strong claim survives where it is true: a positive operator
is never rejected *as ill-conditioned*. The three sentences now say that, and each points at the
remedy. `normalize!` removes the failure completely, because positivity confines a normalized
identity coefficient to `[d^(-N/2), 1]`:

| normalized `⊗_j e^{-β Z_j}`, N = 400 | `c_I` | `<Z>` returned | exact |
|---|---|---|---|
| β = 2.2 | 1.6e-58 | -0.975743 | -0.975743 |
| β = 20 | 6.2230152778609e-61 | -1.0 | -1.0 |

The `β = 20` row sits on the `2^-200 = 6.2230152778611e-61` floor to 3e-14 relative — the rank-1
saturation of the same bound — and the expectation is still exact there. That interval stays
inside `Float64` for `N ≲ 2048 / log2(d)`, the ceiling the manual already documents for the
`(√d)^N` prefactor. `operator_gibbs_state` calls `normalize!` internally (measured: `‖ρ‖ = 1.0`,
`c_I = 8.5e-9` for XXZ at `β = 8, L = 60`) and so cannot meet this; an operator built directly
with `operator_product_state` / `operator_basis_state` should be normalized before measurement,
which the manual now says in a `!!! warning` block.

**Test output is pristine again.** `lognorm`'s realness check emitted two
`log(norm²) is NaN + NaN*im` warnings on the deliberately poisoned fixtures. The four assertions
are now wrapped in `Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger())` — `Base`, so no
new test dependency. Warning count in a full `test_operator_space.jl` run: **0**.

**Re-verification** (foreground; executable code unchanged, so the full suite was not re-run):

- `test/test_operator_space.jl` + `test/test_docstrings.jl` — **passed**, 2m04s, 0 warnings.
  The new testset is now 81 assertions (was 75), the six added ones being the overflow crossover
  and the normalized-operator immunity above.
- One assertion had to be loosened during this pass and it is worth recording why: `d^(-N/2)` is
  an *exact-arithmetic* floor, and the `β = 20` state saturates it to 2 ulp *below*, so
  `2.0^-200 <= abs(c_I)` fails. The test now allows `1e-9` relative slack — enough for the ulps,
  far too little to hide a coefficient that had actually fallen through the floor.
