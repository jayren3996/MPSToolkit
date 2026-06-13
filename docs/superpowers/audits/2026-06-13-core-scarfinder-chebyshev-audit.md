# Core / ScarFinder / Chebyshev / Dependency Audit

**Date:** 2026-06-13
**Scope:** The parts of MPSToolkit.jl the two prior operator-space audits did *not* reach —
core evolution (`src/evolution/`), ScarFinder (`src/scarfinder/`), observables
(`src/observables/`), Chebyshev/KPM (`src/chebyshev/`), the dense models (`src/models/`), the
DAOE/FDAOE *algorithm* (`src/operator_space/daoe.jl`, previously checked only for numerics), the
DMT *truncation algorithm* vs literature, the cross-cutting `nstep` abstraction, and
dependency/compat hygiene.
**Method:** evidence-first. Every bug claim has a runnable reproduction or an exact-formula/ED
oracle; every non-trivial finding was adversarially re-checked (try to refute, keep what
survives). Four parallel investigators (Chebyshev, DAOE-vs-literature, dependency coupling, plus
the author's core/ScarFinder thread) cross-checked each other. Every executable change was gated
on the full suite (`julia --project=. test/runtests.jl`).
**Prior work built on (not repeated):**
`docs/superpowers/audits/2026-06-13-operator-space-transport-audit.md` (37 findings) and
`docs/superpowers/audits/2026-06-13-transport-review.md` (benchmarking review).

---

## Summary

| | count |
|---|---|
| Confirmed defects | 2 (both fixed) |
| Architectural recommendations (surfaced, not silently applied) | 2 |
| Dependency/sustainability recommendations | 3 |
| Coverage gaps closed this pass | 2 (DAOE weight oracle; Chebyshev factor-of-2 oracle) |
| Subsystems verified correct (no defect found) | core TEBD/TDVP, DMT algorithm, DAOE/FDAOE, Chebyshev/KPM, observables, models |
| Prior-audit conclusion refuted | 1 (docs-8: DAOE "changes the trace") |

**Baseline:** full suite green at start (1175 assertions, exit 0). **After this pass:** green,
**1213 assertions** (+38: +3 ScarFinder regression, +32 DAOE oracle, +3 Chebyshev oracle), exit 0.

**Verdict.** The never-audited core (TEBD/TDVP) is the strongest part of the package: its
correctness is already pinned by independent ED oracles (Trotter sign, 2nd-order Strang, the TDVP
`t = -i·time` convention, periodic-wraparound ordering). The DMT, DAOE/FDAOE, and Chebyshev
algorithms are faithful to their source papers. The two real defects are both low-severity,
latent, and one-line: a ScarFinder field-drop and a documentation inaccuracy. The most useful
*non-bug* outputs are (a) the `nstep` "effective step count" abstraction is genuinely inconsistent
and should be reworked, and (b) the package has no CI test job, leaving the unexported-symbol
coupling and the ITensorMPS 0.3/0.4 straddle unverified.

---

## Confirmed defects (both fixed, test-gated)

### D1 (low) — ScarFinder silently drops `normalize` when upgrading a single-step `DMTGateEvolution`

`_scarfinder_evolution` (`src/scarfinder/algorithm.jl:144-153`) treats an evolution whose
effective step count is `1` as "unsupported", warns, and rebuilds it with `nstep=10` via
`_scarfinder_rebuild_evolution`. The `DMTGateEvolution` branch
(`src/scarfinder/algorithm.jl:92-105`) was written (commit `e8a7422`, Mar) **before** the
`normalize` field was added to `DMTGateEvolution` (commit `b8e3155`, Jun) and was never updated,
so the rebuild omitted `normalize` and it reverted to the constructor default `true`. The sibling
`LocalGateEvolution` and `TDVPEvolution` rebuilds carry every field.

Consequence: a `DMTGateEvolution(...; normalize=false)` — built precisely to track the conserved
unnormalized trace of a traceless operator — passed through **any** ScarFinder entry point at the
default `nstep=1` had its trajectory silently re-normalized, destroying the absolute-trace
diagnostic the field exists to preserve. This is the exact failure mode the prior review's A2-1 fix
addressed on `constrained_dmt_evolve!`; this is the second instance of the same root cause (stale
rebuild constructor).

**Reproduction (before fix):**
```julia
using MPSToolkit
evo = DMTGateEvolution(rand(ComplexF64,4,4), 0.01; schedule=[1], maxdim=6, connector_buffer=6, normalize=false)
MPSToolkit._scarfinder_evolution(evo; warn=false).normalize   # => true   (DROPPED; should be false)
# TDVP control, same scenario:
tevo = TDVPEvolution(rand(2,2), -0.1im; nsteps=1, normalize=false)
MPSToolkit._scarfinder_evolution(tevo; warn=false).normalize  # => false  (correct)
```

**Fix:** add `normalize=evolution.normalize` to the DMT rebuild. Regression test in
`test/test_core.jl` ("scarfinder single-step upgrade preserves evolution fields") asserts both the
DMT and TDVP rebuilds preserve the field. Commit `d6b760f`.

### D2 (low, docs) — `daoe.md` claims the DAOE/FDAOE projectors "change the trace" — they do not

`docs/src/manual/daoe.md:156` stated the projectors are "one-sided and change the trace … not
trace-preserving." This is wrong, and it inherited the wrong rationale from the prior audit's
docs-8. The DAOE/FDAOE channels damp each Pauli string by its operator *size*; the trace lives
entirely in the **weight-0 identity component**, which is never damped, so `tr(O)` is exactly
invariant. The genuine distinction (and the part of the guidance that was correct) is that these
are non-idempotent *size-damping channels*, not Hermitian `PρP` sandwiches, so they reduce the
Hilbert–Schmidt norm of the high-weight tail and are **not** drop-in constraint projectors for
`constrained_dmt_evolve!`.

**Reproduction:**
```julia
using MPSToolkit, ITensorMPS
sites = pauli_siteinds(4)
rho = pauli_basis_state(sites, [:I,:I,:I,:I]; coefficient=3.0)
for (l,c) in (([:X,:X,:I,:I],1.5),([:Z,:Y,:X,:I],0.8),([:X,:Y,:Z,:X],0.4))
  global rho = add(rho, pauli_basis_state(sites,l; coefficient=c); cutoff=0.0)
end
P = pauli_daoe_projector(sites; lstar=1, gamma=0.7)
pauli_trace(rho), pauli_trace(apply(P, rho; cutoff=0.0))   # => 12.0, 12.0  (|delta| = 0)
```

**Fix:** reworded the bullet (size-damping channel, trace preserved, not idempotent → not a `P_G`
drop-in). Backed by a new `pauli_trace`-invariance regression for both channels. Commit `168c77d`.

---

## Architectural recommendations (surfaced — your call, not silently applied)

### R1 (medium) — the `nstep` "effective step count" abstraction is inconsistent

The brief flagged three independent smells around `nstep`; investigation confirms they share one
root cause and one is the D1 bug. The remaining design issue is genuine and worth reworking:

- **Every evolution constructor defaults `nstep`/`nsteps` to `1`** (`types.jl:49,125,219`) and
  treats `1` as a perfectly valid count — `dmt_evolve!`, `tebd` `evolve!`, and
  `constrained_dmt_evolve!` all run correctly at `nstep=1`.
- **ScarFinder unilaterally redefines `1` as "unsupported"** (`algorithm.jl:146-152`) and rewrites
  it to `10`, emitting `@warn` at `algorithm.jl:150`. Because the override fires on the *default*
  configuration, the warning is **expected noise in normal use** — it fires repeatedly in the
  passing suite (5 sites: `scarfinder step order`, `scarfinder loop refinement`, the NaN-warning
  test setup, `finite tdvp evolve`), and the behavior is itself codified by the
  `scarfinder upgrades single-step evolutions` testset. So the warning is *not* a misconfiguration
  signal; it is the design firing as written.
- This conflates two distinct counts: `niter` (macro evolve+project cycles in `scarfinder!`) and
  `nstep` (schedule traversals per `evolve!`). Total work is `niter × nstep`; the override makes
  even `niter=1` silently do 10 traversals, i.e. **10× the evolution time the caller asked for**.

**This is not a correctness bug** (it does what it says, with a warning, and is tested) **but it is
a modeling inconsistency**: a shared field defaulted to `1` everywhere is declared invalid by one
consumer. The D1 field-drop is a direct consequence of the rebuild machinery this override needs.

**Recommendation (architectural — do not apply blind):** either (a) drop the override and let
`niter` control the cycle count (most consistent — `nstep=1` is a valid request), or (b) if a
minimum evolution-per-projection is genuinely wanted, make it an explicit ScarFinder keyword
(`min_substeps`) rather than a silent rewrite of the caller's evolution object. Either removes the
default-config warning and the rebuild-drop footgun class. The two benign smells —
`project_every` inert at `nstep==1` (`constrained.jl:56-58`, already documented as info A2-2) and
the `nstep=1` default itself — resolve naturally under (a).

### R2 (low) — ScarFinder has no real-MPS integration coverage; selectors tested only on mocks

All ScarFinder tests in `test_core.jl` drive mock state types (`DummyState`, `TraceState`,
`StepCountState`, …) with hand-defined `evolve!`/`project!`/`score` methods. This pins the
*control flow* (step order, refinement selection, step-count handling) precisely and cheaply, but
nothing runs `scarfinder!`/`trajectory_refine!` on an actual `MPS` end-to-end, and the real
`EntropySelector`/`FidelitySelector` scoring is exercised only by the standalone observable tests,
never inside a refinement loop. **Recommendation:** add one small real-`MPS` ScarFinder test (e.g.
a few projected TFIM steps with an `EntropySelector`) so the integrated path has a smoke oracle.
Low priority; no defect implied.

---

## Dependency & sustainability (one fixed, rest recommendations)

Installed/tested at audit time: **ITensors v0.9.30, ITensorMPS v0.4.1, KrylovKit v0.10.3, Julia
1.12**; suite green incl. the `_DMT_VERIFY_ENVS` env-oracle.

- **DEP1 (fixed) — missing `julia` compat bound.** Added `julia = "1.10"` (commit `8253d7b`). Both
  hard deps declare `julia = "1.10"` (verified by reading the installed Project.tomls), the source
  uses no `>=1.11` syntax, and 1.10 is the current LTS — the safest defensible floor.
- **DEP2 (recommendation, medium) — no CI test job.** The only GitHub workflow builds docs
  (`documentation.yml`, Julia `"1"`). There is no `Pkg.test()` job and no version matrix, so the
  perf-1 cache's correctness (which depends on the **unexported** `ITensorMPS.leftlim`/`rightlim`
  convention — `dmt.jl:176-177,535-536`, the single highest-risk coupling) and the ITensorMPS
  **0.3/0.4 straddle** (`Project.toml` allows `"0.3.44, 0.4"` but only 0.4.1 is installed/tested)
  are unverified automatically. Recommend a CI test matrix pinning one job to 0.3.x and one to
  0.4.x; consider narrowing to `"0.4"` unless a concrete consumer needs 0.3.x.
- **DEP3 (note) — the brief's "exact pin" premise is incorrect.** `ITensors = "0.9.25"` is a
  caret floor (`>=0.9.25, <0.10`), not an exact pin; v0.9.30 resolves and tests green. No action.
- **Internal-API inventory.** Genuinely non-exported symbols used: `leftlim`, `rightlim`
  (`dmt.jl`, perf-1-load-bearing, qualified `ITensorMPS.`), and `set_nsite!`
  (`chebyshev/moments.jl:222`, in a hand-rolled `ProjMPO`+`position!` Lanczos sweep — the most
  fragile non-DMT coupling). The package consistently module-qualifies the fragile symbols, which
  is good practice and keeps the coupling auditable. The 0.4.0 sole breaking change (dropping ODE
  helpers) does not affect MPSToolkit.

---

## Per-subsystem correctness & coverage assessment

### Core evolution — `src/evolution/` (never previously audited) — **CORRECT, well-oracled**

The foundation is in better shape than its line count suggests, because `test/test_oracles.jl`
already pins the hard parts against independent dense references:
- **TEBD** matches `exp(-iHt)` for TFIM N=4 and shows the 4× error drop under dt-halving (2nd-order
  Strang), validating the `-im` Trotter sign and the `_dense_local_operator` reshape convention.
  The span≥3 reshape ordering is independently pinned by `test_finite_tebd.jl`
  (`kron(id2,x,id2)`, `kron(projector_up,x,projector_up)`), and the periodic-wraparound ordering
  (`kron(op_N, op_1)`) by `test_oracles.jl`.
- **TDVP** matches `exp(-iHt)` with the opposite sign explicitly *not* matching, pinning the
  `t = -i·time` convention. I additionally verified the documented `nsteps=nothing &&
  nsweeps=nothing` fallback (derives the count from `time_step`) runs without error.
- **`energy_density`** (dense path) I verified against a statevector oracle for a generic Hermitian
  operator, single- and multi-window (N=2,3): exact match. The MPO path is tested in
  `test_finite_tdvp.jl`.

Coverage gaps (low/enhancement): TDVP's `nsweeps`-only path, custom `updater`, and the
`reserved_solver_keys` throw are thinly covered; `nsweeps` is notably **absent** from the
`reserved_solver_keys` list (`types.jl:228`) although it is a dedicated field — a passed
`solver_kwargs=(nsweeps=…,)` would forward to `tdvp` unguarded (very minor; the resolved count
already uses the dedicated field).

### ScarFinder — `src/scarfinder/` (largest algorithm) — **D1 fixed; R1/R2 surfaced**

Control flow is precisely tested via mocks (see R2). The energy-matching loops (`_match_energy_*`,
proportional rollback on overshoot/no-progress) are covered by real-MPS tests with sensible
oracles. The one defect (D1) and the design issue (R1) are above. The `nstep→10` magic constant
(`_scarfinder_default_steps`) is arbitrary and undocumented as a physical choice.

### DMT truncation algorithm — `src/operator_space/dmt.jl` — **faithful to White's DMT**

Read the full truncation core (`_mat_trunc!`, `_complete_orthonormal_basis`,
`_dmt_bond_truncate!`, `_dmt_window_truncate!`). The structure is a faithful connector-preserving
DMT (White, Zaletel, Mong, Refael, PRB 97, 035127): orthogonalize to the bond; build left/right
identity (trace) environments; rotate into a connector-aligned basis whose leading columns are the
trace/connector singular directions; preserve the rank-1 identity connector **exactly** (subtract,
truncate the trailing block, add back — with the relative-tolerance guard from the prior audit's
numerics-1 fix); fill the remaining bond budget by SVD-truncating the non-connector block. The
rigorous invariant (trace/identity preservation, sector preservation) is ED-oracle-tested
(`test_dmt.jl` 97 assertions; sector/energy oracles in `test_operator_space.jl`). I did not
line-by-line re-derive every contraction, but the algorithm shape and the ED invariants are sound.

**Load-bearing conventions, pinned down (the brief's "4·maxdim" is outdated):**
- `gate_maxdim` default is `max(maxdim*16, 64)` (`types.jl:127`, `dmt.jl:35,525`; DMTOptions
  reconciled to match, prior convention-3), **not** `4·maxdim`. Its role is the temporary bond
  budget that lets the raw gate apply *exactly* (`cutoff=0`) before DMT truncates back to `maxdim`.
  A single 2-site gate on dimension-4 operator-space sites can grow a bond by at most ×4 (= d), so
  `4·maxdim` is the tight single-gate bound; the `16×` (= d²) default is conservative headroom for
  the 3-site PXP gates (which span two bonds) and successive growth. It is an accuracy↔cost knob
  (the per-gate product/SVD at the inflated bond dominates cost — see the transport review), not a
  theorem-derived constant, and the docstring correctly says to co-tune it with `maxdim`.
- `connector_buffer` default `8` is the number of protected connector operator-directions on the
  boundary *beyond* the special rank-1 trace. For Pauli space a one-site boundary needs ≥4
  (I,X,Y,Z); `8` gives headroom for a broader boundary. Heuristic budget, co-tuned with `maxdim`
  (`connector_buffer ≤ maxdim` is enforced). The rigorous protection is the rank-1 trace
  subtraction, not the buffer size.

perf-1 env cache: not re-verified here (two prior audits confirmed it bit-for-bit; guarded by
`_DMT_VERIFY_ENVS`). The dependency-coupling risk it carries is DEP2.

### DAOE / FDAOE — `src/operator_space/daoe.jl` (algorithm previously unexamined) — **CORRECT**

Verified the channel construction against the source papers (Rakovszky et al. arXiv:2004.05177;
Kuo et al. arXiv:2311.17148) and against an independent Majorana-weight oracle to machine
precision: DAOE reproduces `e^{-γ·max(ℓ-ℓ*,0)}` with `ℓ` the Hamming weight (number of
non-identity sites, **independent of span** — verified directly: contiguous, gapped, and
edge-separated weight-2 strings all damp identically); FDAOE reproduces `e^{-γ(w-w*)}` with the
fermionic/Majorana weight over every Pauli string up to N=5 (JW Z-tails cost 0, interior I's cost
2, parity tracked past saturation). MPO bond dims are `ℓ*+1` (DAOE) and `w*+2` (FDAOE) as
documented; edge cases (γ=0, cutoff>chain, single-site, large γ, off-diagonal=0) all correct;
the `γ≥0` guard (prior numerics-5) is present. The "projector" name is a documented misnomer
(`P_γ² = P_{2γ} ≠ P_γ`); docstrings consistently say "damping MPO".

**Coverage gap closed:** the shipped tests were 4 hand-picked strings. Added a weight-swept DAOE
oracle (sizes 0..5 × `lstar∈{0,1,2}` × `γ∈{0.3,0.9}`, with gapped/edge-separated strings) and the
trace-preservation regression (D2). `operator-space projector MPOs` 46 → 78.

### Chebyshev / KPM — `src/chebyshev/` (never audited) — **mathematically CORRECT**

Verified against Weiße et al. (RMP 78, 275): the three-term moment recursion matches the exact
ED spectral sum `μ_n = Σ_k w_k T_n(E_k)` to ~1e-15 (entangling H); the Jackson kernel is
bit-for-bit the closed form with `g_0=1`; the reconstruction has the correct `1/(π√(1-x²))`
weight, the single-count `g_0·μ_0`, the factor-of-2 on `n≥1`, and the energy→x Jacobian
`1/halfwidth`; the `H̃=(H-center)/halfwidth` rescaling is consistently inverted in
`SpectralFunction`; the variational fit `_optimize_chebyshev_vector!` converges monotonically in
`maxdim` and respects `|μ_n|≤μ_0`. No correctness defect.

**Coverage gap closed:** the only ED oracle was the `m=0` KPM sum rule, which is blind to every
`n≥1` term (`∫T_n/(π√(1-x²))=0`), so a wrong factor on the `2·Σ` would still pass — confirmed by a
factor-5 corruption that integrates to `μ_0`. Added a `T_m`-projection oracle (`m∈{1,2,5}`)
recovering `μ_m` via the same Chebyshev-Gauss quadrature, which the factor-5 corruption fails.
`exact-diagonalization oracles` 8 → 11. Robustness footguns (unrescaled-H silently blows up
moments; non-Hermitian H has its imaginary moment part silently dropped) are all real but
documented; the doubling trick is not used (a ~2× perf enhancement, not a bug).

### Observables / models — `src/observables/`, `src/models/` — **CORRECT**

`energy_density` (verified vs statevector oracle, above); `bond_entropy`/`entanglement_spectrum`
compute the standard von Neumann entropy of the squared Schmidt values with correct
normalization and non-mutation (tested); `fidelity_distance` edge cases tested. TFIM bond
Hamiltonian field-splitting is validated by the ED oracle (it builds *both* the dense ED H and the
TEBD gates and they match); PXP terms by `test_pxp.jl`.

---

## Could not verify / open questions

- **DMT algorithm — line-by-line vs proof.** I verified the algorithm *structure* matches White's
  connector-preserving DMT and that its rigorous invariants (trace/identity/sector preservation)
  hold under the ED oracles, but I did not re-derive every contraction symbolically. A residual
  possibility is a subtle basis-rotation error that happens to preserve all currently-tested
  invariants; unlikely given the breadth of `test_dmt.jl`, but not a closed proof.
- **ITensorMPS 0.3.x is untested** (DEP2). Only 0.4.1 is installed; `leftlim`/`rightlim`/
  `set_nsite!`/`ProjMPO`/`position!` signatures and behavior on 0.3.44–0.3.45 are unconfirmed
  (the same-day 0.3/0.4 release timing makes drift unlikely but unproven; no tracked Manifest, no
  CI matrix).
- **`_optimize_chebyshev_vector!` and `energy_cutoff!` at large N / high entanglement.** Both are
  correct and convergent in the dense-ED-reachable regime tested (≤10 sites); their accuracy floor
  under heavy truncation, and the quantitative filtering quality of `energy_cutoff!` on
  strongly-entangled Chebyshev vectors (which the code itself documents as only "weakly removed"),
  are numerics questions beyond what dense ED can oracle.
- **Secondary items confirmed as documented limitations, not pursued as defects:** `(√2)^N`
  unnormalized-trace overflow (prior numerics-4; far outside current run sizes); the single-pass
  classical Gram-Schmidt in `_complete_orthonormal_basis` (safe only because the protected block
  is hard-bounded ≤4 columns — fragile if the site dimension ever changes from 4); the spin-1/2
  hard-coding of `pauli_state_from_mpo` (dim 2/4). Each is a real latent constraint worth a code
  comment but not a present bug.

---

## Changes landed in this pass (each its own commit, full suite green after each)

| Commit | Change |
|---|---|
| `d6b760f` | **D1 fix** — DMT rebuild preserves `normalize`; ScarFinder regression test |
| `8253d7b` | **DEP1** — `julia = "1.10"` compat bound |
| `168c77d` | **D2 fix** — corrected daoe.md trace claim; DAOE weight + DAOE/FDAOE trace oracles |
| `d98b213` | Chebyshev `T_m`-projection oracle (factor-of-2 sensitive) |

Held for decision (not applied): **R1** (rework the `nstep` override), **R2** (real-MPS ScarFinder
test), **DEP2** (CI test matrix + ITensorMPS straddle), and the minor `nsweeps`-reserved-key gap.
