# Operator-space transport review - 2026-06-13

## Scope

A read-only review of the operator-space transport work on branch `feat/pxp-energy-transport`, consolidating a multi-lens audit (goal-completion, bug-hunt, XXZ-sufficiency, PXP-sufficiency) with adversarial verification of every finding. Three artifact classes are in scope:

1. **Code under review** — the perf-1 DMT environment-cache optimization (commit 3e8a30b: `src/operator_space/dmt.jl` `_DMTEnvCache`/`_left_env_at!`/`_right_env_at!`/`_invalidate_env!`/`_orthogonalize_env!`/`_DMT_VERIFY_ENVS`) and its regression guard (1674dad: `test/test_dmt.jl`), plus the audit-fix commits since the audit baseline 799cfe6 (basis-completion threshold 683c569; normalize field b8e3155; near-zero-trace guard 69582cd; maxdim validation 6463f30; reverse-schedule hoist fd92182; rename 8bbfe31).
2. **Benchmark scripts** — `examples/operator_space/xxz_transport_regimes.jl`, `pxp_energy_transport.jl`, `pxp_energy_correlator.jl`.
3. **Objectives** — PXP plan Tasks 1-6 (`docs/superpowers/plans/2026-06-12-pxp-energy-transport.md`, all checkbox-marked done) and perf-1.

Verification spanned re-reading the cited source at file:line, re-running the shipped XXZ script and reduced PXP runs under `julia --project=.`, recomputing the cited slopes/drifts from the production CSVs in `/tmp/pxp_runs/`, and checking `git ls-files` / `git grep` for committed evidence. **No repo file was modified.**

Verdict tally: of 17 findings across four lenses, **1 was refuted and dropped** (A1-1, a misreading of a Task 6 doc-layout clause), and 16 were confirmed or downgraded. Exactly **one confirmed bug survives at severity >= low** (A2-1, low); the two perf-1/audit-fix bug probes (A2-3, A2-4) found NO defect.

## Baseline

Full suite `julia --project=. test/runtests.jl` is **GREEN**: 1165 Pass / 1165 Total, zero Fail / Error / Broken. Highlights: operator-space DMT 97/97; Pauli vectorization 608/608; PXP model helpers 34/34; expectations 30/30; imaginary-time thermal 12/12; constrained DMT 8/8; energy-correlator-matches-dense-ED 4/4; normalize=false-preserves-scale 3/3; ED oracles 8/8; projector MPOs 46/46; docstring coverage 81/81; docstring examples 3/3; portable docs 17/17. The only stderr noise is a benign ScarFinder informational warning from an unrelated test set. Foundation is solid.

## Goal completion

**Status: complete.** All seven objectives are genuinely implemented, exported, and tested by dense-ED oracles — not merely checkbox-marked.

- Tasks 1-5: every named function (`pxp_term_*`, `pxp_constraint_mpo`, `pauli_state_from_mpo`, `pauli_superoperator_mpo`, `pauli_pxp_constraint_state`/`_projector`, `pauli_trace`/`_expectation`/`_expectation_profile`, `pauli_gate_from_imaginary_time`, `pauli_gibbs_state`, `constrained_dmt_evolve!`) exists in `src/`, is exported in `src/MPSToolkit.jl`, and has a dense-oracle test in `test/test_pxp.jl`.
- Task 6: the example script, the dmt.md/operator-space.md worked-example sections, and the api.md index all exist. **Finding A1-1 (the one goal finding) was REFUTED**: Task 6's "api.md @docs blocks for all new exports" is satisfied by rendered `@docs` blocks living in the manual pages (`docs/src/manual/operator-space.md:191-197`, `docs/src/manual/dmt.md:492,498-502`), indexed from `api.md:75-101`. That is the package's own self-documented convention (`api.md:10,97` explicitly state the manual pages carry the expanded docstrings), and `test/test_docstrings.jl:8-94` asserts `Docs.hasdoc` on every new export. No deliverable is missing.
- perf-1: the env cache is real (struct + `_DMT_VERIFY_ENVS` toggle threaded through `dmt_evolve!`), and its regression test (`test/test_dmt.jl:332-379`) asserts genuine bit-for-bit equality against the from-scratch rebuild (index identity + norm-of-difference under verify-on, plus `_link_dims` canonical-form equality), which I independently confirmed holds on adversarial schedules the test does not exercise.

The one thing NOT version-controlled is the chi=48 production transport **data** that backs the PXP physics conclusion — but that is benchmark evidence, not a code/test objective, so it is recorded under Benchmarking sufficiency (2b), not as a goal gap.

## Bugs

After adversarial verification, **one** confirmed bug survives at severity >= low.

### A2-1 (low) — `constrained_dmt_evolve!` ignores `evo.normalize`, silently re-normalizing a `normalize=false` trajectory

`src/operator_space/constrained.jl:47` defaults the keyword `normalize::Bool=true` (a literal), then forwards it to `dmt_evolve!` (line 68) and applies `normalize && normalize!(rho)` per checkpoint (line 70). Both sibling drivers default to the field instead — `dmt_evolve!` (`dmt.jl:654`) and `evolve!` (`dmt.jl:704`) use `normalize::Bool=evo.normalize`. A user who builds `DMTGateEvolution(...; normalize=false)` to track the conserved unnormalized trace of a traceless operator (the b8e3155 field) and calls `constrained_dmt_evolve!(rho, evo, projector)` **without** an explicit `normalize=false` kwarg gets the field silently overridden to `true`, re-normalizing to norm 1.0 each checkpoint and destroying the absolute-trace diagnostic.

Reproduced (`julia --project=.`, correlator-protocol traceless operator, `maxdim` below the initial bond so truncation sheds Hilbert-Schmidt norm, `evo.normalize=false` in the field): `dmt_evolve!` no-kwarg honors the field (norm != 1.0); `constrained_dmt_evolve!` no-kwarg ignores it (norm = 1.0); `constrained` with explicit `normalize=false` honors it. The path is **unguarded** — every `normalize=false` test (`test/test_pxp.jl:561,608`) passes the kwarg explicitly *and* builds `evo` with the field at its default `true` (`test_pxp.jl:590-600`); the only field-inheritance test (`test_dmt.jl:99-124`) covers `evolve!`/`dmt_evolve!`, not the constrained driver.

**Severity is low, not medium**, because: (a) it is **latent** — both shipped scripts pass `normalize=false` explicitly at the call site (`pxp_energy_correlator.jl:107`, `pxp_energy_transport.jl:127`), so every shipped caller and test is unaffected; (b) the constrained driver is internally consistent with its **own** docstring, which documents a literal `true` default with no `evo.normalize`-inheritance claim (`constrained.jl:32-33`), unlike `dmt_evolve!`'s docstring (`dmt.jl:641`) — so this is a cross-API least-surprise gap, not a contract violation; (c) the scientific deliverable `M2(t)` (`correlator:112`) is a **ratio** `sum((j-c)^2 p[j])/sum(p)`, invariant to overall rescale, and `test_pxp.jl:613-619` asserts both paths give identical profiles — only the secondary absolute-trace / conservation-drift diagnostic is affected.

**Fix sketch (not applied):** change `constrained.jl:47` to `normalize::Bool=evo.normalize`; update the docstring (lines 32-33) to say the default is taken from `evo.normalize`; add a regression test that builds `evo` with `normalize=false`, calls the constrained driver with no kwarg at a truncating `maxdim`, and asserts `norm(rho) != 1`.

### Non-bugs and disclosed nuances (no action / info)

- **A2-2 (info, was low):** `project_every` is silently inert when `evo.nstep==1` — the loop pattern both scripts use (`min(project_every, 1) == 1` always, `constrained.jl:54,56`). Produces no incorrect output (both scripts hardcode `project_every=1`, the correct value), and the chunking semantics are already documented (`constrained.jl:6-7,21`). API-ergonomics nuance, mis-filed under the bug axis.
- **A2-3 (info):** perf-1 env cache — **NO bug found**. Under `_DMT_VERIFY_ENVS=true` (which throws on any cached-vs-rebuilt env mismatch, `dmt.jl:397-400`), `dmt_evolve!` completed with zero assertion trips on far-jumping, bond-revisiting, span-3, scrambled-reverse, alternate-edge, and `maxdim=2`/`maxdim=1` schedules. The `-2/+2` watermark margin (`_invalidate_env!`) is conservatively correct: a 6->9 orthocenter move has true stale range left {6,7,8}/right {7,8,9}, strictly inside the cleared sets. Cross-chunk safety holds (fresh cache per chunk; `rho[:] = apply(projector, rho)` fully replaces the state).
- **A2-4 (info):** `_complete_orthonormal_basis` threshold fix (683c569) — **NO bug found**. The fix corrects the linear-dependence test from `eps(candidate_norm)` (collapses to ~eps^2) to `ambient_dim*eps(1.0)`; orthonormality holds (worst `||B'B - I||` ~ 1e-12) on near-duplicate, random, and standard-subspace protected blocks. The single-pass classical Gram-Schmidt is safe here only because the protected block is hard-bounded to <=4 columns (`dmt.jl:209`) and SVD-orthonormalized first; the supporting sentence calling CGS unconditionally "safe" is overstated but the offending regime is unreachable.

## Benchmarking sufficiency

### 2a — XXZ three-regime script (`xxz_transport_regimes.jl`): demonstration-only

**What it demonstrates (genuine, valuable).** Operator-space DMT separates the three XXZ infinite-T transport regimes qualitatively from ONE short cheap run. The shipped table reproduces digit-for-digit (nsites=30, maxdim=24, t_max=6.8, window [3,6.5]): Delta=0.5 **z=1.078** (ballistic target 1), Delta=1.0 **z=1.412** (KPZ target 3/2), Delta=2.0 **z=2.007** (diffusive target 2), npts=18 each. Targets and citations are textbook-correct and correctly attributed (Ljubotina-Znidaric-Prosen PRL 122 210602; Yi-Thomas arXiv:2310.06886) — the gap is convergence, not wrong physics (A3-1, A3-6).

**What a publishable conclusion needs (not met):**

- *chi-convergence.* Delta=1 drifts monotonically DOWN with bond dimension at the shipped window: z = 1.490 (chi16) -> 1.443 (chi24) -> 1.355 (chi32). chi=24's near-3/2 agreement is not an asymptote; chi=32 lands just below 4/3 by truncation undershoot. Delta=2 is non-monotone in chi (1.732/2.043/1.877 at chi16/24/32) — z=2 is a crossing of a drifting estimator (A3-2, A3-4).
- *fit-window robustness.* Same fixed trajectory, window swept: Delta=1,chi24 = 1.277/1.443/1.522/1.344/1.540; Delta=2,chi24 = 1.712/2.043/2.354/1.854/2.826 (1.7->2.8 swing). Delta=0.5 is the stable exception (1.04-1.13) (A3-3).
- *a scaling plateau.* The Delta=1 local slope is a monotone ramp through 3/2: z_local(t) at chi24 = 1.116(t2)/1.348(t4)/1.495(t5)/1.506(t6)/1.478(t6.8); the window fit averages the ramp to ~1.41 (A3-5).

**Disclosure credit and the one elevated finding.** The script openly admits "tuned so all three z come out close in a short run" (line 47), labels the targets "the predicted physics, not inferred from the fit" (line 55), and describes the early-ballistic/late-truncation squeeze (lines 28-40); the DMT manual repeats the caveat (dmt.md:263-272). Because these are present, A3-2/A3-3/A3-4/A3-6 **downgrade to low/info**. The one finding that stays at **medium is A3-5**: the XXZ script prints only the single window-fit number and affirmatively calls the window a "scaling plateau" recovered "close to good accuracy" (lines 35,133-135), with **no** per-row local-slope column and **no** crossover caveat — whereas both PXP scripts print a per-row local-slope column and narrate the crossover. That disclosure asymmetry (the slope ramp is hidden where the PXP ramp is shown) is a real consistency gap in a shipped benchmark.

### 2b — PXP energy transport (`pxp_energy_correlator.jl` + `pxp_energy_transport.jl`): sufficient-shortterm

**Short-term superdiffusion holds, and only at chi>=48.** At chi=48 (normalize=false, N=64) the infinite-T energy-correlator M2 local slope is robustly superdiffusive: ~1.63 at t=6 descending to 1.45-1.48 at t=9-10, plateauing **1.46-1.50 over t=8-14** (LSQ 2/z=**1.480** over t[8,14], n=31, z_eff=1.35). That is ~48% above the diffusive value 1.0 and approaches KPZ 4/3 strictly from above. The front stays interior (edgeweight 8.96e-9 at t=8 -> 1.86e-7 at t=14, below edge_tol=1e-4) and sector leakage is ~5e-11. Every cited slope was recomputed from `/tmp/pxp_runs/results_l1_chi48_newdefault.csv` and matches to 4 sig figs (A4-1).

**The gap to z=3/2 is not closed.** Even chi=48 plateaus at ~1.46-1.50 and never reaches 4/3 within reach (still ~1.48, weakly descending, at t=14); the asymptote needs t>=20-30. The chi=48 truncation error bar (the conservation drift, which with normalize=false IS the error bar) grows 1.83% (t8) -> 9.41% (t12) -> 12.24% (t14) — so the t=12-14 tail sits at ~10-12% drift, **outside** the "χ-converged to 1-6%" envelope that `dmt.md:369` advertises. Reaching t>=20-30 needs chi LARGER than 48, not merely longer time (A4-3, downgraded to low: slope-plateau half disclosed, but the disclosure understates its own error bar at the top of the window).

**The shipped default (chi=32) reproduces the trap, not the signal.** Both scripts default `maxdim=32` (`pxp_energy_correlator.jl:56`, `pxp_energy_transport.jl:56`). In the matched window t[8,10] the slope is **1.301 at chi=32 vs 1.468 at chi=48** — it moves UP +0.167 with chi (not converged), and chi=32 sits just BELOW 4/3, so a reader fitting the shipped demo reads "z~1.54, KPZ-like" by cancellation of truncation undershoot. The M2 value looks converged (0.04%/0.45% deficit at t4/t8) yet the deficit grows monotonically to 3.86% by t10 — lowering the log-log slope while leaving the value near-converged. Two chi=32 runs with projector_maxdim 32 vs 64 give identical slopes (1.301 vs 1.305), confirming genuine chi=32-DMT runs (A4-2, downgraded to medium).

**Reproducibility status — the binding limitation of the checked-in evidence.** The persuasive chi=48 numbers ("1.46-1.49 over t=8-14") exist ONLY as prose — `git grep` locates them at exactly `dmt.md:369`, `pxp_energy_correlator.jl:39`, and `:147`, nowhere else. There is **no committed dataset, figure, or regression test** pinning the chi=48 plateau (`git ls-files` shows only two logo SVGs, zero CSV; `grep test/` for chi=48/the slope returns nothing). The raw production data lives only in volatile `/tmp/pxp_runs/` (dated 2026-06-12). The N=6 ED oracle validates the correlator MACHINERY bit-for-bit but not the large-N exponent. So an external reader running the repo as-shipped gets only the chi=32 undershoot output and cannot regenerate the headline exponent (A4-5, downgraded to low — reproducibility/hygiene gap on a demo script, additive fix, no code defect).

**Supporting context.** The domain-wall protocol is a weaker probe at reachable times (|dE| slope still ~1.09 at t=6.6, above even diffusive 0.5 — ballistic transient); it is correctly framed as an infrastructure/conservation demo with the correlator carrying the exponent claim (A4-4, info — disclosed). Independent reduced N=32, t<=6 runs reproduced clean mechanics and chi-independence at short times, confirming the code path but not reaching the t>=8/chi>=48 regime (A4-6, info).

**Net.** Short-term superdiffusion (slope ~1.48 >> 1.0, descending toward 4/3) is defensible and chi-converged **at chi>=48 over t<=~12**. z=3/2 is **not established** — approached from above, needs t>=20-30 at chi>48. The repo's prose framing of the physics is honest and well-caveated; the two genuine gaps are (i) the shipped runnable default undershoots, and (ii) the convincing numbers are version-uncontrolled.

## Next actions

**Code (to trust/land):**
1. [low, A2-1] `constrained.jl:47` -> `normalize::Bool=evo.normalize` (match `dmt_evolve!`/`evolve!`); fix the docstring.
2. [low, A2-1] Add the unguarded-path regression test (`normalize=false` field, no kwarg, truncating maxdim, assert `norm != 1`).
3. [info, A2-2] Document that `project_every` chunks `evo.nstep` within one call (inert at nstep=1); optionally warn when `project_every > evo.nstep`.
4. [info, A2-3] Add one far-jumping schedule to the perf-1 regression test under verify-on (belt-and-suspenders).
5. Re-run the full suite after the A2-1 change; confirm 1165/1165.

**Physics (to support the conclusions):**
1. [blocking for any z claim] Commit a chi>=48 reference artifact: a checked-in M2(t) CSV + fit snippet, or a commented chi=48/t>=20 block plus a coarse "slope stays above 4/3 over t=8-12" regression test. The current default (chi=32) reproduces the trap; the convincing numbers live only in volatile /tmp.
2. Correct `dmt.md:369`: chi=48 is "χ-converged to 1-6%" only at t<=~10; drift reaches ~12% by t=14, and t>=20-30 requires chi>48 (a moving target).
3. For an asymptotic measurement, scale chi WITH t and extend to t>=20-30; keep framing z=3/2 as approached-from-above / not-yet-established.
4. [XXZ, A3-5] Add a per-row local-slope-vs-t column (with literature 2/z lines) to the XXZ script, mirroring the PXP crossover honesty.
5. [XXZ] When quoting any XXZ z, report a chi sequence and window-sweep band with an explicit "not converged" statement, not a single window-selected value.

## Verdict

- **Goal completion: complete.** Seven of seven objectives implemented, exported, and ED-oracle-tested; suite green 1165/1165; the lone goal finding refuted.
- **Code: trustworthy, one low-severity fix.** The perf-1 cache is bit-for-bit sound on every adversarial schedule probed; the basis-completion fix is correct. The single confirmed bug (A2-1, low) is a latent cross-API normalize-default inconsistency on the constrained driver, ratio-invariant for the shipped science, fixable in one line plus a guard test.
- **XXZ: demonstration-only.** Recovers the three regimes qualitatively and reproduces the shipped table exactly, but the exponents are not chi-converged or window-robust (Delta>=1), and the script's "scaling plateau" wording overstates a monotone crossover (A3-5, medium) where the PXP scripts disclose theirs.
- **PXP: sufficient-shortterm.** Short-term superdiffusion (M2 slope ~1.48 >> diffusive 1.0, descending toward KPZ 4/3) holds and is chi-converged at chi>=48 over t<=~12. z=3/2 is not established (asymptote needs t>=20-30 at chi>48). The two real gaps are reproducibility (chi=48 numbers are prose-only over volatile /tmp; no committed dataset/test) and that the shipped chi=32 default reproduces the undershoot trap rather than the signal.

## Resolution (applied 2026-06-13)

All four follow-ups selected from this review were implemented on `feat/pxp-energy-transport`, each test-gated (full suite green, 1175 assertions):

- **A2-1 (bug).** `constrained_dmt_evolve!` now defaults `normalize` to `evo.normalize` (`src/operator_space/constrained.jl`); regression test added in `test/test_pxp.jl`.
- **A4-5 / A4-3 (reproducibility).** The chi=48 reference run is committed as `examples/operator_space/data/pxp_energy_correlator_chi48_N64.csv`, re-fit by `examples/operator_space/pxp_energy_correlator_reference.jl` (LSQ 2/z = 1.480, drift 1.8%->12.2%), and pinned in CI by `test/test_transport_reference.jl`. The `dmt.md` reference paragraph is corrected: chi=48 is 1–6% converged only through t≈10 (~9–12% by t=14), and the 4/3 asymptote needs chi>48, not merely longer time.
- **A3-5 (XXZ honesty).** `examples/operator_space/xxz_transport_regimes.jl` now prints a per-row local-z(t) column, exposing the Delta>=1 crossover ramp the single fitted z had averaged away.

A2-2 (`project_every` inert when `evo.nstep==1`) was left as a documented info-level nuance; no other finding required code action.