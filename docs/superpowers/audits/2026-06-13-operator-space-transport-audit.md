# Operator-Space DMT Transport — Audit Report

**Date:** 2026-06-13
**Scope:** `src/operator_space/` (DMT, DAOE/FDAOE, constrained/PXP, vectorize, thermal, expectations, helpers) plus the transport examples, tests, and manual pages. Light pass on `src/bases/pauli.jl`, `src/models/pxp.jl`.
**Method:** multi-agent workflow — per-file inventory → literature reference collection → six-lens audit (correctness · numerics · performance · conventions · docstrings · tests) → adversarial verification (one skeptic per finding, refute-by-default) → low-risk auto-apply. 65 agents.

## Summary

| | count |
|---|---|
| Raw findings | 40 |
| Confirmed (survived adversarial verify) | 37 |
| Refuted | 3 |
| Low-risk fixes auto-applied | 12 |
| Behavior-affecting proposals (held for decision) | 25 |
| Severity: high / medium / low | 0 / 4 / 33 |

**Verification:** full `julia --project=. test/runtests.jl` passes (0 failures, 0 errors) *after* the fixes below, including the regression follow-up.

---

## Resolution (2026-06-13)

The 12 low-risk fixes plus the behavior-affecting proposals were implemented as a test-gated,
commit-per-phase sequence on `feat/pxp-energy-transport` (each phase: TDD where behavioral →
full-suite gate → commit). Final state: **full suite green**.

| Phase | Implemented | Findings |
|---|---|---|
| P0 | Baseline: the 12 low-risk fixes + this report | applied autofixes |
| P1 | Constraint/scar citations; DAOE-vs-`P_G` clarification | docs-7, docs-8 |
| P2 | O(N²)→O(N) reverse-schedule-scan hoist; imaginary-time Hermiticity guard | perf-2, numerics-6 |
| P3 | `maxdim≥1` validator tightening; `gate_maxdim` default reconciliation | convention-2/correctness-4/numerics-3, convention-3 |
| P4 | **`numerics-1`** relative-tolerance trace guard (silent-amplification footgun); **`test-2`** energy-correlator ED oracle; edge tests | numerics-1, test-2, test-4 |
| P5 | `normalize` field on `DMTGateEvolution`; `pauli_fdaoe_projector` canonical + alias; `nstep` alias | convention-6/5/1, test-1/3 |
| P6 | **Bonus correctness bug fix** (see below); coverage tests | test-5/7/8 + `_complete_orthonormal_basis` fix |
| P7 | `(√2)^N` overflow limit documented | numerics-4 |
| P8 | `perf-1` deferred → `perf-4` applied (dead-kwarg footgun removed + design/findings documented in-code); **`perf-1` later implemented — see follow-up below** | perf-1, perf-4 |

**Bonus bug fix (not in the original 40 findings).** Adding `test-8` (rank-deficient
`_complete_orthonormal_basis` coverage) surfaced a genuine correctness bug: the linear-dependence
threshold was scaled by `eps(candidate_norm)` (the post-projection residual's own magnitude),
collapsing to ~`eps²` for roundoff residuals, so a numerically-dependent standard-basis vector was
accepted as a spurious near-duplicate column — producing a **non-orthonormal** "orthonormal" DMT
basis. Fixed to scale by the candidate's unit norm.

**`perf-1` (O(N²) environment rebuild) — implemented (follow-up, 2026-06-13).** Completed on
`perf/dmt-threaded-environment` and merged into `feat/pxp-energy-transport`. `dmt_evolve!` now
builds one `_DMTEnvCache` and threads it through the sweep; each `_dmt_bond_truncate!` extends a
memoized running environment by only the sites mutated since the previous bond (watermark
invalidation with a one-site-beyond margin for the open link index) instead of rebuilding both
identity environments from scratch. `orthogonalize!` and the entire truncation sequence are
unchanged, so a cached environment is *by construction* the same contraction of the same tensors
the rebuild would produce — the truncation output is **bit-for-bit identical**.

This is why the prototype's open problem (threaded forward leaving redundant near-zero bond dims
that broke the reverse sweep) does not recur: the cache memoizes *only* the environment value and
keeps the per-bond re-gauge, so it never changes the gauge or the canonical-form output. The three
"hard" coupled problems dissolve — gauge consistency is automatic (env value is gauge-agnostic, the
per-bond `orthogonalize!` stays), the opaque gate primitive is handled by invalidating the bounded
re-gauge range it touches (via `leftlim`/`rightlim`), and the non-monotonic / mixed-span /
overlapping PXP schedule needs no special casing because invalidation is per-operation and
footprint-based (every gate's footprint is its O(span) window — verified empirically). Correctness
is guarded by a `_DMT_VERIFY_ENVS[]` toggle that asserts every cached env equals the from-scratch
rebuild (index identity + norm of difference): the **full ED-oracle suite passes with it on** (0
failures, 0 assertion trips), and the threaded vs rebuild final states match to overlap 1.0 with
**identical bond dimensions** on the PXP hard case (N=6/8/10), forward and reverse, across multiple
sweeps.

**Measured.** An isolated environment-sweep benchmark confirms the asymptotic fix: rebuild cost
∝ N² (`rebuild/N²` flat at ~13.7), threaded cost ∝ N (`threaded/N` flat at ~39), env-term speedup
growing 4.4× (N=16) → 21× (N=64) → 90× (N=256). **End-to-end**, in the accuracy-preferred regime
(`maxdim=32`, `gate_maxdim=256`), the per-gate `product`/SVD at the *inflated* bond dominates total
cost at reachable N, so the wall-clock speedup is ~1.0–1.05× (break-even, never a regression) —
the threaded env becomes the end-to-end bottleneck, and thus a wall-clock win, only at very large N
or when `gate_maxdim` is closer to `maxdim`. The change is asymptotically correct and zero-risk; it
pays off increasingly as PXP runs are pushed to the long chains needed for the asymptotic z=3/2.

---

## Applied low-risk fixes (verified)

All changes are docstrings/comments, citations, defensive validation, type annotations, or a behavior-identical allocation hoist. Default behavior is unchanged.

| Finding | File(s) | Change |
|---|---|---|
| `docs-1` / `correctness-1` | `docs/src/manual/dmt.md`, `examples/operator_space/pxp_energy_transport.jl`, `examples/operator_space/pxp_energy_correlator.jl` | Corrected the PXP energy-transport reference's fourth author: **Abanin → Papić** (PRX 13, 011033 (2023); actual authors Ljubotina, Desaules, Serbyn, Papić). 5 locations. |
| `numerics-5` | `src/operator_space/daoe.jl` | Added `gamma >= 0` guards to `pauli_daoe_projector` and `fdaoe_projector` (negative `gamma` silently *amplifies* instead of damping). |
| `numerics-2` | `src/operator_space/dmt.jl` | `evolve!(psi, ::DMTGateEvolution; normalize=true)` now forwards `normalize` to `dmt_evolve!` (was hard-coded `true`, silently blocking traceless-trace tracking through the generic interface). Additive, default-preserving. |
| `perf-3` | `src/operator_space/expectations.jl` | Hoisted the single-site Pauli basis out of the per-term loop in `pauli_expectation_profile`. **(See regression note below.)** |
| `convention-7` | `src/operator_space/helpers.jl` | Tightened `pauli_total_sz_state`'s `coefficient` kwarg to `Union{Nothing,Number}`, matching `pauli_basis_state`. |
| `docs-2` | `src/operator_space/dmt.jl` | Fixed a docstring attached to the wrong function (`_reverse_gate_for_step` block sat above `_is_default_reverse_schedule`). |
| `docs-4` | `src/operator_space/dmt.jl` | Documented the previously-undocumented `orthogonalize` / `left_env` / `right_env` kwargs on `_dmt_bond_truncate!`, including the gauge caveat. |
| `docs-5` | `docs/src/manual/operator-space.md` | Disambiguated the `pauli_gate` index notation: `G[α_out, α_in] = tr[P_{α_out}† U P_{α_in} U†]` (row = output, col = input). |
| `docs-6` | `src/operator_space/thermal.jl` | Clarified the `nsteps` docstring of `pauli_gibbs_state` (per-application `dbeta` vs per-side exponent accounting). |

Correctly **skipped** by the apply pass: `convention-4` (targeted a different file), `docs-3` (superseded by `numerics-2`).

### ⚠ Regression caught and fixed

`perf-3` introduced a Julia scoping bug that the test suite caught (the apply agent verified the file *parsed* but not that tests *ran*). The hoisted comprehension `[… for matrix in values(pauli_matrices())]` was inserted into `pauli_expectation_profile`, which has a function-local results array named `values`. Per Julia scoping, any assignment to `values` makes it local to the whole function, so the comprehension resolved `values(…)` to the not-yet-assigned local instead of `Base.values` → `UndefVarError`, breaking all 14 `pauli_expectation*` tests (and anything downstream of `tr(ρ O)`).

**Fix:** renamed the local array `values → results`, which both fixes the error and removes the latent shadowing hazard. Full suite green after the fix.

---

## Open proposals (behavior-affecting — your call)

Grouped by priority. Several low-severity items are different lenses on the same root issue (noted).

### Tier 1 — highest value for the transport work

- **`numerics-1` (medium) — silent 1/ε amplification for traceless operators.** `pauli_expectation_profile` guards the normalized division with an *exact* `iszero(denominator)` check, then divides. For a **traceless** operator (the energy-correlator protocol in `pxp_energy_correlator.jl`), DMT/SVD truncation leaves an `O(ε·‖·‖)` trace residue — `iszero` reports it as nonzero, the guard passes, and the result is divided by ~1e-13, returning garbage amplified by ~1/ε with **no error**. Reproduced: trace `4e-13` → output `9.99e12`. The DMT kernel `_mat_trunc!` already guards the same failure mode with a *relative* tolerance `size·eps·norm`; the expectation layer should match. **Fix:** replace the exact `iszero` with a norm-scaled relative tolerance that throws a clear `ArgumentError` pointing to `normalize=false`.

- **`perf-1` + `perf-4` (medium) — O(N²)-per-sweep identity-environment rebuild (the dominant hot-loop cost). ✓ RESOLVED** (`perf-4` interim footgun removal, then `perf-1` implemented as a threaded environment cache — see the follow-up in the Resolution section above). `_dmt_bond_truncate!` unconditionally rebuilds *both* identity environments from scratch on every bond (`_left_identity_environment` contracts `1:bond-1`, `_right_identity_environment` contracts `bond+2:N`). The in-tree caller never supplies the cached `left_env`/`right_env`, so the fast-path is dead and the per-sweep cost is `Σ_b[O(b)+O(N-b)] = O(N²·χ²)` where a swept running environment gives `O(N·χ²)`. This grows quadratically as runs are pushed to the longer times needed for the asymptotic `z=3/2`. **Fix:** thread a running environment through the sweep (left env for `:R`, right env for `:L`), extend/contract by one site per step, pass it via the *existing* `left_env`/`right_env` kwargs. Gauge bookkeeping is the catch (the inline comment documents a prior `:L` correctness bug) — validate bit-for-bit against the ED energy-conservation / sector-preservation oracles in `test_operator_space.jl` and the connector/identity regressions in `test_dmt.jl`. (`perf-4`: if sweep-threading is deferred, instead delete the dead kwargs to remove the stale-environment footgun.)

- **`perf-2` (low, but a clean safe win) — O(N²) schedule scan per sweep.** `_is_default_reverse_schedule` (and the forward-index map) is recomputed inside `_reverse_gate_index` for every `(bond, index)` on every reverse step, though it is loop-invariant. **Fix:** compute once at the top of `dmt_evolve!` and thread the flag/map in. Behavior-preserving.

- **`test-2` (medium) — energy-correlator protocol is never validated against ED.** `pxp_energy_correlator.jl` runs the load-bearing protocol — `O(0)=P_G h_c P_G` (traceless), `constrained_dmt_evolve!(…; normalize=false)`, measure `C(x,t)=tr(P_G h_x O(t))` — but the tests cover only the *pieces* (the `normalize=false` test uses a large-trace Gibbs ρ, which cannot catch a sign/normalization bug specific to the traceless single-density correlator). **Fix:** add an N=6 testset at exact bond dimension asserting the full `normalize=false` profile matches the dense `C(x,t)=tr(P_G h_x U O0 U†)` site-by-site, plus time-invariance of the conserved `Σ_x C(x,t)`.

### Tier 2 — correctness/robustness hardening

- **maxdim==0 inconsistency** (`convention-2` = `correctness-4` = `numerics-3`, one fix): public constructors require `maxdim≥1`, but `_validate_dmt_step` accepts `maxdim≥0` and `_dmt_bond_truncate!` silently no-ops on `≤0`, so the keyword `dmt_step!` path silently skips all truncation. **Fix:** tighten `_validate_dmt_step` to `maxdim≥1` and drop the `maxdim==0 ||` branch in the connector_buffer guard.
- **`normalize` reachability through the uniform interface** (`convention-6`, with `correctness-2`/`test-1`): `numerics-2` already added the `normalize` keyword to `evolve!`. Remaining optional improvement: carry `normalize` as a `DMTGateEvolution` field (mirroring `TDVPEvolution`) so the choice survives anywhere the struct is passed. Add a regression test pinning the generic-dispatch behavior (`test-1`).
- **`numerics-6` — `pauli_gate_from_imaginary_time`/`pauli_gibbs_state` assume Hermiticity but only check squareness.** Add a tolerant `isapprox(h, h'; atol=…)` guard (or symmetrize) so a non-Hermitian input fails loudly rather than producing a non-positive "Gibbs" slice.
- **`numerics-4` — `(√2)^N` normalization overflows to Inf for N ≳ 2048** in the `normalize=false` `pauli_trace`/profile paths. Either document the `N < ~2000` limit or distribute the `√2` per site during contraction. (Far outside current run sizes — low urgency.)
- *(`correctness-3`/`convention-8` — negative-`gamma` guard — already resolved by the applied `numerics-5` fix.)*

### Tier 3 — conventions, docs, test coverage

- `convention-1`: `nstep` (DMT/LocalGate) vs `nsteps` (Gibbs/TDVP) keyword spelling diverges. Standardize or document; `pauli_gibbs_state`/`TDVPEvolution` are public, so alias rather than hard-rename.
- `convention-3`: `gate_maxdim` default is duplicated as literal `480` (DMTOptions) vs `max(maxdim*16, 64)` (DMTGateEvolution/dmt_step!) — reconcile.
- `convention-5`: `fdaoe_projector` is the only public Pauli-space constructor lacking the `pauli_` prefix; rename with a deprecated alias for symmetry with `pauli_daoe_projector`.
- `docs-7`: add constraint-origin citations (Fendley hard-square; Turner scars) on the PXP constraint helpers.
- `docs-8`: clarify in `daoe.md` that DAOE/FDAOE projectors are one-sided diagonal damping MPOs (they change the trace), unlike the Hermitian `P_G` sandwich — so they are not drop-in `P_G` replacements in `constrained_dmt_evolve!`.
- Test coverage (`test-1`, `test-3`, `test-4`, `test-5`, `test-7`, `test-8`): generic `evolve!(::DMTGateEvolution)` dispatch; DAOE/FDAOE argument validation; `pauli_expectation_profile` error/edge paths (empty terms, zero-trace guard, out-of-range start); `constrained_dmt_evolve!` projector kwargs and `project_every` invariance; non-Pauli-site validation throw; `_complete_orthonormal_basis` "could not complete" / rank-deficient paths.

---

## Refuted findings (transparency)

The adversarial pass rejected 3 plausible-but-wrong findings:

1. **`correctness-5` — "default coefficient `2.0^(N/2-1)` overflows."** Arithmetic was wrong (`2.0^512 ≈ 1.3e154`, finite; overflow only at N≈2050, ~2× the claimed threshold, and unreachable — `2.0^k` is exact up to k=1023). Not actionable.
2. **`convention-9` — "PXP constructors validate dim-4 inconsistently."** Premise backwards: the PXP constructors delegate to a clear dim-4 error one frame down; the genuinely under-guarded cases are the DAOE projectors (out of scope). Cosmetic with a false rationale.
3. **`test-6` — "DMT truncation should preserve Hermiticity; add a test."** Empirically refuted: the skeptic *ran Julia* and showed DMT truncation does **not** preserve Hermiticity (max imag part `0.0 → 0.74` across a sweep) — `_mat_trunc!` operates on a non-symmetric reduced matrix by design. The proposed test would assert a false invariant.

---

## Literature reference map (collect-only)

Canonical primary sources, mapped to the code they underpin. Useful for adding citations (most docstrings/manual pages already cite a subset).

**DMT** — White, Zaletel, Mong, Refael, *Quantum dynamics of thermalizing systems*, PRB 97, 035127 (2018), arXiv:1707.01506 → `dmt.jl` (`_mat_trunc!`, `_complete_orthonormal_basis`, `_dmt_bond_truncate!`, `DMTOptions`).
**DMT benchmark** — Yi-Thomas, Ware, Sau, White, *Comparing numerical methods for hydrodynamics in a 1D lattice spin model*, PRB 110, 134308 (2024), arXiv:2310.06886.
**DAOE** — Rakovszky, von Keyserlingk, Pollmann, *Dissipation-assisted operator evolution…*, PRB 105, 075131 (2022), arXiv:2004.05177 → `daoe.jl` (`pauli_daoe_projector`, `_daoe_transition`).
**FDAOE** — Kuo, Ware, Lunts, Hafezi, White, *Energy diffusion in weakly interacting chains with fermionic DAOE*, PRB 110, 075149 (2024), arXiv:2311.17148 → `daoe.jl` (`fdaoe_projector`, `_fdaoe_*`). (Also benchmarks DMT energy transport.)
**DAOE (open-system)** — Yoo, White, Swingle, *Open-system spin transport and operator weight dissipation*, PRB 107, 115118 (2023), arXiv:2210.06494.
**Vectorized/superoperator MPS** — Prosen, Žnidarič, *Matrix product simulations of NESS of quantum spin chains*, JSTAT (2009) P02035, arXiv:0811.4188 → `vectorize.jl`, `helpers.jl` (`pauli_lindblad_generator`, `pauli_superoperator_mpo`).
**Operator-space entanglement** — Prosen, Pižorn, PRA 76, 032316 (2007), arXiv:0706.2480 → motivates DAOE size-damping & DMT identity-protection.
**PXP energy transport (z=3/2 KPZ)** — Ljubotina, Desaules, Serbyn, **Papić**, *Superdiffusive Energy Transport in Kinetically Constrained Models*, PRX 13, 011033 (2023), **arXiv:2210.01146** → both PXP examples, `constrained_dmt_evolve!`, `pxp.jl`. *(Note: arXiv:2207.05744 is a wrong id for this paper — it resolves to an unrelated graphene-oxide paper.)*
**PXP / scars** — Turner et al., Nat. Phys. 14, 745 (2018), arXiv:1711.03528; PRB 98, 155134 (2018), arXiv:1806.10933. **PXP / Rydberg** — Bernien et al., Nature 551, 579 (2017), arXiv:1707.04344. **PXP constraint origin** — Fendley, Sengupta, Sachdev, PRB 69, 075106 (2004), arXiv:cond-mat/0309438.
**XXZ transport** — Ljubotina, Žnidarič, Prosen, *KPZ Physics in the Quantum Heisenberg Magnet*, PRL 122, 210602 (2019), arXiv:1903.01329; Bulchandani, Gopalakrishnan, Ilievski, *Superdiffusion in spin chains*, JSTAT (2021) 084001, arXiv:2103.01976; De Nardis et al., PRL 125, 070601 (2020), arXiv:2003.13708 → `xxz_transport_regimes.jl`.

---

## Reusing this audit

The workflow script is saved and re-runnable (per-session path printed when launched). Re-run after addressing proposals to confirm no new issues, or re-scope it to other subsystems (`src/evolution/`, `src/chebyshev/`).

**Process note:** the apply agents validated that edited files *parsed*, but a Julia variable-shadowing bug only surfaced at *runtime* (`perf-3`). Future audit workflows that auto-apply should run the test suite inside the workflow (or restrict auto-apply to non-executable files — docs/comments/citations — and route executable edits through a test-gated stage).
