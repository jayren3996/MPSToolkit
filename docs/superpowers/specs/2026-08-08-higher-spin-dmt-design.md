# Higher-Spin DMT: Faithful Kernel and Generic Local Dimension

**Date:** 2026-08-08

## Goal

Make DMT usable — correctly and efficiently — for models with local Hilbert space dimension
`d > 2` (spin-1, spin-3/2, and beyond), and in the process bring the truncation kernel in line
with the guarantee the method is defined by.

Two things block this today. The operator-space layer is hard-wired to spin-1/2 (`dim(siteind) == 4`
validation, the `pauli_*` basis/gate/state/measurement helpers, `_spinhalf_span`). And the
truncation kernel does not actually deliver DMT's defining property even at `d = 2`: the
protected block is reinstated and then immediately clipped away by the final repair SVD. The
second problem is invisible at spin-1/2 in the diagnostics the suite currently runs, and gets
structurally worse as `d` grows, so it must be fixed before generalizing rather than after.

Higher-spin DMT is established practice, not speculation — arXiv:2205.02853 runs DMT on SU(3)
spin-1 and `d = 4` chains at `chi` up to 256-512. What does not exist is a written-down onsite
operator basis for `d > 2`: that paper does not state which it used, and no DMT paper does. So
the generic-`d` contract this work defines has published precedent for the physics but has to be
pinned down here for the implementation.

## Scope

In scope:

- a faithful DMT bond truncation, generic in `d`, replacing the current kernel
- a generic operator-space layer (`operator_*`) with `pauli_*` retained as `d = 2` wrappers
- performance work targeting `chi >= 400` and/or `d >= 4`
- tests that pin the preservation guarantee directly, at `d = 2, 3, 4`
- documentation of the new budget semantics and of EDKit.jl interop

Out of scope:

- spin-S model builders (`spin_matrices`, XXZ/anisotropy Hamiltonians). The user supplies dense
  Hamiltonians from EDKit.jl; MPSToolkit only needs to accept them.
- generalizing DAOE and the PXP constraint path beyond spin-1/2 (they stay `d = 2`, guarded by
  explicit errors)
- new example scripts; higher-spin validation lives in the test suite
- infinite systems, periodic boundary local gates (already unsupported)

## Background: what the current kernel does and does not do

DMT (White, Zaletel, Mong, Refael, PRB 97, 035127; arXiv:1707.01506) writes the density operator
as an MPS over a vectorized onsite operator basis and, at each bond, rotates into a basis whose
leading directions are the local-operator directions on the sites adjacent to the cut. Writing
the bond matrix as `M = Q_L^T s Q_R`, it preserves the first `d^2` rows and the first `d^2`
columns of `M` **exactly** and SVD-truncates only the lower-right block, after subtracting the
rank-one trace connector `M_{alpha 0} M_{0 beta} / M_{00}`. The consequence is the guarantee that
defines the method: every observable of diameter `<= 3` survives truncation exactly, so local
conserved densities and their hydrodynamic tails are not damaged by compression. The resulting
bond matrix has rank at most `2 d^2 + chi'`.

The current implementation follows this shape but breaks it in three places.

**1. The protected block is clipped away (dominant defect).** `_mat_trunc!` protects a
`connector_buffer`-sized corner and truncates the trailing block to `chi' = maxdim - connector_buffer`
singular values. The reinstated matrix therefore has rank up to `2 * connector_buffer + chi' =
maxdim + connector_buffer`, and the final repair SVD then clips it back to `maxdim` — discarding
roughly `connector_buffer` directions of the block that was just protected. Within the current
API there is no setting that avoids this: the budget is spent before the protected block is
accounted for.

Measured on a 7-site chain, one bond truncation `chi` 25 -> 12, maximum relative error over all
observables of diameter `<= 3`:

| kernel | error |
|:--|:--|
| current, `connector_buffer = 0` | 0.50 |
| current, `connector_buffer = 4` | 0.27 |
| current, `connector_buffer = 8` | 0.26 |
| faithful, `chi'` sized so nothing protected is clipped | 2e-13 |

Over a 20-step XXZ domain-wall melt at `maxdim = 24`, total `S^z` drifts by `4.1e-8` under the
current kernel and by `1.3e-13` under the faithful one.

**2. `connector_buffer` is the wrong abstraction.** The protected subspace is not a tunable
count of directions; it is structurally the `d^2`-dimensional local-operator subspace on each
side of the bond (or `d^(2n)` for `n` sites per side, preserving diameter `2n + 1`). At `d = 2`
the default `8` happens to exceed `d^2 = 4`, so it over-protects with four arbitrary complement
directions. At `d = 3` it is *below* `d^2 = 9`, so the local-operator subspace is not even
covered.

**3. A latent conjugation error.** The paper's pairing is `M = Q_L^T s Q_R` — a transpose on the
left, because `tr(x_L sigma)` is linear, not sesquilinear, in the left Schmidt operator. The
implementation builds the left basis from the column space of `P_L` rather than its complex
conjugate. For a Hermitian density operator the coefficient tensor is real and the two coincide,
which is why this has never shown up; for a non-Hermitian operator (including the off-label
correlator use in `examples/operator_space/pxp_energy_correlator.jl`) it produces ~50% errors on
the same probes.

Why the existing suite passes: it checks a conserved *total* (`total S^z`, tolerance `1e-2`) and
the trace. Both are dominated by the rank-one trace connector, which does survive the clipping
approximately. A total conserved charge is a weak probe — it stayed at `1e-6` while individual
diameter-`<=3` observables were 26% wrong. The baseline suite on `main` is fully green
(106 passing DMT assertions among them), so this is a gap in what is asserted, not a known
failure.

## Evidence gathered during design

A prototype of the proposed kernel was written and measured before this spec, so the claims below
are observations rather than expectations:

- The faithful rule restores the guarantee: diameter-`<=3` preservation at `2e-13` versus `0.26`
  for the current kernel on the same input.
- It generalizes: the same kernel, unmodified, preserved diameter-`<=3` observables to `<= 2e-13`
  at `d = 2, 3, 4` for both Hermitian and non-Hermitian operators, and ran a spin-1 Heisenberg
  melt with `|dS^z| ~ 1e-10`.
- The conjugation convention is load-bearing exactly where predicted: dropping `conj` on the left
  protected block leaves Hermitian operators at `7.6e-14` and pushes non-Hermitian ones to `0.51`.
- The generalized Gell-Mann basis reduces to `(I, X, Y, Z) / sqrt(2)` **exactly** at `d = 2`
  (element-wise, to `1e-15`), and is orthonormal and Hermitian for `d = 2..5`.
- The randomized complement truncation agrees with the dense-optimal one to overlap
  `0.998-0.999` while preserving the guarantee to the same `1e-13`.
- Component profiling reordered the performance plan; see section 3.

## Design

### 1. Faithful bond truncation, basis-free and generic in `d`

At bond `b`, with the orthogonality center at `b` and `d = sqrt(dim(siteind))`. `S` below is the
bond matrix: diagonal Schmidt values if the bond tensor is factorized by SVD, or the triangular
`R` factor under the cheaper QR variant of section 3 — the rule is identical either way.

```
S            = Diag(Schmidt values)                     chi x chi, diagonal, never materialized
P_L          = conj(L_env * A_b)                        chi x d^2   (conj: transpose pairing)
P_R          = B_{b+1} * R_env                          chi x d^2
Q_L, Q_R     = thin orthonormal bases of P_L, P_R       k_L, k_R <= d^2
q0, r0       = normalized first columns of P_L, P_R     identity/trace directions
C            = (S r0)(q0^dag S) / (q0^dag S r0)         rank-1 connector (skipped if ill-conditioned)
B            = S - C                                    diagonal minus rank-1, kept structured
D            = P_L^perp B P_R^perp                      via low-rank updates, not dense projectors
M'           = S - D + truncate_{chi'}(D)
```

with `chi' = maxdim - 2 d^(2n)` the complement budget (`n` sites protected per side; see
section 2). `M'` is assembled directly in factored form `M' = F G^dag`, as the sum of four
low-rank pieces

```
F G^dag = C                                     rank 1     (trace connector)
        + Q_L (Q_L^dag B)                       rank k_L   (protected rows)
        + (B Q_R) Q_R^dag                       rank k_R   (protected columns)
        - Q_L (Q_L^dag B Q_R) Q_R^dag           (removes the double count)
        + U_D Sigma_D V_D^dag                   rank chi'  (truncated complement)
```

so `rank(M') <= 2 d^(2n) + chi'` by construction, and the refactorization never sees a dense
`chi x chi` matrix. The connector needs no separate budget: `B` annihilates the identity
directions (`B r0 = 0`, `q0^dag B = 0`), so the connector's range already lies inside the
protected block — matching the paper's rank bound `2 d^2 + chi'` for `n = 1`
(arXiv:1707.01506 Sec. III.3 and App. B). The complement budget is therefore

```
chi' = maxdim - 2 d^(2n)
```

and the final factorization never has to discard anything protected. (`k_L` and `k_R` are `d^(2n)`
whenever the thin QR is taken without rank detection, which is the safe choice: a rank-deficient
protected block then yields a few arbitrary extra orthonormal directions, which over-protects
rather than under-protects.) This is the single change
that restores the guarantee. It is also the published convention: arXiv:1902.01859 writes
`chi = chi_preserve + chi_extra` with `chi_preserve = 2 d^(2n)`, and the one public Julia/ITensor
implementation (`stuartthomas25/ITensorsDMT`) sets `cut = maxdim - rank_offset` with
`rank_offset = 2 d^2`, erroring when the budget is below it.

Correctness properties this design relies on, both verified numerically during design:

- Truncating a bond preserves *every* observable of diameter `<= 3` exactly, at any bond, so
  composing truncations across a sweep preserves them too.
- The guarantee is independent of how well `truncate_{chi'}(D)` approximates `D`, because `D`
  lives entirely in the doubly-orthogonal complement and the protected border and connector are
  added back exactly. Any approximation whose column space lies in `col(D)` is safe — which is
  what licenses the randomized methods in section 3.

### 2. Budget semantics and options

`maxdim` becomes the **total** post-truncation bond dimension, inclusive of the protected block —
the `chi = chi_preserve + chi_extra` split of arXiv:1902.01859. Validated once at the entry to
`dmt_step!` / `dmt_evolve!` (where `d` is known from the sites):

```
maxdim >= 2 * d^(2n) + 1,      n = (preserve_diameter - 1) / 2
```

| `preserve_diameter` | `d = 2` | `d = 3` | `d = 4` |
|:--|--:|--:|--:|
| 3 (`n = 1`) | 9 | 19 | 33 |
| 5 (`n = 2`) | 33 | 163 | 513 |

Failing this raises an `ArgumentError` naming the local dimension, the diameter, and the implied
floor, at the start of evolution rather than mid-sweep. Since a truncation only fires when
`chi > maxdim`, and `maxdim > 2 d^(2n) >= d^(2n)`, the protected block always fits and the
per-bond budget can never be over-subscribed at run time — which also disposes of the
boundary-bond case that the public implementations have to guard explicitly.

`preserve_diameter::Int = 3` (odd; `n = (L-1)/2` sites protected per side) replaces
`connector_buffer`, which is removed. Passing `connector_buffer` raises an `ArgumentError`
explaining the change and the new floor, rather than a bare `MethodError`. The name and the
odd-diameter convention follow the literature (arXiv:1902.01859 SM runs `L = 3, 5, 7`;
arXiv:2311.17148 App. D reports `L = 3`), so a user reading a paper can set the knob directly.

This is a first-class parameter, not a nicety: arXiv:2311.17148 attributes a systematic
underestimate of its energy-current correlator to running at `L = 3` when the current is a 4-site
operator, and expects `L >= 4` to converge faster. arXiv:1902.01859 further reports that varying
`L` at *fixed* total `chi` changes accuracy, despite preserving the same amount of information —
so the split between `chi_preserve` and `chi_extra` is itself a tuning axis. The `d^(2n)` scaling
means this is much more expensive at higher spin (`L = 5` needs `chi >= 163` at `d = 3`), which
the documentation must state.

`gate_maxdim` keeps its meaning (temporary budget during raw gate application) but its default
is retied to the DMT budget rather than the current `16 * maxdim`; see section 3.

### 3. Performance

Measured baseline (8-site synthetic chain, 32 BLAS threads), per bond:

| regime | gate application | DMT truncation |
|:--|:--|:--|
| `d^2 = 9`, `chi = 200` (inflates to 1800) | 3.4 s | 9.0 s |
| `d^2 = 16`, `chi = 100` (inflates to 1600) | 2.3 s | 6.7 s |
| `d^2 = 9`, `chi = 400` | 15.8 s | — |

Component costs at `chi_inf = 1800`: dense `chi^3` SVD 2.06 s; `svd(chi x d^2; full = true)`
0.505 s; thin QR of the same matrix 0.0013 s.

A prototype of the faithful kernel was profiled per component. At `chi = 1600, d = 3,
maxdim = 200`, per bond:

| component | dense | optimized |
|:--|--:|--:|
| `orthogonalize!` | 0.87 s | 0.87 s |
| factorize `psi[bond]` (SVD) | 4.42 s | **1.06 s** (QR) |
| absorb `v` into `psi[bond+1]` | 0.25 s | 0.25 s |
| protected bases | 0.47 s | 0.47 s |
| truncate the complement `D` | 2.23 s | **0.31 s** (randomized) |
| refactorize `M'` | 0.10 s | 0.10 s |
| **total** | **8.34 s** | **~3.1 s** |

The profile overturned the priority I assumed from the micro-benchmarks: after the obvious wins,
the **canonical-form factorization of the bond tensor is 86% of the cost**, not the DMT algebra.
The plan is ordered accordingly.

**Stage A — structural, no approximation.**

- Build protected bases with a thin QR instead of `svd(...; full = true)`. Two calls per bond:
  0.505 s -> 0.0013 s each at `chi = 1800`. This also removes the orthonormality-defect failure
  mode that the current code documents at length, since no explicit unitary completion is formed.
- Keep `S` diagonal and `B = S - C` structured; form `D` by low-rank updates. Eliminates the two
  dense `chi^3` basis-change matmuls and the two `chi x chi` basis matrices.
- Assemble `M'` in factored form and refactorize via QR of the two `chi x r` factors plus an
  `r x r` SVD, `r <= maxdim`. In isolation this is 14.2 s -> 0.60 s at `chi = 3600, r = 400`;
  in situ it is a small share of the bond step, which is exactly why the profile matters.
- **Replace the bond SVD with a QR.** The kernel needs orthonormal left and right bases and the
  bond matrix expressed in them; it does **not** need the Schmidt (diagonal) form. Writing
  `psi[b] = Q R` instead of `U S V^dag` leaves the truncated operator unchanged — the two
  representations differ by unitaries on both sides, and every step of the rule (protected
  subspaces, projectors, truncation of the doubly-projected complement) is covariant under them,
  while the final refactorization recovers the Schmidt form anyway. Measured 4.42 s -> 1.06 s at
  `chi = 1600`. The matrix-free `B` becomes triangular rather than diagonal, so complement
  matvecs cost `O(chi^2)` instead of `O(chi)`; the net is still strongly favorable.

**Stage B — rank-targeted factorizations.**

- `truncate_{chi'}(D)` by a block randomized range finder with power iterations, using matvecs
  against the structured `D` instead of forming it densely. Measured 2.23 s -> 0.31 s at
  `chi = 1600`; overlap with the dense-optimal result 0.998-0.999 on `d = 3, 4` test states.
- Fuse gate application with the DMT step so the two-site block is factorized **once** at the
  target rank, instead of the current sequence (gate SVDs the two-site block at `gate_maxdim`,
  then the kernel factorizes `psi[bond]` again). This targets the 15-19 s/bond gate cost at
  `chi = 400, d = 3`, which is the largest single remaining item; the size of the win must be
  measured rather than assumed.
- With the fused path, the default `gate_maxdim` becomes "exact gate application, rank-targeted
  factorization" rather than a `16x` inflation followed by a plain-SVD pre-truncation. The plain
  pre-truncation is itself a correctness hazard: it discards small singular values *before* DMT
  sees them, which is precisely what DMT exists to avoid.

Stage B is gated on Stage A landing with the property tests green, and each Stage B change is
required to leave the preservation tests at the same tolerance. Both randomized steps are safe by
the argument in section 1: they cannot break the guarantee, only the optimality of the discarded
weight.

### 4. Generic operator-space layer

Local basis: generalized Gell-Mann, Hermitian, orthonormal under the Hilbert-Schmidt inner
product, identity first as `I / sqrt(d)`. At `d = 2` this reduces exactly to `(I, X, Y, Z) / sqrt(2)`
in that order, so the `pauli_*` wrappers are bit-compatible with today's behavior. Basis matrices
are built once per `d` and cached.

The kernel depends on exactly two properties of the basis — Hilbert-Schmidt orthonormality, and
the identity sitting at index 1 — so the basis is supplied through a single provider function
rather than being hard-wired. This matters because no published DMT paper states which onsite
basis it used at `d > 2` (arXiv:2205.02853 ran `d = 3` and `d = 4` without saying), while the one
public Julia implementation deliberately uses a *non-Hermitian, charge-definite* band basis so
that ITensor U(1) quantum numbers work and the tensors become block-sparse. Hermiticity is not
required by the algorithm. Block-sparse QN-conserving operator space would likely be a larger
performance win at `chi >= 400` than any of the dense linear-algebra work in section 3, but it
changes the tensor type throughout the package and is deliberately **out of scope here**; keeping
the basis pluggable is what makes it a later, localized change rather than a rewrite.

```julia
operator_siteinds(N; d = 2)                      # dim d^2 sites; d recovered as sqrt(dim) elsewhere
operator_basis_state(sites, labels)              # by integer basis label
operator_product_state(sites, ops)               # tensor product of dense local matrices
operator_local_sum_state(sites, O, coeffs)       # sum_j c_j O_j, bond dimension 2
operator_gate(U; d), operator_gate_from_hamiltonian(h, dt; d)
operator_lindblad_generator, operator_gate_from_lindbladian, operator_gate_from_imaginary_time
operator_trace, operator_expectation, operator_expectation_profile
operator_state_from_mpo, operator_superoperator_mpo
operator_gibbs_state
```

Every `pauli_*` name survives as a one-line `d = 2` wrapper, so existing scripts, notebooks and
examples keep working. `pauli_total_sz_state` and `pauli_domain_wall_state` become thin cases of
`operator_local_sum_state`.

`operator_local_sum_state` is load-bearing for higher spin: at `d = 3`, `S^z = diag(1, 0, -1)` is
**not** proportional to any single Gell-Mann element, so an `S^z` domain wall cannot be written
with basis labels and must be built from a dense local matrix.

Interop with EDKit.jl: `spin(...; D = d)` returns a sparse `ComplexF64` matrix and composes
multi-site strings with `kron` left-to-right, first character = leftmost site — the same
convention MPSToolkit uses. Gate builders accept `AbstractMatrix` and densify internally, so

```julia
h = spin((1, "xx"), (1, "yy"), (1, "zz"); D = 3)      # EDKit, sparse, no Array() needed
g = operator_gate_from_hamiltonian(h, dt; d = 3)
```

works directly. This is documented in the DMT and operator-space manual pages.

Trace conventions generalize `2^(N/2)` to `d^(N/2)`; the `Float64` overflow bound becomes
`N ~ 2048 / log2(d)` and is documented alongside the existing note.

### 5. Internal structure

`src/operator_space/dmt.jl` (688 lines and growing) splits into focused units:

- `dmt/types.jl` — `DMTOptions`, `DMTGateEvolution` validation, budget floor checks
- `dmt/environments.jl` — identity/trace environment cache (existing `_DMTEnvCache`, generalized
  and renamed; its amortization argument and verification flag are unchanged)
- `dmt/bond.jl` — the faithful bond truncation
- `dmt/lowrank.jl` — factored refactorization and the randomized range finder
- `dmt/driver.jl` — `dmt_step!`, `dmt_evolve!`, schedules, reverse-sweep gate mapping

New: `src/bases/operator_basis.jl` (generalized Gell-Mann, cached), `src/operator_space/gates.jl`
(gate builders split out of `helpers.jl`).

## Testing strategy

The suite's current weakness is that it never tests the property DMT exists to provide. The new
tests do so directly.

1. **Preservation property (the decisive test).** For `d in (2, 3, 4)`, and for both real
   (Hermitian) and complex (non-Hermitian) coefficient tensors: build an operator-space MPS with
   a genuine bond dimension, truncate one bond, and assert that *every* observable of diameter
   `<= 3` — all one-, two-, and three-site dense probes across the chain — is preserved to
   `1e-12` relative, and that the bond respects `maxdim`. The complex case is what pins the
   conjugation convention.
2. **Sweep-level conservation.** Full forward/reverse DMT evolution of a domain-wall melt at
   `d = 2` and `d = 3`: total charge and total energy conserved to `~1e-12` over many steps, with
   the bond dimension pinned at the budget.
3. **ED oracle.** `N = 6`, `d = 3`: exact Liouvillian evolution of the vectorized density
   operator versus DMT with a budget large enough to be exact, and versus a truncated run where
   local densities must still match exactly.
4. **Budget semantics.** `maxdim` below the `2 d^(2n) + 1` floor raises with a message naming `d`
   and the diameter; `connector_buffer` raises with the migration message; `preserve_diameter = 5`
   preserves diameter-`<=5` observables (and diameter-6 ones are then demonstrably *not*
   preserved, pinning the boundary rather than just the inclusion).
5. **Purity monotone.** `Z = tr(rho) / sqrt(tr(rho^2))` must satisfy `Z >= 1` under DMT while
   plain Frobenius truncation drives it below 1 (arXiv:1707.01506 Sec. IV.B). A cheap sharp
   smoke test that distinguishes the two truncations directly.
6. **Backward compatibility.** `pauli_*` equals `operator_*(; d = 2)` on every wrapped function;
   the `d = 2` Gell-Mann basis equals `(I, X, Y, Z) / sqrt(2)` exactly.
7. **Stage B equivalence.** The randomized complement truncation satisfies the same preservation
   tolerance as the dense path, and its discarded weight is within a few percent of the dense
   optimum on the same inputs. The fused gate+DMT step matches the unfused one on small cases.
8. **Performance guard.** A benchmark script under `dev/` plus a coarse allocation bound in the
   test suite, to catch reintroduction of `chi x chi` temporaries or `full = true` SVDs.

Existing DMT tests that encode the old semantics (`connector_buffer` behavior, `_mat_trunc!`
corner protection, the `1e-2` charge tolerance) are rewritten against the new contract rather
than deleted; the environment-cache regression tests carry over unchanged.

## Migration and breaking changes

- `connector_buffer` removed from `DMTOptions`, `DMTGateEvolution`, `dmt_step!`. Raises an
  `ArgumentError` with the replacement. Call sites to migrate: `src/evolution/types.jl`,
  `src/operator_space/dmt.jl`, `src/operator_space/constrained.jl`, `src/scarfinder/algorithm.jl`
  (ScarFinder constructs a `DMTGateEvolution`), `examples/dmt/*.jl` (3 files),
  `examples/operator_space/pxp_energy_correlator.jl`, `test/test_{dmt,core,pxp}.jl`, and
  `docs/src/manual/dmt.md`.
- `maxdim` now includes the protected block, so a run at fixed `maxdim` keeps fewer complement
  directions than before — and is exactly conserving. Documented prominently.
- `gate_maxdim` default changes.
- Numerical output of the three `examples/dmt/*.jl` scripts shifts. They are updated and re-run;
  the transport exponents they report are re-measured rather than assumed to be unchanged.
- Docs updated: `docs/src/manual/dmt.md` (budget semantics, `d`-generic contract, EDKit interop,
  removal of the "Pauli sites only" restriction), `docs/src/manual/operator-space.md`,
  `docs/src/manual/transport-methods.md` where it references the knobs.

## Success criteria

1. Every observable of diameter `<= 3` is preserved to `1e-12` under bond truncation, at
   `d = 2, 3, 4`, for Hermitian and non-Hermitian operators.
2. A spin-1 transport run executes end to end from an EDKit-supplied Hamiltonian, conserving its
   charge to `~1e-12`, with the bond dimension held at the budget.
3. The `d = 2` public API is unchanged in name and signature except for the documented
   `connector_buffer` / `maxdim` migration.
4. Measured speedup at `chi >= 400` and `d >= 3` versus the current kernel, reported in the
   benchmark script, with Stage A and Stage B contributions separated.
5. Full test suite green, including the existing spin-1/2 tests carried over.

## Risks

- **Scope of the rewrite.** The kernel change touches the one file every DMT workflow depends on.
  Mitigated by landing Stage A behind the property tests before any Stage B optimization, and by
  keeping the environment cache and driver/scheduling logic intact.
- **Randomized truncation quality.** Cannot break the preservation guarantee (argued above), but
  can degrade the discarded-weight optimality. Mitigated by the equivalence test and by keeping
  the dense path selectable.
- **Re-measured example outputs.** The published transport exponents in the docs were produced
  with the old kernel. They must be re-run, and the new numbers reported honestly even if they
  differ.
- **Higher-spin cost floor.** At `d = 3` a `maxdim = 24` run has only 6 complement directions.
  Users must be told that higher spin needs a substantially larger budget; this is a property of
  the method, not of the implementation, and belongs in the documentation. For calibration,
  arXiv:2205.02853 ran `d = 3` at `chi = 128-256`.
- **Transcription risk.** Eq. (33) of arXiv:1707.01506 as printed is wrong (`B' = V Q*_L B_{j+1}`
  should be `V Q^dag_R B_{j+1}`), with no published erratum. This is precisely the class of error
  the preservation property test exists to catch, and a reason not to implement from the paper's
  equations alone.

## References

- C. D. White, M. Zaletel, R. S. K. Mong, G. Refael, *Quantum dynamics of thermalizing systems*,
  PRB 97, 035127 (2018); arXiv:1707.01506. The algorithm; Eqs. (13)-(19) for the transpose
  convention, Sec. III.3 and App. B for the `2 d^2 + chi'` rank bound, Sec. IV.B for the purity
  monotone.
- B. Ye, F. Machado, C. D. White, R. S. K. Mong, N. Y. Yao, *Emergent hydrodynamics in
  nonequilibrium quantum systems*, PRL 125, 030601 (2020); arXiv:1902.01859. SM: preservation of
  `l`-site operators, `chi = chi_preserve + chi_extra` with `chi_preserve = 2 d^(2n)`, and the
  observation that `l` matters at fixed `chi`.
- B. Ye, F. Machado, J. Kemp, R. B. Hutson, N. Y. Yao, *Universal KPZ dynamics in integrable
  quantum systems*, PRL 129, 230602 (2022); arXiv:2205.02853. The precedent for DMT at `d = 3`
  (SU(3) spin-1, Izergin-Korepin) and `d = 4`, at `chi` up to 256-512.
- E.-J. Kuo, B. Ware, P. Lunts, M. Hafezi, C. D. White, arXiv:2311.17148. App. D: preservation
  diameter 3 is insufficient for a 4-site energy current.
- S. Yi-Thomas, B. Ware, J. D. Sau, C. D. White, arXiv:2310.06886. App. C.2: the traceless /
  Heisenberg variant that skips the connector subtraction.
- Public implementations consulted for convention (all research code, mostly unlicensed):
  `stuartthomas25/ITensorsDMT` (Julia/ITensorMPS; total-`maxdim` convention, `rank_offset = 2 d^2`,
  charge-definite qudit basis), `christopherdavidwhite/DMT.jl` (first author's own; the
  `chi_max = 1` marginal-preservation test), `Jack-Kemp/dmt` (C++; `PresRadius`),
  `sajantanand/tenpy` branch `DMT` (asymmetric `k_local`, QR redundancy pruning).
