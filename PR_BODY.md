Makes DMT usable at local Hilbert space dimension `d > 2`, and fixes the truncation kernel, which did not deliver DMT's defining guarantee even at spin-1/2.

## The defect

DMT (White/Zaletel/Mong/Refael, PRB 97 035127) is defined by preserving *exactly* every observable of diameter `<= 3` around a truncated bond. The old kernel reinstated the protected block and then let the final repair SVD clip it back to `maxdim`, so the protection was discarded. Measured on a single bond truncation, diameter-`<=3` observables came out **26-48% wrong**.

The suite was green throughout, because it only ever asserted a conserved *total* — which the rank-one trace connector keeps small even when the local data is badly wrong. A conserved total is a weak probe: it stayed at `1e-6` while individual observables were 26% off.

| | before | after |
|:--|--:|--:|
| `d = 2` preservation error | 0.34 | **1.0e-14** |
| `d = 3` | *blocked* | **4.3e-13** |
| `d = 4` | *blocked* | **2.2e-12** |
| `preserve_diameter = 5` | *unsupported* | **8.1e-14** |

Verified against an exact-diagonalization oracle at `d = 3`, on random dense probes outside the committed probe set, and at the budget floor itself. The guarantee's *edge* is asserted both ways: windows one site wider break at `1e-2`-`1.6e-1` while covered windows of the same width and same state hold at `1e-15`.

## What changed

**Kernel.** Basis-free truncation sizing the complement budget as `chi' = maxdim - (2 d^(2n) - 1)`, so the protected block fits inside the budget instead of being clipped. `conj` on the left protected block, matching the paper's `M = Q_L^T s Q_R` transpose pairing — silently correct for Hermitian operators, ~50% wrong otherwise.

**API.** `connector_buffer` removed (raises an `ArgumentError` naming its replacement); `preserve_diameter` (odd, default 3) replaces it. `maxdim` is now the **total** bond dimension inclusive of the protected block. Floors: 8 / 18 / 32 at `d = 2 / 3 / 4`.

**Generic operator space.** `operator_*` layer generic in `d`, with every `pauli_*` name preserved as exactly the `d = 2` case. Generalized Gell-Mann onsite basis, reducing exactly to `(I, X, Y, Z)/sqrt(2)` at `d = 2`. Hamiltonians from `EDKit.spin(...; D = d)` pass straight in, sparse included.

**Performance.** QR instead of SVD for the bond factorization — identical truncated operator, since the kernel needs orthonormal bases and the bond matrix in them, not the Schmidt form. That is worth **1.2-1.9x per bond step against the unoptimized new kernel**; measured **against the OLD kernel the faithful one still costs 1.3-1.4x** (28.1 -> 37.2 s at L = 64; 35.8 -> 50.7 s at L = 96), buying 5-7 orders of magnitude on conserved-charge drift. `gate_maxdim = 0` (exact gate) as default. Optional randomized complement truncation, left opt-in because `dmt_evolve!` sweeps `:R` then `:L` and independent sketches would cost run-to-run reproducibility.

**Trace guard.** `operator_expectation_profile`'s `normalize = true` guard compared the identity *coefficient* `tr(rho)/d^(N/2)` against `||rho||_HS`, carrying an implicit `d^(N/2)`, and so rejected legitimate low-temperature thermal states as the chain grew — confirmed throwing on positive bond-dimension-1 product Gibbs states at `beta = 0.5, L = 240` and `beta = 1.0, L = 120`. It also silently passed when `norm(rho)` overflowed to `NaN`. The comparison is now `|tr(rho)| > sqrt(eps) * ||rho||_HS` in physical units, evaluated in log space: since `tr(rho) >= ||rho||_HS` for any positive operator, the guard cannot fire on a physical density matrix, while a traceless operator's cancellation residue still throws. Non-finite scales now raise instead of passing. The defect predates this branch. (Reported by an external reviewer; the `_dmt_connector` sibling guard was checked and does **not** share it — `q0`/`r0` are unit-normalized first, so the environment norms divide out, verified `L`-independent to 9-16 significant figures.)

## Breaking changes

- `connector_buffer` removed. **Old `d >= 3` results have to be re-run, not reinterpreted** — raising `connector_buffer` to cover `d^2` does not rescue the old kernel, because it clipped the reinstated protected block regardless of the buffer size, and the dependence on the buffer is not even monotonic.
- `maxdim` now includes the protected block, so a run at fixed `maxdim` keeps `2 d^(2n)` fewer complement directions than before. This is why the manual's XXZ demo exponents moved; the docs explain it rather than hiding it.
- `gate_maxdim` default changed from `max(maxdim * 16, 64)` to `0` (no cap). At `d = 2` that is a no-op, since a two-site gate inflates a bond by at most `d^2 = 4`. At `d >= 3` ported code that passed a finite cap will now cost more time and gain accuracy, because the old cap actually bound.

## Benchmarks

Two, both honest about what they do and do not show.

`examples/dmt/spin1_melt.jl` — spin-1 melt exponents. ULS/SU(3) gives `z = 1.507` against KPZ 1.500; spin-1 Heisenberg gives 1.514 against the diffusive 2.0. **The two do not separate**, so this run does not discriminate the models and the file says so plainly. The crossover needs `t ~ 20-25`; bond dimension is demonstrably not the bottleneck.

`examples/dmt/spin1_semiexact_validation.jl` — the decisive one. Three arms differing only in truncation (cutoff-only reference, DMT, plain SVD TEBD at equal `maxdim`), with guaranteed triangle-inequality bounds against the reference's own measured error floor.

Its result is not the expected one. **Near its budget floor, DMT is measurably worse than the plain SVD it replaces** — `>= 11x` at the last time all diagnostics are valid, `>= 39x` at the last floor-covered time, robust across cutoffs, wall amplitudes, metrics, models and `d`. What DMT buys unconditionally is *invariance*, not raw accuracy: immunity to the cutoff-induced error floor that traps plain SVD in linear response (`>= 29x` at `d = 3`, `>= 35x` at `d = 2`), and `tr(rho)` held to `~1e-13` where plain SVD drifts to `1.3e-4`.

And no budget-only rule predicts the verdict: at identical `d = 3`, `chi' = 22`, `t = 1.2`, spin-1 Heisenberg is a resolved loss (SVD/DMT = 0.504) while ULS is a resolved win (1.775). Overhead fraction does not predict it either — ULS at 64.3% overhead wins while `d = 2` Heisenberg at 44.4% loses. The practical guidance is therefore to run the equal-`maxdim` check on your own Hamiltonian, which is what that script is for.

Both benchmarks are limited to `t <= 1.8` (the exact reference becomes unaffordable beyond it), so neither speaks to asymptotics. They constrain DMT's claim rather than refuting it.

## Not in scope

DAOE and the PXP constraint path remain spin-1/2 only, guarded by explicit errors. Fusing gate application with the DMT step was specified and not implemented; measurements point to implicit-Householder `Q` (~60% of a bond step) as the larger win.

## Verification

37 commits, 44 files, +9683/-1036.  Full suite green (79 testsets). Every task passed an individual spec-and-quality review with a fix loop, followed by a whole-branch review and its own fix wave. Reviews independently verified the load-bearing claims rather than trusting reports — mutation-testing the `conj`, re-deriving the gate convention against five wrong-convention variants, stress-testing the interval arithmetic against 1.6e7 adversarial draws, and re-running both benchmarks' headline numbers.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01RzGfeGxG81hUzxLrx8GGxP
