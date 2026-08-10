# MPSToolkit: improvements derived from Jack-Kemp/dmt

## Context

`feat/higher-spin-dmt` is merged. DMT is correct and generic in `d` (preservation error ~1e-13 at
`d = 2/3/4`), but two things are unresolved: it is **1.3–1.4× slower than the old (broken) kernel**,
and our headline benchmark says DMT is *worse* than plain SVD at equal `maxdim` — a result we
published in the PR body and the manual.

Reading Jack Kemp's C++/ITensor reference implementation (`github.com/Jack-Kemp/dmt`; Kemp is a
co-first-author on the KPZ paper) settled both. It independently confirms our kernel fix — their
`subMaxDim = maxDim - (sdimL + sdimR - 1)` (`DMT.h:719-720`) is our budget subtraction, and their
`MaxDim` is likewise the total inclusive of the protected block — and it exposes four things we are
leaving on the table. Two are performance, one is a latent defect in our own code, one is a
capability, and one is a misreading of our own benchmark.

**Outcome:** ~2.5× faster DMT at `d ≥ 3`, half the peak memory, `preserve_diameter = 5` reachable at
`d = 4`, survivable multi-hour cluster runs, and an honest answer to whether DMT actually beats
plain SVD.

## The evidence base

Measured, not assumed (a `T`-generic transcription of the bond step was run against the shipped
kernel):

| | result |
|:--|:--|
| complex path vs shipped kernel | rel diff **0.0** — the design reduces to today's code expression-for-expression |
| real vs complex path | 1.78e-15 |
| diameter-3 preservation, real vs complex | 2.53e-15 vs 3.70e-15 (real slightly better conditioned) |
| whole bond step, `d=3, χ=800` | **1.92× / 1.68×** (median / min) |
| gate application, `d=3, χ=200` | **2.30× / 2.80×** |
| composite per DMT sweep | **~1.9×** |
| `operator_gate` imaginary residue, all builders, `d=2/3/4` | 1e-18 … 1.2e-16 |
| Hermitian melt-state MPS tensors | `max|imag|` **exactly 0.0** |

Two facts that reframe the work:

1. **Today's kernel already pays the mixed-eltype penalty.** `bond.jl:275`/`278` contract a real
   `left_isometry` against a `ComplexF64` tensor; `NDTensors/.../dense/tensoralgebra/contract.jl:197-211`
   materializes a promoted copy (184 MB per bond at `χ=800`). A real MPS through the shipped kernel
   is therefore **not faster** (0.86–1.16×) and comes back as a **mixed `[Float64, ComplexF64]` MPS**.
   The `bond.jl:111-113` comment states the all-or-nothing constraint and the code violates it.
2. **`operator_gate` is real by theorem, for arbitrary `A`.** `conj(G[a,b]) = tr((P_a A P_b A†)†) =
   tr(A P_b A† P_a) = G[a,b]` because every basis element is Hermitian. Hermiticity of `A` is not
   required. The downcast recovers an exact property rather than approximating one.

**Benchmarking caveat:** this box ran loadavg 121–136 on 64 cores with four other `julia` processes
during measurement. Speed numbers must be re-taken on a quiet box; accuracy numbers are unaffected.

---

## Phase 0 — Benchmark hygiene and the saturated-χ question

Runs in the background while Phase 1 proceeds. No `src/` changes.

**0.1 `dev/bench_dmt.jl`: add an eltype axis.** Parameterize `bench_state` (`:48`), `bench_gate`
(`:149`), `bench_operator_state` (`:161`) and `bench_sweep`'s `rho` (`:203`) on `elt`, sweep
`elt in (Float64, ComplexF64)`. Lower `BLAS.set_num_threads` (`:34`) to ~8 and print
`/proc/loadavg` in the header — the current 32 threads oversubscribe and corrupt small cells. Add a
peak-memory column for `:dense`. Capture the baseline before any Phase 1 edit.

**0.2 Saturated-χ arm in `examples/dmt/spin1_semiexact_validation.jl`.** Our published negative
result was measured at `chi' = 22` against a protected block of `2d² = 18` — the protected block was
**45% of the budget**. Kemp's default is `Cutoff = 1e-16, AbsoluteCutoff = true`, documented as
"the DMT algorithm will immediately saturate to MaxDim": their designed regime is complement budget
≫ protected block, which we never tested.

Sizing matters, and both ends bind. The uncapped reference reaches `χ = 337` at `t = 1.2` and
`χ = 873` at `t = 1.6`, so `maxdim` must sit well above `18` (else we are back at the floor) **and**
well below the exact rank (else no truncation occurs and every arm equals the reference, giving a
vacuous ratio of 1). Use `maxdim ∈ {60, 100, 160}` — overhead fraction 30% / 18% / 11% versus
today's 45% — and **keep the existing floor cells**, because the scientific claim is
floor-versus-saturated, not a replacement number. Report overhead fraction per cell.

**0.3 Act on the result.** If DMT wins at low overhead fraction and loses only near the floor, both
results stand and the guidance becomes sharp: *budget the protected block as overhead, not as your
budget.* Update `PR_BODY.md`'s Benchmarks section, `docs/src/manual/dmt.md`, and the header comment
of the validation script. If DMT still loses at 11% overhead, that is a much stronger negative
result and must be said plainly in the same places.

---

## Phase 1 — Real (Float64) arithmetic path

The core. Two independently shippable sub-phases; **kernel first**, because Phase 1A is a no-op for
production (production states are still complex → `T = ComplexF64` → identical code) while Phase 1B
alone would *hurt* (real state through today's kernel = 0.86×).

Design rule throughout: **`T` is a promotion over the whole chain, never a narrowing.** Every
`Matrix{T}(x)` inside the kernel widens, so silent loss of an imaginary component is structurally
impossible there. Narrowing happens at exactly three audited boundaries, each guarded.

### 1A — Kernel (no user-visible change)

**1A.1 Failing tests first.** `test/dmt_test_helpers.jl`: add
`mps_eltype(psi) = mapreduce(eltype, promote_type, ITensorMPS.data(psi))`. Assert
`mps_eltype(rho) === Float64` after `_dmt_bond_truncate!` in `test/test_dmt_kernel.jl`, and inside
the existing `"Hermitian (real coefficients)"` arms of `test/test_dmt_preservation.jl:17` and
`test/test_dmt_higher_spin.jl:69` — **before and after** the truncation. These fail today. Without
the "after" assertion those arms keep silently exercising the complex path, which is exactly the
present situation.

**1A.2 `_mps_eltype`** in `src/operator_space/dmt/bond.jl`:

```julia
function _mps_eltype(psi::MPS)
  T = eltype(psi[1])
  for n in 2:length(psi); T = promote_type(T, eltype(psi[n])); end
  T <: Number || return ComplexF64   # NDTensors.EmptyNumber: an unset tensor
  return float(T)                     # Int-valued tensors -> Float64
end
```

Call it right after the short-circuit at `bond.jl:191`. It must scan the **whole chain**, not just
`psi[bond]`: the step also consumes `psi[bond+1]`, the interior protected sites (`bond.jl:228-235`)
and both identity environments (`bond.jl:208-212`), which contract over every tensor. Cost is ~100 ns
against `O(χ³)`. Both the `EmptyNumber` and `float` branches need a unit test.

**1A.3 `_dmt_bond_factorize`** (`bond.jl:104`) grows `::Type{T}=ComplexF64` as a fourth positional
argument; `bond.jl:114` becomes `convert(Matrix{T}, ...)`. The `ComplexF64` default keeps
`dev/bench_dmt.jl:62` working unchanged. Rewrite the `bond.jl:111-113` comment: its premise stops
being true, and the reason the conversion must be all-or-nothing is now the NDTensors promoted-copy
behaviour.

> **Trap:** `bond.jl:123` stays `Diagonal(real.(diag(...)))` and must **not** participate in
> inferring `T`. Singular values are real and `_bond_mul(::Diagonal, x)` is elementwise, never BLAS,
> so a real `Diagonal` against complex `x` is free. Inferring `T` from `eltype(bond_matrix)` instead
> makes `_dmt_bond_truncate!(complex_psi, b; factorize=:svd)` pick `T=Float64` and throw on
> `protected_right`. Add a regression test for that exact call.

**1A.4 `lowrank.jl` helpers — defaulted `T`, so `test/test_dmt_kernel.jl` needs zero edits.**

| helper | change |
|:--|:--|
| `_protected_basis` (`:13`) | `(protected, ::Type{T}=eltype(protected))` |
| `_unit_direction` (`bond.jl:58`) | drop the `ComplexF64.` coercion: `x ./ s`, zero case `zeros(float(eltype(x)), …)` |
| `_dmt_connector` (`:46`) | `T` defaulted to `promote_type` of its three arguments |
| `_dmt_complement_ops` (`:70`) | unchanged — already generic, closures inherit from captures |
| `_truncated_svd` (`:113`) | `(mul, adj, chi, rank, ::Type{T}=ComplexF64; mode, …)` — the only signature that *must* grow a parameter, since it cannot infer through closures |
| `_dmt_refactor` (`:144`) | unchanged; `promote_type(eltype(left), eltype(right))` internally |
| `_bond_mul` / `_bond_adj` (`:28-31`) | unchanged |

**1A.5 Function barrier.** Extract `bond.jl:254-269` into
`_dmt_bond_solve(::Type{T}, bond_matrix, protected_left, protected_right, d, radius, maxdim, cutoff, truncation)`
returning `(new_u, new_s, new_v)`. Contains zero ITensors → two clean specializations, and the DMT
stage becomes unit-testable in isolation. Move the `_dmt_complement_ops` call **inside** it so its
closures capture concretely-typed values — removes a per-`mul` dynamic dispatch the current code pays.

**1A.6 `_DMT_VERIFY_ELTYPE = Ref(false)`** beside `_DMT_VERIFY_ENVS` (`dmt.jl:160`), asserting
`eltype(factor_left) === eltype(factor_right) === T` at `bond.jl:267-268`. **Do not skip this.**
`hcat(Matrix{Float64}, Matrix{ComplexF64}) === Matrix{ComplexF64}`, so a partial conversion
re-promotes silently at exactly that line and you pay the complex price believing you are real. This
is the only cheap detector. Enable it in the oracle tests, as `_DMT_VERIFY_ENVS` already is.

`bond.jl:275-279` then needs no change: `left_isometry` (eltype ≤ `T`) against a `T` tensor is
homogeneous once `T` comes from the chain — which is what removes the 184 MB promoted copy.

### 1B — Boundaries (production flips to real)

**1B.1 `_real_operator_matrix`** — downcast when the imaginary part is noise, using the repo's own
threshold (`sqrt(eps(Float64))` relative, as at `gates.jl:119`, `lowrank.jl:49`,
`expectations.jl:191`, `dmt.jl:231`):

```julia
norm(imag(gate)) <= sqrt(eps(Float64)) * max(norm(gate), one(Float64)) || return gate
return real(gate)
```

Applied at exactly **two** returns: `operator_gate` (`gates.jl:75`) and
`operator_lindblad_generator` (`gates.jl:165`). Every other builder routes through one of them and
inherits it free (`exp(::Matrix{Float64})::Matrix{Float64}`). Physical-space intermediates
(`gates.jl:66,71,97,118,147,152`) stay `ComplexF64` — `exp(-im dt h)` is genuinely complex. **Rule:
coerce physical-space intermediates as now; downcast only the operator-space output.**

Return the complex matrix rather than throwing when the residue is large: `operator_basis_matrices`
explicitly sanctions a substituted non-Hermitian basis, for which a complex gate is correct. Kemp's
absolute `1e-12` is not scale-invariant and would misfire on a large-norm generator.

Add `const _OPERATOR_GATE_FORCE_COMPLEX = Ref(false)`. Setting it must reproduce today's numbers
**exactly** across the whole suite — that is how 1B proves "the complex path is unchanged"
mechanically instead of by inspection.

**1B.2 Container re-widening: `src/operator_space/thermal.jl:97`**, `gates = Matrix{ComplexF64}[]`.
Non-obvious and load-bearing — without it the whole `operator_gibbs_state` path stays complex after
the gate downcast. Build the vector from a comprehension so its eltype follows the gates. After 1B,
grep for `Matrix{ComplexF64}[`, `Vector{ComplexF64}(undef`, `ComplexF64[` to catch any sibling.

**1B.3 State builders.** `_operator_coefficients` (`states.jl:77`) is the pivot — measured exactly
`0.0` imaginary for a Hermitian local operator, 0.89–0.98 relative for a non-Hermitian one. Compute
complex, downcast through the same helper, and every builder inherits it. Then one `T` per builder:
`operator_product_state` (`states.jl:135,136`), `_local_sum_state` (`states.jl:179-203`), and
`operator_basis_state` (`states.jl:101-102`) — the last has no `ComplexF64` literal but is a
**mixed-MPS producer**: with a complex `coefficient` only tensor 1 promotes. Allocate all tensors at
one `T`.

This makes the production melt real: `examples/dmt/spin1_melt.jl:291-292` adds a product identity
state to a local-sum state, both from Hermitian operators, and `add(real, real) → Float64`.

Add `operator_real_state(rho::MPS; rtol=sqrt(eps(Float64)))`, throwing unless
`norm(imag) <= rtol*norm` — Kemp's guard, placed at the one boundary that genuinely narrows (an old
checkpoint or a pre-fix complex state the user knows is real).

**1B.4 `tebd_evolve!`** (`tebd.jl:303`) — one uniform promotion rule, covering both directions:

```julia
T = float(promote_type(_mps_eltype(psi), eltype(gate)))
gate_tensor = _dense_local_operator(sites, convert(Matrix{T}, gate))   # cheap: d^(2 span) square
_promote_mps_eltype!(psi, T)                                           # no-op when already T
```

This kills every mixed contraction in the gate stage, including the reverse case (real gate on a
complex state) which pays the promoted-copy penalty today. **`_promote_mps_eltype!` must save and
restore `ITensorMPS.leftlim`/`rightlim`** around its `psi[i] = …` writes, since `setindex!(::MPS,…)`
invalidates the gauge while an eltype conversion does not change it. Those accessors are already
used at `dmt.jl:414-415`.

**1B.5 Leave `operator_expectation_profile` returning `ComplexF64`** — recommended, and it means
**zero tests break** (`test/test_operator_space.jl:816` asserts `isa ComplexF64`). It is not on the
hot path, its `O(N)` environment sweep already runs in the MPS's own eltype, and the trace guard
consumes `abs(...)`/`lognorm(...)` which are `Float64` either way. Instead fix the two wrong
docstrings: `operator_trace` already returns `Float64` for a real MPS while `expectations.jl:14,36`
claim `ComplexF64`.

*Optional, independent:* make `_operator_window_cap` (`expectations.jl:59-82`) return a real cap for
a Hermitian probe while keeping the public return type. `expectations.jl:258` records that this
contraction took a preservation pass from 205 s to 12 s, so 2–4× there is a real test-suite win.
Do it last or never.

---

## Phase 2 — Fuse gate application with the bond factorization

Kemp applies the gate to the two-site product and hands the raw tensor straight to the truncation
(`applyTwoSiteSuperOp` → `truncate_at_site` → `svdBond(AA)`), which performs the **one and only**
factorization. We factorize twice: `product()` inside `tebd_evolve!` splits the gated tensor, then
`_dmt_bond_factorize` re-factorizes `psi[bond]`. Both are at the *inflated* bond (`χ → χd² = 9χ` at
`d=3`), which is why `_dmt_bond_factorize` is 71% of a bond step.

**Design.** The kernel needs a left-orthonormal isometry, a `χ_L × χ_R` bond matrix, and a
right-orthonormal block. From the gated two-site tensor `AA` that is: **QR from the left** (`AA = Q R`,
`Q` left-orthonormal) then **LQ from the right** on the order-3 `R` (`R = L W`, `W` right-orthonormal,
`L` the bond matrix). Two QRs replace one SVD plus one QR, and there is no intermediate
materialization of `psi[bond]`/`psi[bond+1]` at the inflated bond.

- New internal `_dmt_bond_truncate_fused!(psi, bond, AA; kwargs...)`. Keep
  `_dmt_bond_truncate!` untouched for direct callers and `test_dmt_kernel.jl`.
- `dmt_step!` builds `AA = gate_tensor * psi[bond] * psi[bond+1]` itself and takes the fused path
  **only when `gate_maxdim == 0`** (already the default). A finite `gate_maxdim` needs the
  intermediate truncation, so it keeps the current two-step path.
- Falls back to the current path for span-1 gates and the periodic wrap at `n == length(psi)`, which
  `tebd_evolve!` handles specially (`tebd.jl:311-327`).
- `cutoff` semantics are preserved: `product()`'s truncating SVD currently applies it, and
  `_dmt_refactor` applies it at the end either way.

**Gate this on measurement.** Estimated ~1.3× (the gate stage is 39% of a step per
`dmt.jl:386-387`, and fusion removes its factorization but not the gate contraction). If
`dev/bench_dmt.jl` does not show a win after Phase 1, drop it — Phase 1 changes the arithmetic
underneath this estimate.

---

## Phase 3 — Budget off-by-one

Their `subMaxDim = maxDim - (sdimL + sdimR - 1)` is **exactly tight**; ours reserves one more. The
`-1` is real: the connector `C[i,j] = D[i,1] D[1,j] / D[1,1]` satisfies `C[1,j] = D[1,j]` and
`C[i,1] = D[i,1]`, so subtracting it zeroes row 1 and column 1 exactly. The protected part of
`D - C` therefore spans rows `2..sdimL` and columns `2..sdimR` — rank `≤ sdimL + sdimR - 2` — and
adding the rank-one connector back gives `≤ sdimL + sdimR - 1`.

Ours currently lands at `maxdim - 1` total, keeping one fewer complement direction at equal
`maxdim`. Change `_dmt_complement_budget` (`bond.jl:10-12`) to `maxdim - (2 d^(2 radius) - 1)` and
lower the floor in `_validate_dmt_budget` to `2 d^(2 radius)`. Re-run
`test/test_dmt_preservation.jl` and `test/test_dmt_higher_spin.jl` — the at-the-floor preservation
cells are the ones that must still hold at ~1e-13. Update the floors quoted in `PR_BODY.md`
(9/19/33) and `docs/src/manual/dmt.md`.

Worth 0.25% at `χ=400`; do it because the arithmetic should be exactly right, and because the PR
body currently claims our convention matches the reference implementations when it is off by one.

---

## Phase 4 — Selective operator preservation

Kemp's `addPresOperator` / `OnlyPreserveEnergyDensity` (`DMT.h:159-206`, `presAll_ = false`):
protect a *chosen list* of operators instead of the full `d^(2 radius)` block. This is the only route
to `preserve_diameter = 5` at `d = 4`, where our floor is 513 and the blocks are ~0.5 GB. Protecting
~10 operators per side instead of 256 makes that cell reachable.

Upstream it is marked experimental and has `PrintData` debug spam inside the construction loop, so
**design it, don't port it.**

- New keyword `preserve_operators` on `dmt_step!` / `DMTOptions` / `_dmt_bond_truncate!`, defaulting
  to `nothing` = today's full-block behaviour (`presAll_`).
- Accept a vector of `(support, matrix)` pairs in physical space; convert to operator-space
  coefficient vectors via the existing `_operator_coefficients` (`states.jl:77`), then pad with
  identity on the remaining protected sites — the same construction as `DMT.h:168-187`.
- Build the protected block as a QR of `[vec(I) | selected ops]` instead of the full local block.
  `_protected_basis` already does the thin QR; only the input matrix changes. **The identity must
  stay column 1** — `bond.jl:257-259` and the `@assert` at `bond.jl:251-252` both depend on it, and
  the latter must become `size(protected_left, 2) == n_left` rather than `d^(2*left_count)`.
- Budget becomes `maxdim - (n_left + n_right - 1)` — Phase 3's formula with the counts generalized.
- Reject `support > radius` with a clear error, as `DMT.h:219-221` does.

**Verification is the whole point here**, and it is a *weaker* guarantee, so it must be tested as
such: assert that the listed operators are preserved to ~1e-13 **and** that unlisted operators of
the same diameter are *not* (they should break at 1e-2-ish, exactly as
`test_dmt_preservation.jl` already asserts the diameter edge both ways). Then demonstrate the
target cell: `d = 4, preserve_diameter = 5` with energy density + `Sz` preserved, which is
unreachable today.

---

## Phase 5 — Checkpoint / restart

At `χ ≥ 400, d ≥ 4` a single sweep already takes minutes; the melt runs take hours. Kemp checkpoints
`rho` plus the observable time series every N hours and resumes from the same input file, with
`checkpoint_cat.py` concatenating the segments.

- **First verify what already exists**: whether ITensors' HDF5 writer is reachable from this project
  (check `Project.toml` for `HDF5`/`JLD2`, and whether any `src/` code already serializes an MPS).
  Do not add a heavy dependency without checking.
- Save: the MPS, the current time, the sweep index, the accumulated observable series, and the run
  parameters (so a resume can validate it is resuming the same run rather than silently splicing two
  different ones — that validation is the part worth more than the I/O).
- Resume: wall-clock-interval triggered, writing `<name>_t_<time>.{dat,mps}`; on restart, detect the
  latest checkpoint for the given parameters and continue, appending rather than overwriting.
- Scope this to the `examples/dmt/` driver layer plus one small `src/` serialization helper. It is
  infrastructure, not algorithm — it must not touch the kernel.

---

## Verification

**Per phase.** Full suite (`julia --project=. test/runtests.jl`, 79 testsets) green after every
phase. No phase lands on a red suite.

**Phase 1 specifically — four independent proofs:**

1. **Real ≡ complex on a Hermitian operator.** New testset in `test/test_dmt_kernel.jl`, `d in (2,3)`:
   truncate the same base state as `Float64` and as `ComplexF64`, assert relative difference
   `< 1e-12` (measured 1.78e-15 — not bit-identical, since `dgeqrf`/`dgesdd` ≠ `zgeqrf`/`zgesdd`).
   > **Trap:** with `orthogonalize=true` the two runs re-gauge independently and mint different
   > `Index` ids, and the subtraction throws "indices are not permutations of each other". Gauge
   > **once** before copying and pass `orthogonalize=false` — or, better, compare only
   > gauge-invariant quantities: `operator_expectation_profile` over the diameter probes (reuse
   > `preservation_error`, `test/dmt_test_helpers.jl:98`) and `operator_trace`. Ship the
   > gauge-invariant version.
2. **Complex path unchanged.** `test/test_dmt_kernel.jl` must pass with **zero edits** — every
   `lowrank.jl` helper is unit-tested there against complex inputs, and the defaulted-`T` signatures
   exist precisely so those calls are untouched. If a signature change forces an edit there, the
   design is wrong. Then `_OPERATOR_GATE_FORCE_COMPLEX[] = true` must reproduce today's numbers
   exactly across `test_dmt.jl`, `test_dmt_preservation.jl`, `test_dmt_higher_spin.jl`,
   `test_oracles.jl`, `test_transport_reference.jl`.
3. **Preservation guarantee.** Both preservation suites already run the real branch and assert
   `< 1e-11`; measured on the projected real path, 2.53e-15. What is missing is proof the branch runs
   *real* arithmetic — that is 1A.1 plus `_DMT_VERIFY_ELTYPE`.
4. **No mixed MPS anywhere.** After Phase 1, assert in the oracle tests that no code path produces a
   mixed-eltype MPS. Three producers exist today (the kernel itself, `product` with a complex gate,
   `operator_basis_state` with a complex coefficient).

**End to end.** `test/test_oracles.jl` (ED oracle) and `test/test_transport_reference.jl` are the
regression nets. Then run `examples/dmt/spin1_melt.jl`, `examples/dmt/domain_wall_melting.jl` and
`examples/dmt/spin1_semiexact_validation.jl` and confirm the reported `z(t)` ladder and drift figures
are unchanged within their stated floors (1e-13…2e-12).

**Speed.** `dev/bench_dmt.jl` with the eltype axis, **on a quiet box** (`/proc/loadavg` in the
header, BLAS threads pinned). Targets: per-bond kernel ~1.7–1.9×, gate ~2.3–2.8×, full sweep
~1.8–2.1× from Phase 1; a further ~1.3× from Phase 2 if it materializes. Report the `:dense` peak
memory halving — it is what matters for the `d=4, preserve_diameter=5, maxdim=513` ~7 GB case
documented at `dmt.jl:43`.

## Principal risks

| risk | mitigation |
|:--|:--|
| **Partial conversion falls off BLAS, silently.** `hcat` re-promotes at `bond.jl:267-268` with no error, and NDTensors materializes promoted copies on mismatched eltypes. | Single `T` promoted over the whole chain; the `_dmt_bond_solve` barrier; `_DMT_VERIFY_ELTYPE`; the 1A.1 eltype assertions. |
| Silent loss of a genuine imaginary component | `T` is a widening, so impossible inside the kernel. Three narrowing boundaries only, each guarded by a relative tolerance and each backed by the Hermitian-basis theorem. `convert(Float64, ::Complex)` throws `InexactError` on a nonzero imaginary part. |
| Mixed-eltype MPS from an existing producer | `_mps_eltype` is correct on mixed input; the kernel writes both mutated tensors at `T`; `tebd_evolve!` promotes both ways; state builders allocate uniformly. |
| ITensors promoting behind our back in a future version | Verified it does not today for `qr`, `svd`, `product`, `orthogonalize!`, `add`, `random_mps`, `ITensor(data, inds...)` on ITensors 0.9.30 / ITensorMPS 0.4.1. `Project.toml` floats `ITensorMPS = "0.3.44, 0.4"`, so the 1A.1 assertions are the tripwire. |
| `:random` truncation reproducibility | `randn(T, …)` consumes half as many RNG draws, so a seeded run differs between paths. Not a regression — `:random` is already documented non-deterministic (`dmt.jl:31-38`) — but say so in `lowrank.jl:101-111`. |
| Phase 2 regresses accuracy by moving where `cutoff` applies | Compare against the Phase 1 result on the same state; `_dmt_refactor` must remain the single place `cutoff` binds. |
| Phase 4 weakens the guarantee by design | Test the negative direction explicitly: unlisted same-diameter operators must break. A test that only checks the listed ones would pass on a no-op implementation. |

## Critical files

Kernel: `src/operator_space/dmt/bond.jl`, `src/operator_space/dmt/lowrank.jl`,
`src/operator_space/dmt.jl`
Boundaries: `src/operator_space/gates.jl`, `src/operator_space/states.jl`,
`src/operator_space/thermal.jl`, `src/evolution/tebd.jl`
Measurement: `src/operator_space/expectations.jl` (docstrings only)
Tests: `test/test_dmt_kernel.jl` (**must pass unedited**), `test/dmt_test_helpers.jl`,
`test/test_dmt_preservation.jl`, `test/test_dmt_higher_spin.jl`, `test/test_oracles.jl`
Benchmarks: `dev/bench_dmt.jl`, `examples/dmt/spin1_semiexact_validation.jl`
Docs: `docs/src/manual/dmt.md`, `PR_BODY.md`

Reference source, cloned for consultation:
`/tmp/claude-639688/-data-djxg096/045483b4-bab1-4ed3-9ecd-2c1d21abc057/scratchpad/kemp-dmt/`
(`DMT.h:719-720` budget, `:702-703,758` connector, `:154,201,693-694` protected QR, `:538-544`
real cast, `:159-206` selective preservation), with line-cited notes in `kemp-siteset-findings.md`
and `kemp-gates-findings.md` alongside it.
