# PXP Energy-Transport Infrastructure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** PXP model helpers, constraint MPOs in state and operator space, imaginary-time domain-wall preparation, expectation sweeps, and a checkpointed (constraint-projected) DMT evolution driver, with an end-to-end energy-transport example.

**Architecture:** Generic Pauli-vectorization converters (`pauli_state_from_mpo`, `pauli_superoperator_mpo`) lift the bond-dimension-2 physical constraint MPO into operator space; thermal preparation and measurement are plain operator-space TEBD plus identity-environment sweeps; the constrained driver chunks `dmt_evolve!` between projector applications. Spec: `docs/superpowers/specs/2026-06-12-pxp-energy-transport-design.md`.

**Tech Stack:** Julia, ITensors/ITensorMPS, existing MPSToolkit TEBD/DMT machinery.

**Conventions (used throughout):** local basis index 1 = `|0⟩` (ground), 2 = `|1⟩` (excited); `P = (I+Z)/2`, `n = (I−Z)/2`; kron factor order = chain site order; normalized Pauli strings `P_α = σ_α/√2` per site; vectorized coefficients `c_α = tr(P_α† O)`; evolution `ρ(t) = e^{−iHt} ρ e^{+iHt}`.

Run tests per-file during development with:
`julia --project=. test/test_pxp.jl` (file carries its own `using` block, like `test_dmt.jl`).

---

### Task 1: PXP model helpers (`pxp_term_support`, `pxp_term_hamiltonian`, `pxp_constraint_mpo`)

**Files:**
- Create: `src/models/pxp.jl`
- Modify: `src/MPSToolkit.jl` (include after `models/spinhalf.jl`; extend `Models` namespace + top-level exports)
- Create: `test/test_pxp.jl` (add to `test/runtests.jl` after `test_dmt.jl`)

- [ ] **Step 1: Failing tests** — dense oracles in `test/test_pxp.jl`:

```julia
using ITensors, ITensorMPS, LinearAlgebra, MPSToolkit, Test

const _ID2 = Matrix{ComplexF64}(I, 2, 2)

function _embed_term(op, start, nsites)
  span = round(Int, log2(size(op, 1)))
  left = start == 1 ? Matrix{ComplexF64}(I, 1, 1) : foldl(kron, fill(_ID2, start - 1))
  right_count = nsites - (start + span - 1)
  right = right_count == 0 ? Matrix{ComplexF64}(I, 1, 1) : foldl(kron, fill(_ID2, right_count))
  return kron(left, ComplexF64.(op), right)
end

_pxp_dense(nsites; omega=1.0) = sum(
  _embed_term(pxp_term_hamiltonian(nsites, j; omega=omega), first(pxp_term_support(nsites, j)), nsites)
  for j in 1:nsites)

# diagonal dense P_G; site 1 = most significant bit
function _pxp_projector_dense(nsites)
  d = 2^nsites
  diagvals = ones(Float64, d)
  for k in 0:(d - 1)
    bits = digits(k; base=2, pad=nsites) |> reverse   # bits[1] = site 1
    any(bits[i] == 1 && bits[i + 1] == 1 for i in 1:(nsites - 1)) && (diagvals[k + 1] = 0.0)
  end
  return Matrix{ComplexF64}(Diagonal(diagvals))
end

function _mpo_dense(op::MPO, sites)
  contracted = op[1]
  for n in 2:length(op); contracted *= op[n]; end
  arr = Array(contracted, prime.(sites)..., sites...)
  n = length(sites)
  perm = vcat(reverse(1:n), reverse((n + 1):(2n)))
  return reshape(permutedims(arr, perm), 2^n, 2^n)
end
```

Tests: `pxp_term_hamiltonian` equals hand-built `kron(X,P)` / `kron(P,X,P)` / `kron(P,X)`;
`omega` scales; `pxp_term_support` ranges; ArgumentErrors for `nsites<2`, `site` out of range;
`_mpo_dense(pxp_constraint_mpo(sites), sites) ≈ _pxp_projector_dense(N)` for N=2..6;
projector idempotent; `norm(H*P − P*H) ≈ 0` and per-term `[h_j, P] = 0`.

- [ ] **Step 2: Run, verify failure** (`UndefVarError: pxp_term_hamiltonian`).

- [ ] **Step 3: Implement `src/models/pxp.jl`:**

```julia
function pxp_term_support(nsites::Integer, site::Integer)
  nsites >= 2 || throw(ArgumentError("PXP requires at least two sites"))
  1 <= site <= nsites || throw(ArgumentError("PXP term site must lie in 1:nsites"))
  return max(Int(site) - 1, 1):min(Int(site) + 1, Int(nsites))
end

function pxp_term_hamiltonian(nsites::Integer, site::Integer; omega::Real=1.0)
  support = pxp_term_support(nsites, site)
  paulis = pauli_matrices()
  ground = (paulis.I + paulis.Z) / 2
  factors = [s == site ? paulis.X : ground for s in support]
  return omega * foldl(kron, factors)
end

_pxp_constraint_transition(state::Int, local_state::Int) =
  local_state == 1 ? (1, 1.0) : (state == 1 ? (2, 1.0) : (2, 0.0))

function pxp_constraint_mpo(sites)
  all(s -> dim(s) == 2, sites) || throw(ArgumentError("PXP constraint MPO requires dimension-2 sites"))
  return _diagonal_projector_mpo(sites, 2, _pxp_constraint_transition)
end
```

(`_diagonal_projector_mpo` is the existing automaton-MPO builder in `src/operator_space/daoe.jl`;
it is site-dimension generic. Docstrings on all three functions, repo style.)

- [ ] **Step 4: Wire include + exports, run `julia --project=. test/test_pxp.jl` → PASS.**
- [ ] **Step 5: Commit** `feat(models): add PXP term and constraint-MPO helpers`.

---

### Task 2: Pauli vectorization converters + PXP operator-space wrappers

**Files:**
- Create: `src/operator_space/vectorize.jl`, `src/operator_space/pxp.jl`
- Modify: `src/MPSToolkit.jl`
- Test: `test/test_pxp.jl`

- [ ] **Step 1: Failing tests:**

```julia
const _PAULI_NORM = [m / sqrt(2) for m in values(pauli_matrices())]
_pauli_string(labels) = foldl(kron, (_PAULI_NORM[l] for l in labels))

# pauli_state_from_mpo: amplitudes are tr(P_α† O)
# (a) PXP constraint MPO, N = 4; (b) non-Hermitian OpSum MPO ("Sz" 1; 0.7 "S+" 2), N = 2.
for labels in Iterators.product(ntuple(_ -> 1:4, N)...)
  @test inner(pauli_basis_state(psites, collect(labels)), vecO) ≈ tr(_pauli_string(labels)' * Odense) atol=1e-10
end

# pauli_superoperator_mpo: column-by-column dense superoperator oracle, N = 2 and 3:
# apply(S, |α⟩) amplitudes β must equal tr(P_β† M P_α M†).
# pauli_pxp_constraint_state: bond dim 2; pauli_trace equals count of allowed configs (Fibonacci F(N+2)).
# pauli_pxp_constraint_projector: idempotent on a random vectorized operator; fixes the constraint state;
# annihilates pauli_state_from_mpo of |11⟩⟨11|-supported operator.
```

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement.** `vectorize.jl`:

```julia
function _mpo_unprimed_site_index(op::MPO, n::Int)
  links = ITensors.Index[]
  n > 1 && push!(links, commonind(op[n], op[n - 1]))
  n < length(op) && push!(links, commonind(op[n], op[n + 1]))
  remaining = filter(i -> plev(i) == 0 && all(l -> i != l, links), collect(inds(op[n])))
  length(remaining) == 1 || throw(ArgumentError("could not infer the unprimed MPO site index at site $(n)"))
  return only(remaining)
end

function pauli_state_from_mpo(op::MPO, sites)
  length(op) == length(sites) || throw(ArgumentError(...))
  tensors = ITensor[]
  for n in 1:length(op)
    phys = _mpo_unprimed_site_index(op, n)
    dim(phys) == 2 || throw(ArgumentError(...)); dim(sites[n]) == 4 || throw(ArgumentError(...))
    conv = ITensor(ComplexF64, sites[n], prime(phys), phys)   # conv[α,s',s] = conj(P_α[s',s])
    ... fill from _pauli_basis_operators(1) ...
    push!(tensors, op[n] * conv)
  end
  return MPS(tensors)
end

function pauli_superoperator_mpo(op::MPO, sites)
  # per bond: bra link = sim(ket link); per site:
  #   ket = op[n]  (indices l, r, s', s)
  #   bra = replaceinds(dag(op[n]), (s', s, l, r) -> (t', t, lb, rb))
  #   Cin[α, s, t]  = P_α[s, t]          (input leg = sites[n])
  #   Cout[α', s', t'] = conj(P_α'[s', t'])  (output leg = prime(sites[n]))
  #   tensor = ket * bra * Cin * Cout
  # then fuse (r, rb) on every bond with one shared combiner applied to both neighbors.
  ...
  return MPO(tensors)
end
```

`pxp.jl` (operator space): wrappers creating internal `Index(2, "PXPSite,n=…")` physical sites:

```julia
pauli_pxp_constraint_state(sites)     = pauli_state_from_mpo(pxp_constraint_mpo(phys), sites)
pauli_pxp_constraint_projector(sites) = pauli_superoperator_mpo(pxp_constraint_mpo(phys), sites)
```

- [ ] **Step 4: Run tests → PASS.**
- [ ] **Step 5: Commit** `feat(operator-space): add Pauli vectorizers and PXP constraint objects`.

---

### Task 3: Expectation helpers (`pauli_trace`, `pauli_expectation`, `pauli_expectation_profile`)

**Files:** Create `src/operator_space/expectations.jl`; modify `src/MPSToolkit.jl`; test in `test/test_pxp.jl`.

- [ ] **Step 1: Failing tests.** Build a random dense Hermitian `ρ_d` on N=3 (`R + R†` over the
normalized Pauli strings), assemble the matching MPS as `foldl(+, c_α · pauli_basis_state(α))`,
then for windows (1-site at 3, 2-site at 1 and 2, 3-site at 1) with random Hermitian `op`:
`pauli_expectation(rho, op, start)` ≈ `tr(ρ_d · O_embed)/tr(ρ_d)`; `normalize=false` variant ≈
`tr(ρ_d · O_embed)`; `pauli_trace(rho)` ≈ `tr(ρ_d)`; `pauli_expectation_profile` over the PXP
term list equals per-term `pauli_expectation` values and accepts unsorted windows.

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement.** Cap tensor `cap[α…] = tr(P_α O)` over the window (loop
`Iterators.product`, `foldl(kron, …)` per label, skip zeros); identity caps elsewhere via the
existing `_pauli_identity_env`; one O(N) pass with cumulative left environment and a
precomputed right environment array (same pattern as the XXZ profile sweep);

```
tr(ρ O)   = 2^((N − span)/2) · ⟨caps⟩,   tr(ρ) = 2^(N/2) · ⟨identity caps⟩,
normalized value = ⟨cap sweep⟩ / (2^(span/2) · ⟨identity sweep⟩)
```

Span inference from the dense op via existing `_spinhalf_span` (physical 2^span, not 4^span).
Values returned as `ComplexF64`; Hermitian inputs give real values (documented; tests assert
vanishing imaginary part).

- [ ] **Step 4: Run tests → PASS.**
- [ ] **Step 5: Commit** `feat(operator-space): trace and local-expectation sweeps for Pauli MPS`.

---

### Task 4: Imaginary-time thermal/domain-wall preparation

**Files:** Create `src/operator_space/thermal.jl`; modify `src/MPSToolkit.jl`; test in `test/test_pxp.jl`.

- [ ] **Step 1: Failing tests.**
(a) `pauli_gate_from_imaginary_time(h, db)` ≈ `pauli_gate(exp(-(db/2)h))` and equals the
identity matrix at `db = 0`.
(b) N=4 PXP, uniform `weights = fill(β, 4)`, `β = 0.6`, `initial_state = pauli_pxp_constraint_state`:
profile from `pauli_gibbs_state(...; nsteps=20, maxdim=256, cutoff=0.0)` matches dense
`ρ = e^{−βH/2} P_G e^{−βH/2}` expectations `tr(ρ h_j)/tr(ρ)` to `atol=1e-3`, and `nsteps=20`
beats `nsteps=3` (Trotter convergence).
(c) Domain wall N=4, wall 2|3, `weights = [β, 0, 0, −β]`: matches dense
`e^{−K/2} P_G e^{−K/2}`, `K = β h_1 − β h_4`.
(d) Identity default: `initial_state=nothing` reproduces dense `e^{−K/2} 𝟙 e^{−K/2}`.

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement:**

```julia
pauli_gate_from_imaginary_time(h, dbeta) = pauli_gate(exp(-(dbeta / 2) * Matrix{ComplexF64}(h)))

function pauli_gibbs_state(sites, terms, weights; nsteps=4, maxdim=64, cutoff=1e-12, initial_state=nothing)
  # gates sandwich e^{-(w_j/(4 nsteps)) h_j} per application; forward+reverse sweep per step
  # → 2·nsteps applications → e^{-(w_j/2) h_j} per side in total.
  gates = [pauli_gate_from_imaginary_time(h, w / (2 * nsteps)) for ((_, h), w) in zip(terms, weights)]
  rho = isnothing(initial_state) ? pauli_basis_state(sites, fill(1, length(sites))) : copy(initial_state)
  for _ in 1:nsteps
    for i in eachindex(gates);          tebd_evolve!(rho, gates[i], first_site(terms[i]); maxdim, cutoff); end
    for i in reverse(eachindex(gates)); tebd_evolve!(rho, gates[i], first_site(terms[i]); maxdim, cutoff); end
  end
  normalize!(rho)
  return rho
end
```

with validation (`length(terms) == length(weights)`, `nsteps ≥ 1`, site-index agreement when
`initial_state` is passed). `terms` entries are `(start, h)` pairs.

- [ ] **Step 4: Run tests → PASS.**
- [ ] **Step 5: Commit** `feat(operator-space): imaginary-time Gibbs/domain-wall preparation`.

---

### Task 5: Constrained (checkpointed) DMT driver

**Files:** Create `src/operator_space/constrained.jl`; modify `src/MPSToolkit.jl`; test in `test/test_pxp.jl`.

- [ ] **Step 1: Failing tests.** N=6 PXP ED oracle:
domain wall `β=0.4`, wall 3|4, weights `[β, β, 0, 0, −β, −β]` (terms fully inside a half keep
±β; straddlers get 0); prepare with `pauli_gibbs_state` (`nsteps=12`, `maxdim=256`, `cutoff=0`);
dense reference `ρ(t) = e^{−iHt} ρ_0 e^{+iHt}` with `ρ_0 = e^{−K/2} P_G e^{−K/2}`.
Evolve MPS with per-term gates `pauli_gate_from_hamiltonian(h_j, dt)`, `dt=0.05`, schedule =
term starts `[1,1,2,3,4,5]`, `nstep=5` (one forward+reverse sweep advances `2·dt`, so `t=0.5`),
`maxdim=64` (exact at N=6), via
`constrained_dmt_evolve!(rho, evo, projector; project_every=2)`.
Assert: energy profiles match dense to `atol=5e-3`; leakage
`1 − ⟨ρ, 𝒫ρ⟩/⟨ρ,ρ⟩` ≈ 0 after the run; `project_every` validation throws on 0.

- [ ] **Step 2: Run, verify failure.**
- [ ] **Step 3: Implement:**

```julia
function constrained_dmt_evolve!(rho::MPS, evo::DMTGateEvolution, projector::MPO;
    project_every::Integer=1, projector_maxdim::Integer=evo.gate_maxdim, projector_cutoff::Real=evo.cutoff)
  project_every >= 1 || throw(ArgumentError(...))
  length(projector) == length(rho) || throw(ArgumentError(...))
  remaining = evo.nstep
  while remaining > 0
    chunk = min(Int(project_every), remaining)
    chunk_evo = DMTGateEvolution(evo.gate, evo.dt; schedule=evo.schedule,
      reverse_schedule=evo.reverse_schedule, nstep=chunk, maxdim=evo.maxdim, cutoff=evo.cutoff,
      gate_maxdim=evo.gate_maxdim, connector_buffer=evo.connector_buffer)
    dmt_evolve!(rho, chunk_evo)
    rho[:] = apply(projector, rho; maxdim=Int(projector_maxdim), cutoff=projector_cutoff)
    normalize!(rho)
    remaining -= chunk
  end
  return rho
end
```

Default `projector_maxdim = evo.gate_maxdim` is intentional: the projector may inflate the
bond like a raw gate; the next DMT sweep performs the transport-aware compression.

- [ ] **Step 4: Run tests → PASS (tune atol to observed Trotter error with margin).**
- [ ] **Step 5: Commit** `feat(operator-space): checkpointed constraint-projected DMT driver`.

---

### Task 6: Example script, docs, full-suite verification

**Files:**
- Create: `examples/operator_space/pxp_energy_transport.jl`
- Modify: `docs/src/manual/dmt.md` (new worked-example section), `docs/src/api.md`,
  `docs/src/manual/operator-space.md` (links/API lists), `test/runtests.jl` (already done in Task 1),
  README example/feature lists if they enumerate scripts.

- [ ] **Step 1: Example script.** N=40, wall 20|21, `β=0.25`, `dt=0.1`, `maxdim=32`,
`gate_maxdim=128`, `connector_buffer=8`, `project_every=1`, ~50 measurement sweeps
(`t_max ≈ 10`): build terms via `pxp_term_hamiltonian`/`pxp_term_support`; weights `±β` for
terms fully inside a half, `0` for the two straddlers (state factorizes exactly as
`e^{−βH_L} ⊗ e^{+βH_R}`); prepare on top of `pauli_pxp_constraint_state`; evolve with
`constrained_dmt_evolve!`; measure `real.(pauli_expectation_profile(...))`, transferred energy
`ΔE(t) = Σ_{j>wall} [e_j(t) − e_j(0)] ∼ t^{1/z}`, residual leakage, and a one-off
"unprojected-sweep leakage" diagnostic on a copied state; log-log fit of `|ΔE|` on an
intermediate window with an edge-contamination guard; print measured `z` against the
superdiffusive reference `z ≈ 3/2` (PRX 13, 011033 (2023)) with finite-size caveats.
- [ ] **Step 2: Run the example with quick parameters** (small N/maxdim override block at top,
e.g. N=16, 10 sweeps) to confirm it executes end to end.
- [ ] **Step 3: Docs.** `dmt.md`: "Worked example: constrained energy transport in the PXP
chain" — constraint physics, why truncation leaks out of the sector even though `[H, P_G]=0`,
the checkpoint pattern, the domain-wall boundary choice, condensed code; API entries for all
new exports in `api.md` symbol lists + `@docs` blocks (`pxp-*` under Models; `pauli_*`,
`constrained_dmt_evolve!` under OperatorSpace); cross-link from `operator-space.md`.
- [ ] **Step 4: Full verification.** `julia --project=. -e 'using Pkg; Pkg.test()'` → all green
(docstring coverage test included). Re-run `julia --project=. test/test_pxp.jl` standalone.
- [ ] **Step 5: Commit** `feat(examples,docs): PXP energy-transport example and manual section`.

---

## Self-review

- Spec coverage: model helpers (Task 1), constraint MPO state+operator space (Tasks 1–2),
  domain wall via imaginary-time TEBD with boundary treatment (Task 4), checkpointed
  Heisenberg evolution (Task 5), measurement (Task 3), example/docs (Task 6). ✓
- All steps carry concrete code or exact oracles; remaining `...` markers are docstring/
  validation boilerplate spelled out in the spec. ✓
- Naming consistent across tasks (`pxp_term_hamiltonian`, `pauli_gibbs_state`,
  `constrained_dmt_evolve!`, `(start, h)` term pairs). ✓
