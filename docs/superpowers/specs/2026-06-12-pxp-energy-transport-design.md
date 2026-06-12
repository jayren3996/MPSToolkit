# PXP energy-transport infrastructure (DMT) — design

Date: 2026-06-12
Status: approved for implementation

## Goal

Add the infrastructure needed to compute energy transport in the kinetically constrained
PXP chain with the existing operator-space DMT backend:

1. PXP model helpers in `src/models/`.
2. The PXP constraint projector as a physical-space MPO, and the same constraint as an
   **operator-space** superoperator MPO (acting on vectorized density matrices).
3. An energy domain-wall density matrix `ρ ∝ exp(-β H_L) ⊗ exp(+β H_R)` prepared by
   imaginary-time TEBD directly in operator space, with an explicit, controllable treatment
   of the domain boundary.
4. A DMT evolution driver that Heisenberg-evolves the density matrix and periodically
   re-applies the operator-space constraint projector ("checkpoints") to remove
   truncation-induced leakage out of the constrained sector.
5. Operator-space expectation helpers (trace, local expectation, O(N) profile sweep) so the
   energy-density profile can be measured efficiently, plus an end-to-end example script.

## Physics conventions

- Local basis: index 1 = `|0⟩` (ground), index 2 = `|1⟩` (Rydberg excitation).
  `P = |0⟩⟨0| = (I + Z)/2`, `n = |1⟩⟨1| = (I − Z)/2` with `Z = diag(1, −1)`.
- Open-chain PXP Hamiltonian with the standard boundary terms,

  ```
  H = Ω [ X_1 P_2 + Σ_{j=2}^{N−1} P_{j−1} X_j P_{j+1} + P_{N−1} X_N ] ,
  ```

  one term per site `j`; bulk terms are 3-site, the two edge terms are 2-site.
- Constrained subspace: `P_G = Π_{j=1}^{N−1} (1 − n_j n_{j+1})` (no adjacent excitations).
  Every PXP term commutes with `P_G` individually (flipping site `j` requires both
  neighbors in `|0⟩`, which can neither create nor destroy a violation), so exact evolution
  preserves the sector. Truncation does not — hence the checkpoints.
- Density-matrix evolution convention: `ρ(t) = e^{−iHt} ρ e^{+iHt}`, which is exactly the
  superoperator produced by the existing `pauli_gate_from_hamiltonian(h, dt)`.
- The infinite-temperature state of the **constrained** sector is `ρ_∞ ∝ P_G` (not the
  identity), so PXP transport states must be built on top of the vectorized `P_G`.

## Components

### 1. `src/models/pxp.jl` — model helpers

- `pxp_term_hamiltonian(nsites, site; omega=1.0)` → dense matrix of the PXP term assigned
  to `site` (4×4 at the edges, 8×8 in the bulk). Requires `nsites ≥ 2`,
  `1 ≤ site ≤ nsites`.
- `pxp_term_support(nsites, site)` → `UnitRange` of physical sites the term acts on
  (`max(site−1,1):min(site+1,nsites)`); `first(support)` is the TEBD/DMT gate start.
- `pxp_constraint_mpo(sites)` → bond-dimension-2 diagonal MPO of `P_G` on physical
  (dimension-2) site indices, built from the two-state automaton
  ("was the previous site excited"): blocks `A_11 = |0⟩⟨0|`, `A_12 = |1⟩⟨1|`,
  `A_21 = |0⟩⟨0|`, `A_22 = 0`.

The sum convention matches `spinhalf_tfim_bond_hamiltonian`: summing
`pxp_term_hamiltonian` over all `site` values embedded at their supports reproduces the
open-chain `H` exactly. Dense full-chain builders stay in the tests (repo convention).

### 2. `src/operator_space/vectorize.jl` — generic vectorization

Two generic converters (reusable beyond PXP) keep the physical MPO and its operator-space
forms provably consistent:

- `pauli_state_from_mpo(op::MPO, sites)` → Pauli-basis MPS with coefficients
  `c_α = tr(P_α† O)`, where `P_α` are the normalized Pauli strings (`σ/√2` per site).
  Reuses the MPO link structure, so a bond-dimension-k MPO vectorizes at bond dimension k.
- `pauli_superoperator_mpo(op::MPO, sites)` → operator-space MPO implementing
  `ρ ↦ M ρ M†`. Local tensor: `S[α', α] = tr(P_α'† M_loc P_α M_loc'†)` contracted over the
  ket copy and the conjugated bra copy of each MPO tensor; the (ket, bra) link pairs are
  fused with combiners, so bond dimension k lifts to k². Output prime convention matches
  `_diagonal_projector_mpo` (`prime(site), site`), so `ITensorMPS.apply` works unchanged.

Alternatives considered: (a) hand-built automaton for each operator-space object — less
code per object but duplicates the constraint logic in two places with no cross-check;
(b) dense superoperator construction — exponential in N. The generic converters are ~80
lines total, testable against dense oracles at small N, and make the PXP wrappers
one-liners. Chosen: (generic converters + thin wrappers).

### 3. `src/operator_space/pxp.jl` — PXP constraint in operator space

- `pauli_pxp_constraint_state(sites)` → vectorized `P_G` as a Pauli MPS (bond dimension 2):
  the constrained infinite-temperature initial state, and the natural `initial_state` for
  thermal preparation. Implemented as `pauli_state_from_mpo(pxp_constraint_mpo(...))` over
  internally created physical sites.
- `pauli_pxp_constraint_projector(sites)` → operator-space MPO for `ρ ↦ P_G ρ P_G`
  (bond dimension 4), via `pauli_superoperator_mpo`. Idempotent up to truncation; leaves
  any operator supported in the constrained sector exactly invariant.

### 4. `src/operator_space/thermal.jl` — imaginary-time preparation

- `pauli_gate_from_imaginary_time(h, dbeta)` → dense Pauli superoperator for
  `ρ ↦ e^{−(dbeta/2) h} ρ e^{−(dbeta/2) h}` (i.e. `pauli_gate(exp(−(dbeta/2) h))`; works
  because `pauli_gate(A)` implements conjugation `A · A†` for arbitrary `A`, and
  `e^{−τh}` is Hermitian).
- `pauli_gibbs_state(sites, terms, weights; nsteps=4, maxdim=64, cutoff=1e-12, initial_state=nothing)`
  → Pauli MPS for `ρ ∝ e^{−K/2} ρ_0 e^{−K/2}` with `K = Σ_j w_j h_j`, built by Trotterized
  imaginary-time TEBD in operator space:
  - `terms`: vector of `(start, h)` pairs (dense local Hermitian blocks),
  - `weights`: one real weight per term (`+β`, `−β`, `0`, or any profile),
  - per Trotter step the per-term gates `pauli_gate_from_imaginary_time(h_j, w_j/nsteps)`
    are applied in a forward-then-reverse sweep (symmetrized splitting), truncated by
    ordinary TEBD SVD (`maxdim`, `cutoff`) — DMT bias is unnecessary here because the
    state stays low-entanglement at small β,
  - `initial_state=nothing` means the identity product state (unconstrained infinite-T);
    PXP callers pass `pauli_pxp_constraint_state(sites)`,
  - the result is `normalize!`d (2-norm); physical observables normalize by trace anyway.

**Domain boundary.** The wall is encoded entirely in `weights`, so the caller controls the
boundary treatment. The example uses, for a wall between sites `N/2` and `N/2+1`:
`w_j = +β` for terms supported entirely left of the wall, `−β` entirely right, and `0` for
the (two) wall-straddling terms. With that choice `ρ(0) = e^{−βH_L} ⊗ e^{+βH_R}` exactly
(no term couples the halves, so the two exponentials factorize), which is the literal
domain-wall product state; the alternative (assign straddling terms by their center site)
differs only by an `O(β)` local distortion at the wall and gives the same hydrodynamics.
Because each `h_j` commutes with `P_G`, preparation starting from the vectorized `P_G`
stays in the constrained sector up to truncation error.

### 5. `src/operator_space/expectations.jl` — measurement

- `pauli_trace(rho)` → `tr(ρ) = (√2)^N · c_{I…I}` via one identity-cap contraction sweep.
- `pauli_expectation(rho, op, start; normalize=true)` → `tr(ρ O)/tr(ρ)` (or unnormalized
  `tr(ρ O)`) for a dense local `op` whose span is inferred from its size. The window cap
  carries `tr(P_α O)` per multi-site Pauli label; contraction uses explicit identity-cap
  environments (same pattern as the XXZ example's profile sweep), not `inner`, to keep
  conjugation conventions out of the picture.
- `pauli_expectation_profile(rho, terms; normalize=true)` → vector of `tr(ρ h_j)/tr(ρ)`
  for a list of `(start, h)` pairs, computed in one O(N) pass with cumulative left/right
  identity environments. This is the energy-density profile measurement.

### 6. `src/operator_space/constrained.jl` — checkpointed DMT driver

- `constrained_dmt_evolve!(rho, evo::DMTGateEvolution, projector::MPO; project_every=1, projector_maxdim=evo.gate_maxdim, projector_cutoff=evo.cutoff)`
  Runs `evo.nstep` forward-plus-reverse DMT sweeps in chunks of `project_every`; after each
  chunk (and at the end) applies `projector` with `ITensorMPS.apply` and `normalize!`s.
  Default `projector_maxdim = evo.gate_maxdim` is deliberate: the projector application is
  allowed to inflate the bond (like a raw gate application), and the **next DMT sweep**
  performs the transport-aware compression, instead of an ordinary SVD making the
  truncation decision at the checkpoint. Implementation rebuilds a `DMTGateEvolution` with
  `nstep=chunk` from `evo`'s fields per chunk.

Alternatives considered: example-script-only loop (not reusable/testable as a unit — the
user asked for this as functionality) and a new dispatchable evolution type wrapping
`DMTGateEvolution` (more ceremony than a function warrants today). Chosen: driver function.

## Example: `examples/operator_space/pxp_energy_transport.jl`

End-to-end script mirroring `xxz_transport_regimes.jl`:

1. Build PXP terms via `pxp_term_hamiltonian` / `pxp_term_support`.
2. Prepare the β domain wall on top of the vectorized `P_G` with `pauli_gibbs_state`
   (small β, linear-response regime; wall-straddling weights zero).
3. Evolve with `constrained_dmt_evolve!` (3-site bulk gates, 2-site edge gates, sequential
   forward/reverse schedule), checkpointing the constraint projector.
4. Measure the energy profile with `pauli_expectation_profile`; track the transferred
   energy `ΔE(t) = Σ_{j > wall} [⟨h_j⟩(t) − ⟨h_j⟩(0)] ∼ t^{1/z}` and the constraint
   leakage `1 − ⟨ρ, P_G ρ P_G⟩/⟨ρ, ρ⟩`.
5. Fit `1/z` on an intermediate window (same plateau logic as the XXZ example) and print it
   against the literature value `z ≈ 3/2` (superdiffusive energy transport in PXP;
   Ljubotina, Desaules, Serbyn, et al., PRX 13, 011033 (2023)), with honest commentary
   that small `N`/`maxdim` give a rough exponent.

## Testing (`test/test_pxp.jl`, wired into `runtests.jl`)

Dense exact-diagonalization oracles at small N throughout, following `test_oracles.jl`:

1. `pxp_term_hamiltonian`/`pxp_term_support` against hand-built `kron` references; sum of
   embedded terms equals the dense open-chain `H`; argument validation.
2. `pxp_constraint_mpo` contracts to the dense `Π(1 − n n)` projector for N = 2…6;
   idempotent; `[H_dense, P_dense] = 0`.
3. `pauli_state_from_mpo`: amplitudes equal `tr(P_α† O)` for the constraint MPO and for a
   non-Hermitian random-ish MPO; `pauli_trace` matches the dense trace.
4. `pauli_superoperator_mpo`: applying to every Pauli basis state matches the dense
   superoperator matrix `S[β, α] = tr(P_β† M P_α M†)` (N = 2, 3); covers complex blocks.
5. `pauli_pxp_constraint_projector`: idempotent; annihilates vectorized operators with
   adjacent-excitation support; fixes `pauli_pxp_constraint_state` exactly.
6. `pauli_gibbs_state`: uniform-β PXP Gibbs state at N = 4 matches dense
   `e^{−βH/2} ρ_0 e^{−βH/2}` expectations within Trotter error and converges as `nsteps`
   grows; domain-wall weights reproduce the dense factorized state.
7. `pauli_expectation`/`pauli_expectation_profile` against dense traces on random Pauli
   MPS, including edge windows and unnormalized mode.
8. `constrained_dmt_evolve!` on an N = 6 PXP domain wall: energy profile after short
   evolution matches dense `ρ(t) = U ρ U†` within Trotter+truncation tolerance; leakage
   stays at numerical zero; `project_every` chunking equals an unchunked run when the
   projector is a no-op on constrained states.

## Documentation

- New worked-example section in `docs/src/manual/dmt.md` (constrained energy transport in
  PXP) describing the constraint-leakage problem and the checkpoint pattern.
- API entries for all new exported symbols in `docs/src/api.md` and the operator-space
  manual page; example listed alongside the existing operator-space scripts.
- New exports surfaced through the `Models` and `OperatorSpace` namespaces in
  `src/MPSToolkit.jl`, mirroring the existing export style.

## Out of scope

- Periodic boundary conditions (DMT explicitly rejects PBC local gates today).
- Other constrained models (East, Fredkin, higher-spin PXP generalizations) — the generic
  vectorizers make these easy follow-ups but none is built here.
- Current-current (Kubo) formulations of transport; the domain-wall protocol is what the
  user asked for.
