# Operator Space

*Vectorized operator dynamics over a generalized Gell-Mann basis of any local dimension `d` (Pauli at `d = 2`): the Heisenberg-picture and open-system foundation that DMT, DAOE, and Lindblad TEBD are built on.*

## Background

Most of `MPSToolkit.jl` evolves a *state* ``|\psi\rangle`` forward in time. The operator-space layer instead evolves an *operator* — a Heisenberg observable, a density matrix, or a Lindbladian dissipator — by encoding it as an MPS over a local operator basis of dimension ``d^2`` for a physical local Hilbert space of dimension ``d`` (``d = 2S+1`` for spin ``S``). This makes autocorrelation functions, transport coefficients, and open-system steady states accessible with the same finite-MPS machinery used for pure states.

The API comes in two layers: `operator_*` is generic in `d` (`operator_siteinds`, `operator_basis_state`, `operator_gate`, ...), and every `pauli_*` name (`pauli_siteinds`, `pauli_basis_state`, `pauli_gate`, ...) is a thin `d = 2` wrapper around the matching `operator_*` function, kept because it is the common spin-1/2 case and because it predates the generic layer. The two are numerically interchangeable at `d = 2` — same basis, same conventions — so this page documents `operator_*` first and calls out each `d = 2` specialization alongside it.

### Heisenberg vs Schrödinger picture

In the Schrödinger picture the state carries the time dependence, ``\rho(t) = e^{-iHt}\,\rho(0)\,e^{iHt}``, and observables are fixed. In the Heisenberg picture the roles swap: the state is fixed and an observable ``O`` evolves as

```math
O(t) = e^{iHt}\,O\,e^{-iHt}, \qquad \dot O = i[H, O].
```

For computing a single correlator ``\langle A(t)\,B\rangle`` the Heisenberg picture is often the more economical choice — you evolve one operator rather than re-evolving the full state for every observable. The cost is that ``O(t)`` is now an object in operator space, with its own notion of "size" and its own entanglement structure.

### Vectorization

To put an operator inside an MPS we vectorize it. Writing ``\rho \mapsto |\rho\rangle\!\rangle`` for the column-stacking (or, here, operator-basis-component) isomorphism between the operator ``\rho`` acting on ``(\mathbb{C}^d)^{\otimes N}`` and a vector in the ``d^{2N}``-dimensional Liouville space, left/right multiplication becomes matrix multiplication on the doubled space:

```math
A\,\rho\,B \;\mapsto\; (A \otimes B^{\mathsf T})\,|\rho\rangle\!\rangle.
```

Under vectorization the Heisenberg and Lindblad equations of motion are *linear* in ``|\rho\rangle\!\rangle``, so the generator is an ordinary (super)operator that can be Trotterized into local gates exactly as in state-space TEBD.

### The local operator basis: generalized Gell-Mann, Pauli at `d = 2`

Rather than column-stacking, `MPSToolkit.jl` vectorizes each site in a fixed, cached, orthonormal
onsite basis — [`operator_basis_matrices`](@ref)`(d)`, the **generalized Gell-Mann basis** for a
local Hilbert space of dimension `d`. It is Hermitian, traceless apart from the identity
direction, and orthonormal under the Hilbert–Schmidt inner product
``\operatorname{tr}(P_\alpha^\dagger P_\beta) = \delta_{\alpha\beta}``. The ordering is fixed:
**the identity comes first**, as ``P_1 = I/\sqrt d``; then, for every pair ``j < k`` with
``1 \le j < k \le d``, the real symmetric and imaginary antisymmetric off-diagonal generators
built from ``|j\rangle\langle k|``; then the ``d - 1`` diagonal (Cartan) generators. At `d = 2`
this ordering reduces **exactly** — element-wise, not merely up to a change of basis — to

```math
(P_1, P_2, P_3, P_4) \;=\; (I, X, Y, Z)/\sqrt2 ,
```

the normalized Pauli basis every `pauli_*` helper uses; that identity is what makes `pauli_*`
the `d = 2` case of `operator_*` rather than merely a same-named alternative.
`operator_basis_matrices` is cached per `d`, so repeated calls — including the ones every
`operator_*` helper below makes internally — cost nothing after the first.

An ``N``-site operator becomes an MPS whose physical leg on each site is a dimension-``d^2``
index, and a product of local basis elements (e.g. the Pauli string ``X_1 Z_3`` at `d = 2`) is a
product state in this basis. This ordering is shared by every operator-space helper on this
page, so states and gates built by different helpers stay mutually consistent.

!!! note "Two normalizations"
    `pauli_components` decomposes a ``2\times2`` matrix against the **unnormalized** basis ``\{I, X, Y, Z\}``, whereas the operator-space *state* helpers use the **normalized** basis ``\sigma/\sqrt2`` (equivalently `operator_basis_matrices(2)`). The two differ by ``(\sqrt2)^N`` over ``N`` sites, so do not feed `pauli_components` output straight into `pauli_basis_state`.

### Operator entanglement and compressibility

The reason this representation is useful is that a vectorized operator has its own bipartite entanglement — the *operator entanglement* — defined through the singular values of ``|O\rangle\!\rangle`` across a cut. A local operator (a short Pauli string) starts with essentially zero operator entanglement. Under chaotic dynamics it spreads, and the operator entanglement generically grows, eventually saturating the available bond dimension. But for many physically important cases — transport of a conserved density, weakly dissipative open systems, short-time autocorrelators — the operator MPS stays compressible far longer than the corresponding Schrödinger-picture state, which is exactly the regime where operator-space TEBD, DMT, and DAOE pay off.

### Open-system Lindblad dynamics

For an open system the density matrix obeys the Lindblad master equation

```math
\dot\rho = \mathcal{L}\rho = -i[H, \rho] + \sum_j \Big( L_j\,\rho\,L_j^\dagger - \tfrac{1}{2}\{L_j^\dagger L_j,\, \rho\} \Big),
```

with Hamiltonian ``H`` and jump operators ``L_j``. Vectorized, ``\mathcal{L}`` is a (non-Hermitian) superoperator; the formal solution ``\rho(t) = e^{t\mathcal{L}}\rho(0)`` is again generated by a local object that can be Trotterized. `pauli_lindblad_generator` builds the dense local ``\mathcal{L}`` in the Pauli basis, and `pauli_gate_from_lindbladian` exponentiates it into a single-step gate ``e^{\,dt\,\mathcal{L}}`` ready for the TEBD scheduler. Setting all ``L_j = 0`` recovers the purely Hamiltonian (commutator) generator.

## Basis and state helpers

`operator_siteinds(N; d=2)` builds the `N` dimension-``d^2`` site indices that define operator
space; everything else is expressed over these indices. `pauli_siteinds(N)` is the `d = 2` case
(dimension-4 sites, basis `(I, X, Y, Z)`).

```julia
using MPSToolkit
using ITensors
using ITensorMPS

sites = operator_siteinds(4; d = 2)                       # == pauli_siteinds(4)
op    = operator_basis_state(sites, ["I", "X", "I", "Z"])  # the Pauli string X_2 Z_4, == pauli_basis_state
```

`operator_basis_state(sites, labels; coefficient=1.0)` builds a product-state operator MPS by
selecting **one basis element per site**. Each label is an integer index `1:d^2`, or `"I"` (any
`d`), or `"X"`/`"Y"`/`"Z"` when `d == 2`; `coefficient` is an overall scalar stored on the first
tensor. `pauli_basis_state` is the `d = 2` case, accepting the same letter labels.

Two builders exist for operators that are **not** a single basis element:

- `operator_product_state(sites, ops; coefficient=1.0)` — the tensor product ``\bigotimes_j O_j``
  of dense local matrices `ops`, one `d x d` matrix per site. This is the **literal** operator:
  its amplitude at basis label `alpha` on site `j` is `tr(P_alpha' * ops[j])`, so an identity
  argument contributes a factor `sqrt(d)` (since `I = sqrt(d) * P_1`).
- `operator_local_sum_state(sites, op, coeffs)` — the literal local-density sum ``\sum_j c_j O_j``
  (dense `op` on site `j` with weight `coeffs[j]`, the literal identity `I` elsewhere), at bond
  dimension 2. Uses the same literal convention as `operator_product_state`, so the two compose
  correctly, e.g. `add(operator_product_state(sites, fill(identity_matrix, N)),
  operator_local_sum_state(sites, O, coeffs))` is the literal `I + sum_j c_j O_j` — exactly the
  domain-wall construction transport runs use (see [DMT](dmt.md)).

!!! warning "`operator_basis_state` is normalized; `operator_product_state` is literal — they differ by `d^(N/2)`"
    `operator_basis_state(sites, fill(1, N))` places amplitude `1` on the all-identity label, so
    it represents the **normalized** operator ``(I/\sqrt d)^{\otimes N}`` — amplitude `1`, not
    the literal identity. `operator_product_state(sites, fill(identity_matrix, N))` (a dense
    `d x d` identity on every site) represents the **literal** identity ``I^{\otimes N}`` —
    amplitude ``\sqrt d^{\,N}``. The two differ by exactly `d^(N/2)`; verified numerically: at
    `d = 3, N = 5` the ratio of the two amplitudes is `15.588457... == sqrt(3)^5` to machine
    precision. Mixing the two conventions without rescaling (e.g. adding a normalized-basis state
    to a literal one) silently produces the wrong operator. `operator_local_sum_state` and
    `operator_product_state` share the literal convention, so they always compose correctly with
    each other, as in the domain-wall example above.

!!! warning "At `d >= 3`, a physical `S^z` is not a single basis element"
    `pauli_basis_state(sites, ["Z", ...])` works at `d = 2` because the physical spin operator is
    proportional to basis element 4 there. This stops holding at `d >= 3`: the physical
    `S^z = diag(J, J-1, ..., -J)` (`J = (d-1)/2`) is generally a **combination** of more than one
    diagonal Gell-Mann generator. At `d = 3`, `S^z = diag(1, 0, -1)` decomposes onto *two* of the
    9 `operator_basis_matrices(3)` directions, not one (`tr(P' * Sz)` is nonzero for exactly 2 of
    the 9 basis matrices `P`, verified numerically) — so no integer label selects it, and
    `operator_basis_state` cannot build it, since it only ever selects a single basis direction.
    Build the physical `S^z` — or any operator that is not literally a basis element — from its
    dense matrix with `operator_product_state` / `operator_local_sum_state` instead. This is why
    the [DMT higher-spin example](dmt.md) builds its domain wall that way rather than with a
    basis label.

`pauli_total_sz_state(sites; coefficient=nothing)` and `pauli_domain_wall_state(sites; kink,
coefficient=1.0)` are `d = 2`-only conveniences with their own pre-existing bond-dimension-2
normalization (see their docstrings for the exact factor); they require dimension-2 sites and
throw `ArgumentError` otherwise. For any other local dimension, build the equivalent state with
`operator_local_sum_state` directly, using the dense `S^z` matrix as in the warning above.

```julia
sites = pauli_siteinds(6)
Sz    = pauli_total_sz_state(sites)        # MPS for ∑_j S_j^z, bond dimension 2
```

The dense single-site primitives back these helpers and are useful for building or inspecting
custom local objects:

- `operator_basis_matrices(d)` — the cached generalized Gell-Mann basis, `d^2` dense `d x d`
  matrices, identity first (see **The local operator basis** above).
- `local_dimension(x)` — recover the physical dimension `d` from an operator-space `Index` (or
  its dimension `d^2`); every `operator_*` helper uses this internally, so `d` never has to be
  passed twice once the sites are built.
- `pauli_matrices(; include_identity=true)` — a named tuple of the dense ``2\times2`` Pauli matrices in `(I, X, Y, Z)` order (or `(X, Y, Z)` when the identity is dropped); the `d = 2` case of `operator_basis_matrices`, returned as a named tuple rather than a vector.
- `pauli_basis(; include_identity=true)` — the same matrices as a vector of `Symbol => Matrix` pairs.
- `pauli_components(operator; include_identity=true)` — the Hilbert–Schmidt coefficients of a ``2\times2`` `operator` against the **unnormalized** basis, normalized so that `operator == sum(c[name] * pauli_matrices()[name] for name in keys(c))`.

```julia
mats = pauli_matrices()                       # (I, X, Y, Z)
c    = pauli_components([0 1; 1 0])           # (I=0, X=1, Y=0, Z=0)  -> X
```

## Local Hamiltonian and Lindblad maps

These helpers turn a dense *physical-space* operator (acting on one or a few sites of local
dimension `d`) into the corresponding dense *operator-space* gate. The number of sites is
inferred from the matrix dimension, so a ``d^2 \times d^2`` matrix is read as a two-site gate;
sparse input (e.g. from `EDKit.spin(...; D = d)`) is accepted and densified internally.

- `operator_gate(op; d=2)` — the superoperator induced by conjugation,
  ``G[\alpha_{\text{out}},\alpha_{\text{in}}] = \operatorname{tr}\!\big[P_{\alpha_{\text{out}}}^\dagger\, A\, P_{\alpha_{\text{in}}} A^\dagger\big]``
  with the basis operators ``P_\alpha`` of `operator_basis_matrices(d)` (row = output index,
  column = input index). This is the operator-space image of ``\rho \mapsto A\rho A^\dagger``.
  `pauli_gate(unitary)` is the `d = 2` case.
- `operator_gate_from_hamiltonian(h, dt; d=2)` — the conjugation gate for the unitary
  ``e^{-i\,dt\,h}``, i.e. one Heisenberg/commutator TEBD step; the function to pass as
  `map_hamiltonian` to the TEBD scheduler. `pauli_gate_from_hamiltonian(h, dt)` is the `d = 2`
  case.
- `operator_gate_from_imaginary_time(h, dbeta; d=2)` — the two-sided imaginary-time slice
  ``\rho \mapsto e^{-(d\beta/2) h}\rho\, e^{-(d\beta/2) h}`` for Hermitian `h`, used by
  `operator_gibbs_state` below. `pauli_gate_from_imaginary_time(h, dbeta)` is the `d = 2` case.
- `operator_lindblad_generator(h, jumps; d=2)` — the dense local Lindbladian ``\mathcal{L}`` for
  Hamiltonian `h` and jump operator(s) `jumps` (a single matrix or a collection, each matching
  the dimension of `h`). `pauli_lindblad_generator(h, jumps)` is the `d = 2` case.
- `operator_gate_from_lindbladian(h, jumps, dt; d=2)` — the dissipative one-step gate
  ``e^{\,dt\,\mathcal{L}}``. `pauli_gate_from_lindbladian(h, jumps, dt)` is the `d = 2` case.

```julia
using MPSToolkit

# Two-site XXZ bond -> Heisenberg-picture conjugation gate (d = 2, operator_* and pauli_* agree)
h    = spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=1.0, Jz=0.8)
gate = operator_gate_from_hamiltonian(h, 0.05; d=2)   # == pauli_gate_from_hamiltonian(h, 0.05); 16 x 16 superoperator

# Single-site dephasing channel -> dissipative gate
sx   = pauli_matrices().X
sz   = pauli_matrices().Z
L    = pauli_gate_from_lindbladian(sx, [sqrt(0.3) * sz], 0.05)   # 4 x 4 gate
```

`operator_gate_from_hamiltonian` (and its `pauli_gate_from_hamiltonian` `d = 2` specialization)
matches the `(h, dt) -> gate` signature expected by `local_gates_from_hamiltonians` and
`tebd_strang_evolution`, which is what makes operator-space TEBD a drop-in use of the ordinary
scheduler, at any `d`. Building the dense local `h` for spin `> 1/2` is outside this package's
scope — [DMT's Higher spin section](dmt.md) shows a `d = 3` example sourcing `h` from
[EDKit.jl](https://github.com/jayren3996/EDKit.jl) (not an MPSToolkit dependency).

## Operator-space TEBD

Operator-space TEBD reuses the entire TEBD scheduler. The only change from state-space TEBD is the `map_hamiltonian` keyword: instead of the default ``h \mapsto e^{-i\,dt\,h}`` state gate, you pass `pauli_gate_from_hamiltonian`, which produces the conjugation superoperator. The scheduler (`tebd_strang_evolution` → `local_gates_from_hamiltonians` → `LocalGateEvolution`) is otherwise identical, so the odd–even–odd Strang schedule, the per-step truncation budget, and the [`evolve!`](@ref) driver all carry over unchanged.

The example below evolves a single-site ``Z`` operator under a transverse-field Ising chain in the Heisenberg picture and reports the growth of bond dimension and operator (bipartite) entanglement.

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 8
dt     = 0.05

# Operator-space site indices and an initial local operator: Z on the central site.
sites  = pauli_siteinds(nsites)
labels = fill("I", nsites)
labels[nsites ÷ 2] = "Z"
op     = pauli_basis_state(sites, labels)

# Strang TEBD over operator space: the ONLY operator-space-specific choice is
# map_hamiltonian = pauli_gate_from_hamiltonian, which turns each local Hamiltonian
# into its conjugation (Heisenberg) superoperator instead of a state propagator.
evolution = tebd_strang_evolution(
  nsites,
  dt;
  local_hamiltonian = (bond, weight) ->
    weight * spinhalf_tfim_bond_hamiltonian(nsites, bond; J=1.0, g=0.9),
  map_hamiltonian = pauli_gate_from_hamiltonian,
  maxdim = 64,
  cutoff = 1e-12,
)

# Each evolve! call runs one full Strang sweep (dt). Run several to spread the operator.
nsweeps = 40
for sweep in 1:nsweeps
  evolve!(op, evolution)
  if sweep % 10 == 0
    cut = nsites ÷ 2
    println("sweep $sweep  t = $(round(sweep * dt; digits=3))  ",
            "maxlinkdim = $(maxlinkdim(op))  ",
            "S_op(cut=$cut) = $(round(bond_entropy(op, cut); digits=4))")
  end
end
```

`bond_entropy` here is the von Neumann entropy of the *vectorized operator* across the cut — the operator entanglement — so its growth is a direct diagnostic of operator spreading and of how hard the evolution is to compress. Replacing `pauli_gate_from_hamiltonian` with `(h, dt) -> pauli_gate_from_lindbladian(h, jumps, dt)` (for fixed local `jumps`) turns the same scheduler into a Lindblad propagator; see the open-system notebooks below.

The same pattern works at any local dimension `d`: build `sites = operator_siteinds(nsites; d)` instead of `pauli_siteinds(nsites)`, and pass `map_hamiltonian = (h, dt) -> operator_gate_from_hamiltonian(h, dt; d)` instead of `pauli_gate_from_hamiltonian`. Nothing else in the scheduler changes.

## Vectorizing operators, thermal states, and expectations

Three families of helpers connect physical-space operators to the operator-basis
representation. Each is generic in `d` (`operator_*`), with a `pauli_*` name as the `d = 2`
case:

- **Vectorizers.** `operator_state_from_mpo(op, sites)` turns any local-dimension-`d` `MPO` into
  the operator MPS with amplitudes ``\mathrm{tr}(P_\alpha^\dagger O)`` (bond dimension
  preserved), and `operator_superoperator_mpo(op, sites)` lifts an MPO ``M`` to the
  operator-space MPO implementing the sandwich ``\rho \mapsto M \rho M^\dagger`` (bond dimension
  squared). `pauli_state_from_mpo` / `pauli_superoperator_mpo` are the `d = 2` cases. The PXP
  constraint objects `pauli_pxp_constraint_state` and `pauli_pxp_constraint_projector` remain
  **spin-1/2 only** (they raise `ArgumentError` at any other local dimension) — see
  [DMT](dmt.md) for the constrained-transport workflow.
- **Imaginary-time preparation.** `operator_gate_from_imaginary_time(h, dbeta; d=2)` builds the
  two-sided slice ``\rho \mapsto e^{-(d\beta/2)h} \rho\, e^{-(d\beta/2)h}``, and
  `operator_gibbs_state(sites, terms, weights; ...)` Trotterizes
  ``\rho \propto e^{-K/2} \rho_0 e^{-K/2}`` with ``K = \sum_j w_j h_j`` — uniform weights give
  Gibbs states, sign-split weights give the energy domain walls used in transport runs.
  `pauli_gate_from_imaginary_time` / `pauli_gibbs_state` are the `d = 2` cases.
- **Expectations.** `operator_trace(rho)`, `operator_expectation(rho, op, start)`, and the
  ``O(N)`` batched `operator_expectation_profile(rho, terms)` evaluate ``\mathrm{tr}(\rho O)`` —
  optionally over ``\mathrm{tr}(\rho)`` (`normalize=true`, the default) — for dense local
  windows directly from the operator MPS with cumulative identity environments. `pauli_trace`,
  `pauli_expectation`, `pauli_expectation_profile` are the `d = 2` cases.

!!! note "Trace convention: `d^(N/2)`, and where it overflows"
    The unnormalized trace of an `N`-site operator carries a prefactor `sqrt(d)^N` — `2^(N/2)`
    at `d = 2`, generalizing to `d^(N/2)` at any `d` (verified numerically: `operator_trace` of
    the literal `N = 5`, `d = 3` identity is `243.0 == 3.0^5`, and of the *normalized*-basis
    identity state — `operator_basis_state(sites, fill(1, N))` — is
    `15.588... == sqrt(3.0)^5`). This prefactor overflows `Float64` around
    `N ~ 2048 / log2(d)` — far beyond any feasible operator-space MPS size, but worth knowing
    before reading an absolute (`normalize=false`) trace or expectation value on a very long,
    low-`d` chain. Normalized ratio observables (`normalize=true`, the default) are unaffected,
    since the prefactor cancels between numerator and denominator.

!!! note "When `normalize=true` refuses"
    Dividing by the trace is meaningless when the trace is a cancellation residue rather than a
    trace, so `normalize=true` raises an `ArgumentError` when
    ``|\mathrm{tr}(\rho)| \le \sqrt{\varepsilon}\,\|\rho\|_{HS}`` — the traceless case (a
    transport current, an evolved two-point correlator). Measure those with `normalize=false`.
    The comparison is between the *physical* trace and the Hilbert-Schmidt norm, both evaluated
    in log space so neither has to be representable: a positive `\rho` always satisfies
    ``\mathrm{tr}(\rho) = \sum_i \lambda_i \ge \sqrt{\sum_i \lambda_i^2} = \|\rho\|_{HS}``, so
    **no thermal state is ever rejected**, at any temperature or chain length. (Before
    2026-08-10 the check compared the identity *coefficient* ``\mathrm{tr}(\rho)/d^{N/2}``
    against the same right-hand side, which carried an implicit ``d^{N/2}`` and therefore
    rejected the positive product state ``\rho = \bigotimes_j e^{-\beta Z_j}`` at
    ``\beta = 0.5, N = 240``.)

## Related operator-space tools

This page covers the generic operator-basis representation (any local dimension `d`, Pauli at
`d = 2`) and operator-space TEBD. Two specialized truncation backends live on their own pages
and reuse the same states and gates:

- **[DMT](dmt.md)** (density matrix truncation) — a *transport-specific* truncation that
  protects the identity/trace component and short-range reduced-operator data at every bond
  before discarding long-range correlations. Generic in `d` (see DMT's own **Higher spin**
  section). Documented there: `dmt_step!`, `dmt_evolve!`, `DMTGateEvolution`, and `DMTOptions`.
- **[DAOE](daoe.md)** (dissipation-assisted operator evolution) — diagonal projector MPOs that
  damp operator components by operator *size* rather than bond dimension. **Spin-1/2 only**:
  `pauli_daoe_projector` and `fdaoe_projector` require dimension-4 Pauli sites and raise
  `ArgumentError` at any other local dimension. Documented there.

For general operator-space evolution without a transport assumption, prefer plain
[`LocalGateEvolution`](@ref) TEBD with the helpers on this page (`operator_*` at any `d`,
`pauli_*` at `d = 2`); DMT's identity-preserving truncation rule is only appropriate when the
problem has transport structure.

## API

```@docs
operator_siteinds
operator_basis_state
operator_product_state
operator_local_sum_state
pauli_siteinds
pauli_basis_state
pauli_total_sz_state
pauli_domain_wall_state
operator_basis_matrices
local_dimension
pauli_matrices
pauli_basis
pauli_components
operator_gate
operator_gate_from_hamiltonian
operator_lindblad_generator
operator_gate_from_lindbladian
pauli_gate
pauli_gate_from_hamiltonian
pauli_lindblad_generator
pauli_gate_from_lindbladian
operator_state_from_mpo
operator_superoperator_mpo
pauli_state_from_mpo
pauli_superoperator_mpo
operator_gate_from_imaginary_time
operator_gibbs_state
pauli_gate_from_imaginary_time
pauli_gibbs_state
operator_trace
operator_expectation
operator_expectation_profile
pauli_trace
pauli_expectation
pauli_expectation_profile
```

## Examples

- [examples/operator_space/operator_tebd_helper_apis.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/operator_tebd_helper_apis.ipynb)
- [examples/operator_space/dmt_scheduler.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/dmt_scheduler.ipynb)
- [examples/open_systems/pauli_lindblad_tebd.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/open_systems/pauli_lindblad_tebd.ipynb)
- [examples/open_systems/boundary_driven_xxz_steady_state.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/open_systems/boundary_driven_xxz_steady_state.ipynb)

## References

- Tibor Prosen and Marko Žnidarič, [Matrix product simulations of non-equilibrium steady states of quantum spin chains](https://doi.org/10.1088/1742-5468/2009/02/P02035), J. Stat. Mech. (2009) P02035 — vectorized-MPS treatment of Lindblad dynamics and boundary-driven steady states.
- Tomaž Prosen and Iztok Pižorn, [Operator space entanglement entropy in a transverse Ising chain](https://doi.org/10.1103/PhysRevA.76.032316), Phys. Rev. A 76, 032316 (2007) — operator entanglement and its growth.
- Tibor Rakovszky, C. W. von Keyserlingk, and Frank Pollmann, [Dissipation-assisted operator evolution method for capturing hydrodynamic transport](https://arxiv.org/abs/2004.05177) — operator spreading and size-based truncation.
- Christopher David White, Michael Zaletel, Roger S. K. Mong, and Gil Refael, [Quantum dynamics of thermalizing systems](https://doi.org/10.1103/PhysRevB.97.035127), Phys. Rev. B 97, 035127 (2018) — density matrix truncation in operator space.
