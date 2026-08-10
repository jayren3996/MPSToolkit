# DMT

*Density matrix truncation: a transport-specific operator-space compression that shares TEBD's gate scheduling but replaces the local SVD truncation with one that protects the identity and low-weight operator components at every bond.*

!!! warning "Transport-specific algorithm"
    DMT is a **transport-specific** truncation scheme.  Its core truncation step protects the
    identity component of the vectorized operator at every bond, which is the right bias for
    transport observables (e.g. spin or energy diffusion) but **not** for general
    operator-space tasks.  If your problem does not have a transport structure, use ordinary
    TEBD truncation ([`LocalGateEvolution`](@ref)) in operator space instead.

## Background

DMT (density matrix truncation) is a method for simulating the real-time dynamics of
*local densities* in one-dimensional systems. It was introduced by White, Zaletel, Mong, and
Refael to follow thermalization and hydrodynamic transport in regimes where ordinary matrix
product state (or matrix product operator) truncation fails. The setting is operator space:
the object being evolved is not a wavefunction but a density matrix ``\rho`` (or, more
generally, a mixed operator) written as a matrix product state over a *vectorized* local
basis. MPSToolkit vectorizes onto a normalized, Hermitian, orthonormal onsite basis (a
generalized Gell-Mann basis, identity first) for a local Hilbert space of any dimension ``d``;
at ``d = 2`` this basis is *exactly* the normalized Pauli basis ``(I, X, Y, Z)/\sqrt2``, so
every `pauli_*` helper is the ``d = 2`` case of the generic `operator_*` one. See
[Operator Space](operator-space.md) for the vectorization conventions, `operator_siteinds` /
`pauli_siteinds`, `operator_basis_state` / `pauli_basis_state`, and the gate builders. This
page's worked examples use the spin-1/2 `pauli_*` names throughout, except the dedicated
**Higher spin** section below — every DMT entry point (`DMTOptions`, `DMTGateEvolution`,
`dmt_step!`, `dmt_evolve!`) works unchanged at any local dimension.

### Why ordinary truncation fails for transport

Real-time evolution of a generic operator spreads support over an exponentially growing set of
Pauli strings. The information that matters for *transport* — how a conserved charge such as
``S^z`` or energy flows along the chain — lives in the slowly varying, low-weight part of the
operator: the identity component on most sites, single-site densities, and short-range
correlations. The fast-growing, high-weight tail carries little of this hydrodynamic signal
and is, for transport purposes, effectively scrambled.

Standard SVD truncation keeps the largest Schmidt values at each bond. Operationally that is
the right thing to do for a pure state, where the Schmidt spectrum measures entanglement and
the discarded weight bounds the error in the *state*. In operator space it is the wrong bias:
the largest singular vectors are dominated by the proliferating high-weight strings, so an SVD
truncation preferentially *keeps* the scrambled tail and *discards* the low-weight components
that encode the conserved density. The result is that diffusion constants and other transport
coefficients come out wrong, even when the nominal truncation error looks small.

### The truncation rule

DMT changes which information is protected, and it does so *exactly* rather than by degree. At
a bond with the orthogonality center placed on it, write the bond matrix (the Schmidt/bond
tensor connecting the left and right halves of the operator MPS) as ``M``. DMT expresses ``M``
in a basis whose leading directions span the local-operator subspace on the ``n`` sites
adjacent to the cut on each side, with ``n = (\texttt{preserve\_diameter} - 1)/2`` (``n = 1``,
i.e. the first ``d^2`` rows and the first ``d^2`` columns, at the default
`preserve_diameter = 3`). Those ``d^{2n}`` leading rows *and* the ``d^{2n}`` leading columns of
``M`` — together with a rank-one trace/identity connector ``C`` that lives inside the same
protected block — are reinstated **exactly** after truncation; only the doubly-orthogonal
complement that lies outside both the row and the column protected block is compressed by SVD:

```math
M \;=\; \underbrace{C}_{\text{rank-1 trace connector}} \;+\; \underbrace{\text{(protected rows and columns)}}_{\text{kept exactly}} \;+\; \underbrace{D}_{\text{doubly-orthogonal complement: truncated}} .
```

Because the protected block is reinstated exactly — not merely weighted more heavily by the
SVD — **every observable of diameter at most `preserve_diameter` survives the truncation
exactly** (to floating-point precision), independent of how the complement ``D`` is compressed.
This exact-preservation property, for the leading ``d^2`` rows and columns at diameter 3, is the
actual content of the guarantee DMT is defined by (White, Zaletel, Mong, Refael,
[arXiv:1707.01506](https://arxiv.org/abs/1707.01506), Sec. III). The resulting bond has rank at
most

```math
\operatorname{rank}(M') \;\le\; 2\, d^{2n} + \chi' ,
```

where ``\chi'`` is the number of singular directions kept from the complement — the
``\chi = \chi_{\text{preserve}} + \chi_{\text{extra}}`` split of Ye, Machado, Mong, Yao
([arXiv:1902.01859](https://arxiv.org/abs/1902.01859)). `maxdim` is this total rank bound
(see **Budget semantics** below), so `preserve_diameter` and `maxdim` together fix
``\chi' = \texttt{maxdim} - 2 d^{2n}``. Locally conserved densities and their slow hydrodynamic
tails survive truncation exactly even though they may correspond to small singular values; what
gets compressed is only the high-weight, long-range connected correlations that do not feed
back into transport.

The trace component is treated with care numerically. A traceless operator (for example a
transport *current*) has essentially zero identity overlap, so the connector subtraction is
applied only when the trace component is well-conditioned relative to the operator norm;
otherwise it is skipped while the rest of the block is still truncated, so the bond dimension
target is always enforced. This is implemented in the DMT truncation kernel rather than being
something the user tunes directly.

### Forward and reverse sweeps

DMT is built on the same gate decomposition as TEBD: the generator of the dynamics is split
into local (typically nearest-neighbor) terms, each exponentiated into a gate, and the gates
are applied in a scheduled sweep along the chain. A single DMT update applies one local gate
and then performs the identity-preserving truncation on the bonds inside that gate's support.
Sweeping left-to-right and then right-to-left over a step makes the truncation bias symmetric
in space and keeps the orthogonality center moving with the active bond. In MPSToolkit a sweep
direction of `:R` is used for the forward pass and `:L` for the reverse pass, and one
evolution call runs a complete forward-plus-reverse sweep (repeated `nstep` times).

## DMT in MPSToolkit

DMT is exposed as an operator-space truncation *backend* that plugs into the same scheduling
abstraction as TEBD. The division of labor is:

- **The scheduler decides where gates act.** A forward `schedule` and a `reverse_schedule` of
  bond positions, together with a gate (or per-step gate provider), describe the circuit. This
  is exactly the role the schedule plays for [`LocalGateEvolution`](@ref).
- **The DMT kernel decides how each post-gate truncation is performed.** Instead of an SVD that
  keeps the largest singular values, it applies the identity-preserving rule described above.

Four public symbols cover the workflow:

- `DMTGateEvolution` bundles the gate specification, the forward and reverse schedules, the
  step count, and the truncation budgets into one configuration object. It is the DMT analogue
  of `LocalGateEvolution`.
- `DMTOptions` collects just the truncation settings (`maxdim`, `cutoff`, `gate_maxdim`,
  `preserve_diameter`, `truncation`) so they can be passed to a single step without
  constructing a full evolution.
- `dmt_step!` is the lowest-level entry point: apply one local gate at a bond and then perform
  the associated DMT truncation. It exists both in a keyword form and in an overload that takes
  a `DMTOptions`.
- `dmt_evolve!` is the scheduled driver. It walks the forward schedule applying `dmt_step!`
  with `direction=:R`, then walks the reverse schedule with `direction=:L`, repeats that for
  `evo.nstep` sweeps, and finally normalizes the operator MPS.

`dmt_step!` first applies the raw gate — by default with **no cap at all**
(`gate_maxdim = 0` means "apply the gate exactly") and *no* cutoff — then truncates the
(possibly much wider) bond back down to `maxdim` using the DMT rule. In this two-stage
structure the gate is allowed to expand the bond freely and the transport-aware kernel then
compresses it; a positive `gate_maxdim` pre-truncates the gate-inflated bond with a plain SVD
*before* DMT ever sees it, discarding small singular values before the kernel can decide which
of them are protected — exactly the error DMT exists to avoid — so raise it only to bound peak
memory, never for accuracy.

The truncation budget has a few knobs worth understanding:

- `maxdim`: the **total** post-truncation bond dimension, *inclusive* of the protected block
  (see **Budget semantics** below) — not the complement alone, and not an addition on top of a
  separately-sized buffer.
- `cutoff`: cutoff used only in the final "repair" SVD that re-factorizes the truncated bond
  back into MPS form.
- `gate_maxdim`: temporary bond-dimension ceiling allowed while the raw gate is applied, before
  DMT compresses the bond back to `maxdim`. Defaults to `0`, meaning no cap: the gate is applied
  exactly.
- `preserve_diameter`: the positive odd diameter of the observables preserved exactly (default
  `3`). `radius = (preserve_diameter - 1) / 2` sites are protected on each side of the cut, and
  the protected block has dimension `2 d^(2 radius)`. Replaces the removed `connector_buffer`
  (see the migration note below).
- `truncation`: `:dense` (default) or `:random` complement truncation. `:random` is faster at
  large bond dimension and preserves the guarantee to the same tolerance, but is not
  deterministic — see [`DMTOptions`](@ref) for the measured tradeoff and why `:dense` ships as
  the default. `:dense` also sets peak memory, not just speed: it materializes the `chi x chi`
  complement and factorizes it, where `chi` is the **gate-inflated** bond `d^2 maxdim`, for a
  measured transient of ~6.4 `chi x chi` `ComplexF64` matrices against ~2.2 for `:random` —
  0.33 GB at `d = 3, maxdim = 200`, 1.1 GB at `d = 4, maxdim = 200`, and ~7 GB at
  `d = 4, preserve_diameter = 5, maxdim = 513` — so at `d >= 4`, or at `preserve_diameter = 5`,
  choosing `:random` is a memory decision and not only a speed one.

!!! warning "Migration: `connector_buffer` removed"
    `connector_buffer` no longer exists in `DMTOptions`, `DMTGateEvolution`, or `dmt_step!`;
    passing it raises an `ArgumentError` explaining the change rather than a bare `MethodError`.
    Replace it with `preserve_diameter` (odd, default `3`). The other half of the migration is
    a change in what `maxdim` *means*: it is now the **total** bond dimension, inclusive of the
    protected block, rather than a target added on top of a separately-sized connector buffer.
    A script re-run at the same `maxdim` therefore keeps `maxdim - (2 d^(2n) - 1)` complement
    directions (`n = (preserve_diameter - 1) / 2`), not `maxdim` complement directions plus a
    buffer on top.

`maxdim` must be at least `2 d^(2n) + 1`; an `ArgumentError` naming `d`, `preserve_diameter`,
and the implied floor is raised at the start of the step or sweep, before anything is mutated.

| `preserve_diameter` | protected block `2 d^(2n)` | minimum `maxdim`, `d = 2` | `d = 3` | `d = 4` |
|:--|--:|--:|--:|--:|
| 3 | `2 d^2`  | 9  | 19  | 33  |
| 5 | `2 d^4`  | 33 | 163 | 513 |

The bottom-right cell is a memory decision as well as a budget: at `d = 4, preserve_diameter = 5`
the floor `maxdim = 513` gives a gate-inflated bond of `d^2 maxdim = 8208`, and the default
`truncation = :dense` transiently needs ~7 GB there (see the `truncation` bullet above), against
~2.4 GB for `:random`. Pass `truncation = :random` at that corner of the table unless the machine
has the headroom.

DMT works at any local Hilbert space dimension `d`: the step validates only that every site in
the gate's window shares a common local dimension (`operator_siteinds(nsites; d)` sites for any
`d >= 2`, including the `d = 2` `pauli_siteinds` case). Periodic-boundary local gates are not
supported.

## Worked example

The following is fully self-contained. It prepares a single-site ``Z`` operator on a six-site
Pauli chain, builds a nearest-neighbor TFIM bond gate, lays out a left-to-right forward
schedule (with the default reversal for the backward sweep), and runs one DMT sweep.

```julia
using MPSToolkit
using ITensors
using ITensorMPS

nsites = 6
sites = pauli_siteinds(nsites)
state = pauli_basis_state(sites, ["Z", "I", "I", "I", "I", "I"])

gate = pauli_gate_from_hamiltonian(spinhalf_tfim_bond_hamiltonian(nsites, 1; J=1.0, g=0.6), 0.05)
schedule = collect(1:(nsites - 1))

evolution = DMTGateEvolution(
  fill(gate, length(schedule)),
  0.05;
  schedule=schedule,
  reverse_schedule=collect(reverse(schedule)),
  maxdim=16,
  cutoff=1e-10,
)

dmt_evolve!(state, evolution)
```

What to expect:

- The gate is a ``16 \times 16`` superoperator: the TFIM bond Hamiltonian acts on two spin-1/2
  sites, so its induced Pauli-space gate spans two dimension-``d^2 = 4`` operator sites.
  `dmt_evolve!` infers this span automatically from the gate size and the site dimensions — at
  any `d`, not just `d = 2`.
- `schedule = 1:(nsites-1)` applies the two-site gate on every adjacent bond in the forward
  pass; `reverse_schedule = reverse(schedule)` undoes that ordering for the backward pass. With
  this default reversal the driver maps reverse-sweep gates back to their forward counterparts
  automatically, so passing one gate per bond is enough.
- After the forward and reverse sweeps the operator MPS is normalized in place and returned.
  Because the initial operator is a single low-weight Pauli string, the bond dimension grows
  only as the gate spreads support, and DMT holds the identity/low-weight content while
  trimming the bond back toward `maxdim = 16`.

To run a single update by hand instead of a full sweep, call `dmt_step!` directly — for example
`dmt_step!(state, gate, 1; maxdim=16, cutoff=1e-10, direction=:R)`, or bundle the settings with
`dmt_step!(state, gate, 1, DMTOptions(maxdim=16, cutoff=1e-10); direction=:R)`. The scheduled
driver is simply this step repeated over the forward and reverse schedules.

## Worked example: transport exponents in the XXZ chain

The headline application of DMT is reading off a transport coefficient or dynamical exponent. The
spin-1/2 XXZ chain

```math
H = \sum_j \left( S^x_j S^x_{j+1} + S^y_j S^y_{j+1} + \Delta\, S^z_j S^z_{j+1} \right)
```

has three infinite-temperature spin-transport regimes set by the anisotropy ``\Delta``,
distinguished by the dynamical exponent ``z`` in ``\langle x^2 \rangle(t) \sim t^{2/z}``:

| ``\Delta``     | regime               | ``z``   |
|:---------------|:---------------------|:--------|
| ``\Delta < 1`` | ballistic            | ``1``   |
| ``\Delta = 1`` | superdiffusive (KPZ) | ``3/2`` |
| ``\Delta > 1`` | diffusive            | ``2``   |

Heisenberg-evolve the local operator ``O(0) = \sigma^z`` at the chain center in Pauli operator
space with the DMT backend, and read off the infinite-temperature autocorrelation profile
``p(x, t) = \langle \sigma^z_x \mid O(t) \rangle`` — the coefficient of the single-``Z`` string on
each site. Its second moment ``M_2(t) = \sum_x (x - x_0)^2\, p(x, t) / \sum_x p(x, t) \sim t^{2/z}``
gives ``z``. Forming ``M_2`` as a ratio makes it invariant to the `normalize!` rescaling that
`dmt_evolve!` applies, and the profile is computed in one ``O(N)`` sweep with cumulative identity
(trace) environments — the same contraction DMT uses internally — rather than ``N`` separate inner
products:

```julia
using MPSToolkit, ITensors, ITensorMPS

nsites, maxdim, dt, ncall = 30, 24, 0.1, 34   # one evolve! call advances t by 2*dt
center = (nsites + 1) ÷ 2

# coefficient of the single-Z string on every site, in one O(N) sweep
function single_z_profile(state)
    n = length(state)
    cap(s, k) = (t = ITensor(s); t[s => k] = 1.0; t)   # k = 1 -> identity, k = 4 -> single Z
    right = Vector{ITensor}(undef, n + 1); right[n + 1] = ITensor(1.0)
    for s in n:-1:1
        right[s] = right[s + 1] * cap(siteind(state, s), 1) * state[s]
    end
    left = ITensor(1.0); p = Vector{Float64}(undef, n)
    for x in 1:n
        p[x] = real(scalar(left * (cap(siteind(state, x), 4) * state[x]) * right[x + 1]))
        left = left * cap(siteind(state, x), 1) * state[x]
    end
    return p
end

function transport_trace(Delta)
    sites = pauli_siteinds(nsites)
    O = pauli_basis_state(sites, [j == center ? "Z" : "I" for j in 1:nsites])
    gate = pauli_gate_from_hamiltonian(
        spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=1.0, Jz=Delta), dt)
    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule, reverse_schedule=reverse(schedule),
        maxdim, cutoff=1e-10, gate_maxdim=4 * maxdim)   # == the default (0, no cap) at d=2, since d^2=4
    width2(state) = (p = single_z_profile(state);
        sum((x - center)^2 * p[x] for x in 1:nsites) / sum(p))
    times = [0.0]; M2 = [width2(O)]
    for _ in 1:ncall
        evolve!(O, evo); push!(times, times[end] + 2dt); push!(M2, width2(O))
    end
    return times, M2
end
```

Extract ``z`` by fitting a log-log slope of ``M_2(t)`` — but fit an **intermediate** time window,
not the whole trace:

```julia
function fit_z(times, M2, (tmin, tmax))
    idx = findall(i -> tmin <= times[i] <= tmax && M2[i] > 0, eachindex(times))
    lt, lm = log.(times[idx]), log.(M2[idx])
    n, sx, sy = length(lt), sum(lt), sum(lm)
    slope = (n * sum(lt .* lm) - sx * sy) / (n * sum(lt .^ 2) - sx^2)
    return 2 / slope
end

for Delta in (0.5, 1.0, 2.0)
    times, M2 = transport_trace(Delta)
    println(Delta, "  ->  z = ", round(fit_z(times, M2, (3.0, 6.5)); digits=3))
end
# 0.5  ->  z = 1.093   (ballistic,            z = 1)
# 1.0  ->  z = 1.349   (superdiffusive / KPZ, z = 3/2)
# 2.0  ->  z = 1.658   (diffusive,            z = 2)
```

These values are freshly measured, exactly as configured above (`nsites=30, maxdim=24`,
~2 minutes total), against the current, faithful DMT kernel. This is a small demonstration
configuration — a `maxdim = 24` bond at `d = 2` leaves only `chi' = maxdim - 8 = 16` complement
directions after the protected block — not a converged production run; see the tip below on
sharpening it.

The three values above sit a little further from `1, 3/2, 2` than the old `1.09/1.45/1.90` this
same snippet used to report, at the identical `nsites=30, maxdim=24` — not because the kernel got
less accurate, but because the fixed `maxdim` budget is now spent differently. This snippet used
to pass `connector_buffer = 4`, so the old kernel's complement got `maxdim - connector_buffer =
20` directions. The new kernel's protected block is not a tunable buffer but the structural
`2 d^2 = 8` directions the guarantee requires, leaving `maxdim - 8 = 16` for the complement — four
fewer than before, at the same nominal `maxdim`. That is the price of the guarantee the old
kernel only advertised: its `connector_buffer` mechanism measured diameter-`<=3` preservation
errors of order `0.3-0.5` (it reinstated the protected corner and then let its own repair SVD clip
part of it away), against `~1e-14` for this kernel on the same class of probe.

!!! tip "Fit the intermediate plateau, not the tail"
    A single power-law fit over the whole trace is biased two ways. At **early** times every
    ``\Delta`` spreads at the Lieb-Robinson speed, so the effective exponent looks ballistic
    (``z \to 1``). At **late** times the conserved-charge tail outgrows the kept bond dimension
    `maxdim`, DMT truncates it, and ``M_2`` growth is artificially suppressed — the effective
    exponent stops being a clean power law and breaks down. Between them is a plateau where the
    hydrodynamic scaling is visible and the bond dimension is still adequate; fitting there (here
    ``t \in [3, 6.5]``) recovers all three exponents from one short run. Extrapolating the
    corrupted late-time tail would push ``z`` the *wrong* way. Widen the window and raise
    `nsites`/`maxdim` to sharpen the values further, at higher cost.

A complete, runnable transport script — using the closely related **domain-wall melting** protocol
(the charge ``\mathcal{T}(t)\sim t^{1/z}`` transferred across the wall in place of the
autocorrelation width, with the same front-contamination guard) — is
[`examples/dmt/domain_wall_melting.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/domain_wall_melting.jl).
At production scale (`nsites=80, maxdim=48`, 100 forward+reverse sweeps) that script measures
``z = 1.021, 1.618, 1.792`` for the same three ``\Delta``, with total ``S^z`` drift
``\sim 5\text{-}7\times10^{-11}`` — closer to the ballistic/KPZ/diffusive targets ``1, 3/2, 2``
than the small demonstration above, as expected from the much larger chain and bond dimension.
This is a genuinely **different measurement** — a different quench protocol (domain-wall melt
versus autocorrelation width) at a different `(nsites, maxdim)` — not a recalculation of the
snippet above at higher precision, so do not read the two numbers side by side as before/after.

## Worked example: constrained energy transport in the PXP chain

The PXP chain — the effective model of a Rydberg-blockaded atom array —

```math
H = \Omega \left[ X_1 P_2 + \sum_{j=2}^{N-1} P_{j-1} X_j P_{j+1} + P_{N-1} X_N \right],
\qquad P = |0\rangle\langle 0| ,
```

flips a site only when both neighbors are in the ground state. Every term therefore commutes
with the global constraint projector ``P_G = \prod_j (1 - n_j n_{j+1})``: dynamics never
connects the constrained sector (no adjacent excitations) to blockade-violating
configurations. Energy is the only local conserved density, and its high-temperature
transport shows a slow crossover: bare PXP sits near a proximate integrable point and is
*near-ballistic* (``z \approx 1``) out to ``t \sim 100``, reaching KPZ superdiffusion (``z = 3/2``)
only as an intermediate-time transient, with ordinary *diffusion* (``z = 2``) the believed
ultimate asymptote (Ljubotina, Desaules, Serbyn, Papić, PRX **13**, 011033 (2023); McRoberts &
Moessner, PRL **133**, 256301 (2024)).

Simulating this with DMT needs three constraint-aware ingredients beyond the XXZ workflow
above, all built from the model helpers `pxp_term_hamiltonian` / `pxp_term_support` and the
constraint MPO `pxp_constraint_mpo`:

1. **The right infinite-temperature state.** The reference state of the constrained sector is
   ``\rho_\infty \propto P_G``, *not* the identity: the identity carries weight on
   blockade-violating configurations that the PXP Hamiltonian never mixes with physical
   dynamics but that would pollute traces and energy densities.
   `pauli_pxp_constraint_state(sites)` returns the vectorized ``P_G`` as a bond-dimension-2
   operator MPS (built via the generic vectorizer `pauli_state_from_mpo`).
2. **A sector-respecting initial state.** For a domain-wall quench,
   ``\rho(0) \propto e^{-\beta H_L} \otimes e^{+\beta H_R}`` (left half slightly cold, right
   half slightly hot) is prepared by imaginary-time TEBD in operator space:
   `pauli_gibbs_state(sites, terms, weights; initial_state=pauli_pxp_constraint_state(sites))`
   Trotterizes ``\rho \propto e^{-K/2}\, \rho_0\, e^{-K/2}`` with ``K = \sum_j w_j h_j``.
   The wall is encoded entirely in the `weights`: ``+\beta`` for terms supported in the left
   half, ``-\beta`` in the right half, and ``0`` for the two wall-straddling terms — with that
   choice no term couples the halves and the two exponentials factorize *exactly*. (Assigning
   straddling terms by their center site instead differs only by an ``O(\beta)`` local
   distortion at the wall and gives the same hydrodynamics.) Because each ``h_j`` commutes
   with ``P_G``, the preparation stays in the constrained sector up to truncation error.
3. **Constraint checkpoints during evolution.** Exact gates commute with ``P_G``, but DMT (or
   any) truncation leaks weight out of the sector. The operator-space projector MPO
   ``\rho \mapsto P_G \rho P_G`` — `pauli_pxp_constraint_projector(sites)`, bond dimension 4,
   built by the generic `pauli_superoperator_mpo` — removes the leaked weight, and
   `constrained_dmt_evolve!` interleaves it with the DMT sweeps.

The energy domain-wall melt below is the **only** DMT protocol; it ships as a full script at
[`examples/dmt/pxp_energy_melting.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/pxp_energy_melting.jl).

!!! warning "DMT is for density operators, not two-point correlators"
    DMT truncation protects the **identity/trace component** and nearby local Pauli data of a
    near-infinite-temperature *density operator* (see [`DMTOptions`](@ref)) — that is the regime it
    was designed and validated for, and the domain-wall melt below, where ``\rho(t)`` is exactly
    such a state, is its intended use. Do **not** use DMT (`dmt_evolve!` / `constrained_dmt_evolve!`)
    to Heisenberg-evolve a *traceless* operator such as a two-point energy correlator
    ``O(0) \propto P_G h_c P_G``: a traceless operator has no trace component for DMT to anchor on,
    so the truncation runs outside its design regime. Evolve correlation / Green's functions with
    ordinary TEBD ([`LocalGateEvolution`](@ref) / `tebd_evolve!`) instead.

### Energy domain-wall melt — the DMT protocol

The quench protocol: prepare ``\rho(0) \propto e^{-\beta H_L} \otimes e^{+\beta H_R}`` and
track the transferred energy. From
[`examples/dmt/pxp_energy_melting.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/pxp_energy_melting.jl):

```julia
using MPSToolkit, ITensors, ITensorMPS

nsites, wall, beta, dt = 40, 20, 0.25, 0.1
sites = pauli_siteinds(nsites)
terms = [(first(pxp_term_support(nsites, j)), pxp_term_hamiltonian(nsites, j)) for j in 1:nsites]
weight(j) = (s = pxp_term_support(nsites, j);
             last(s) <= wall ? beta : first(s) > wall ? -beta : 0.0)

rho = pauli_gibbs_state(sites, terms, [weight(j) for j in 1:nsites];
    nsteps=6, maxdim=64, initial_state=pauli_pxp_constraint_state(sites))

gates = [pauli_gate_from_hamiltonian(h, dt) for (_, h) in terms]
schedule = [start for (start, _) in terms]     # 3-site bulk gates, 2-site edge gates
evo = DMTGateEvolution(gates, dt; schedule, reverse_schedule=reverse(schedule),
    nstep=1, maxdim=32, cutoff=1e-12, gate_maxdim=128)
projector = pauli_pxp_constraint_projector(sites)

profile(state) = real.(pauli_expectation_profile(state, terms))   # e_j = tr(rho h_j)/tr(rho)
e0 = profile(rho)
dE = Float64[]
for k in 1:50                                   # one call advances t by 2*dt
    constrained_dmt_evolve!(rho, evo, projector; project_every=1)
    push!(dE, sum(profile(rho)[(wall + 1):end] - e0[(wall + 1):end]))
end
```

The transferred energy ``\Delta E(t) = \sum_{j > \mathrm{wall}} [e_j(t) - e_j(0)] \sim
t^{1/z}`` is fitted on a late window, with the same caveats as the XXZ example plus one more:
the domain-wall quench approaches its hydrodynamic power law *slowly* (quadratic start,
ballistic-looking middle), so short runs measure an **effective, crossover-bound exponent**
that overestimates ``1/z`` — the example prints the local log-log slope so the drift toward
the asymptote is visible. The projector application at a checkpoint compresses back to
`projector_maxdim` (default `2 * maxdim`) immediately: the projected state is only an
``O(\text{leakage})`` perturbation of the state DMT just compressed, so this plain-SVD step
perturbs the protected components at the same negligible order, while letting the bond stay
inflated until the next sweep adds substantial cost per sweep (measured 1.3–2× at moderate
``\chi`` and growing with ``\chi``) for no measured accuracy gain.

One further lesson from the production-scale validation runs: **profile mirror symmetry is
a useful truncation gauge**. The melt is a mirror-symmetric setup, so the profile
asymmetry should vanish in exact arithmetic; in practice it grows with truncation pressure
and *shrinks* with increasing ``\chi`` (it is truncation noise, not a sweep-direction bias).
Monitoring it alongside the conserved total gives a cheap convergence diagnostic.

In the melt the measurement is `pauli_expectation_profile(rho, terms)`, which
evaluates every ``\mathrm{tr}(\rho h_j)`` (over ``\mathrm{tr}(\rho)`` when normalized) in
one ``O(N)`` pass with cumulative identity environments, mixed window sizes included. The full script adds the diagnostics discussed above — the local-slope crossover column, the
conserved-total drift, and the constraint-leakage probe (residual leakage after projection vs.
leakage accrued by one unprojected sweep):
[`pxp_energy_melting.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/pxp_energy_melting.jl).

!!! tip "The checkpoint pattern is generic"
    Nothing in `constrained_dmt_evolve!` is PXP-specific: it interleaves DMT sweeps with any
    operator-space MPO. Other kinetically constrained models need only their own constraint
    MPO — build it in physical space (cf. `pxp_constraint_mpo`) and lift it with
    `pauli_superoperator_mpo`. The same pattern also fits DAOE-style projectors
    ([DAOE](daoe.md)).

## Higher spin

Every DMT entry point is generic in the local Hilbert space dimension `d`: build sites with
`operator_siteinds(nsites; d)` instead of `pauli_siteinds(nsites)`, and `DMTGateEvolution`,
`dmt_step!`, `dmt_evolve!`, and `DMTOptions` all take `d` from the sites without further
configuration. The onsite basis is the generalized Gell-Mann basis (see
[Operator Space](operator-space.md)); at `d = 2` it is exactly the Pauli basis used throughout
the rest of this page, so a spin-1/2 script written against `pauli_*` names is already the
`d = 2` special case of everything below.

Building the dense local Hamiltonian that `operator_gate_from_hamiltonian` needs is outside
this package's scope for spin `> 1/2` (`MPSToolkit` supplies the operator-space machinery, not
spin-`S` model builders — see [Operator Space](operator-space.md)); any source of dense
`d x d` matrices works, sparse included, since the gate builders densify internally. The
snippet below uses [EDKit.jl](https://github.com/jayren3996/EDKit.jl) — **not a dependency of
MPSToolkit** — as one convenient source, because its multi-site string convention (`kron`
applied left-to-right, first character = leftmost site) matches MPSToolkit's own:

```julia
using MPSToolkit, ITensors, ITensorMPS
using LinearAlgebra: I
using EDKit: spin                       # illustrative only: EDKit.jl is not an MPSToolkit
                                         # dependency; any source of dense d x d matrices works

d, nsites, dt, maxdim = 3, 40, 0.1, 64
sites = operator_siteinds(nsites; d = d)

h = spin((1, "xx"), (1, "yy"), (1, "zz"); D = d)      # sparse spin-1 Heisenberg bond
gate = operator_gate_from_hamiltonian(h, dt; d = d)   # sparse input densified internally

sz = Matrix(spin("z"; D = d))
identity_matrix = Matrix{ComplexF64}(I, d, d)
rho = add(operator_product_state(sites, fill(identity_matrix, nsites)),
          operator_local_sum_state(sites, sz, [j <= nsites ÷ 2 ? -0.25 : 0.25 for j in 1:nsites]);
          maxdim = 8, cutoff = 0.0)

schedule = collect(1:(nsites - 1))
evo = DMTGateEvolution(gate, dt; schedule, reverse_schedule = reverse(schedule),
                       maxdim = maxdim, cutoff = 1e-12)
profile(state) = real.(operator_expectation_profile(state, [(x, sz) for x in 1:nsites]))

println("t=0.0  center = ", round.(profile(rho)[19:22]; digits=4))
for k in 1:3
    dmt_evolve!(rho, evo)
    println("t=", round(2dt * k; digits=2), "  center = ", round.(profile(rho)[19:22]; digits=4),
            "  maxlinkdim=", maxlinkdim(rho))
end
# t=0.0  center = [-0.1667, -0.1667, 0.1667, 0.1667]
# t=0.2  center = [-0.1666, -0.1581, 0.1581, 0.1666]  maxlinkdim=59
# t=0.4  center = [-0.1658, -0.1349, 0.1349, 0.1658]  maxlinkdim=63
# t=0.6  center = [-0.1628, -0.1039, 0.1039, 0.1628]  maxlinkdim=63
```

`rho` is the `S^z` domain wall ``I + \sum_j c_j S^z_j`` (``c_j = \mp 0.25``) on top of the
literal identity, exactly analogous to the spin-1/2 melts above. Read back through
`operator_expectation_profile`, the physical ``S^z`` profile starts at the clean value
``\pm 1/6`` (the normalized-basis overlap of ``S^z`` with itself at `d = 3`) and smooths from a
step toward the chain center as the DMT-truncated Heisenberg evolution proceeds — the same
melting behavior as the `d = 2` examples, now at spin 1.

At `d = 3` the floor for `preserve_diameter = 3` is `maxdim >= 2*3^2 + 1 = 19` (see the budget
table above); `maxdim = 64` here leaves `chi' = 64 - 18 = 46` complement directions, a
demonstration-scale budget well below the `chi = 128-256` that arXiv:2205.02853 used for
converged spin-1 SU(3) transport.

For the production version of this melt — spin-1 Heisenberg against the integrable ULS/SU(3)
point, a `maxdim` convergence ladder, the front-containment guard, and the extracted dynamical
exponents — see
[examples/dmt/spin1_melt.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/spin1_melt.jl).

!!! warning "Near the DMT budget floor, DMT is measurably worse than the plain SVD it replaces"
    [examples/dmt/spin1_semiexact_validation.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/spin1_semiexact_validation.jl)
    runs the same melt at `N = 12` against an **uncapped reference** and scores DMT against plain
    SVD `LocalGateEvolution` at the *same* `maxdim`, same gates, same schedule. The reference is
    itself cutoff-limited, so its own error is measured per case and per time and subtracted: every
    figure below is a guaranteed bound from the interval
    `[(svd − floor)/(dmt + floor), (svd + floor)/(dmt − floor)]`, truncated toward 1, and cells
    whose interval contains 1 are reported as unresolved rather than as results.

    **The robust result is negative.** At the tightest budget the kernel admits (`chi' = 2`, i.e.
    `maxdim = 20` at `d = 3`), DMT's profile error at each case's last **floor-covered** time is
    **≥ 38.7x worse** than plain SVD's for spin-1 Heisenberg (`t = 1.2`), **≥ 43.7x worse** for
    ULS/SU(3) (`t = 1.2`), and **≥ 13.0x worse** in a `chi'`-matched `d = 2` control (`t = 2.6`) —
    same sign in every wall-amplitude and cutoff cell. (At Heisenberg's last front-*contained*
    time, `t = 1.6`, the raw loss is 11.2x, but the reference's error is unmeasured there; the two
    clocks are different and the labels matter.) The `d = 2` control reproduces the loss at
    `chi' = 2` and `chi' = 10`, so this is a property of the method, not of higher spin.

    **No kernel-only win is claimed at large `chi'`.** At `chi' = 46` both arms sit at or below the
    reference's resolution wherever that resolution is known.

    **No budget-only rule fits the data.** At the same `d = 3`, same `chi' = 22`, same `t = 1.2`,
    spin-1 Heisenberg is a resolved **loss** (≥ 1.7x) while ULS is a resolved **win** (≥ 1.7x) —
    the model decides, not the budget. The protected-block *fraction* does not predict it either:
    ULS at `chi' = 10` (64% overhead) wins by ≥ 5.0x while `d = 2` Heisenberg at the same `chi'`
    (44% overhead) loses by ≥ 28.9x. So run the equal-`maxdim` comparison for your own model
    rather than trusting any constant; `chi' = 2` is a resolved loss everywhere tested, and the
    crossover moves by at least 2.2x in `chi'` between two models at the same `d`.

    **What is resolved on the positive side is immunity to a norm-relative cutoff.** With an
    infinitesimal wall and the ordinary `cutoff = 1e-10`, the plain-SVD arm hits an error floor
    that more bond dimension does not remove (`1.02e-4` at `chi' = 22` versus `1.03e-4` at
    `chi' = 46`), because `cutoff` bounds discarded weight relative to the *full* Hilbert-Schmidt
    norm while the signal is a few percent of it. DMT is structurally immune: **≥ 2.1x** better at
    `chi' = 22` and **≥ 29.0x** at `chi' = 46` (`d = 3`), **≥ 19.1x** and **≥ 34.7x** at the
    matched `d = 2` rungs. Not universal, though — ULS at `chi' = 46` is a resolved win of ≥ 3.5x
    at `t = 1.0` and a resolved *loss* of ≥ 1.5x at `t = 1.2` in that same cell. DMT is separately
    invariant to the cutoff itself (error moves ≤ 0.02% at `chi' >= 10`, against 34–268% for plain
    SVD) and holds `tr(rho)` to `~1e-13` where plain SVD drifts to `1e-4`. Its *wall-amplitude*
    invariance is **not** unconditional: 272% spread across `eps` at `chi' = 2`.

!!! warning "`S^z` is not a single basis element at `d >= 3`"
    Unlike at `d = 2`, where `pauli_basis_state(sites, ["Z", ...])` selects the single-`Z`
    direction directly, the physical `S^z = diag(1, 0, -1)` at `d = 3` is a **combination** of
    two diagonal Gell-Mann generators, not a scalar multiple of either alone — so no integer
    basis label represents it, and `operator_basis_state` cannot build it (it only ever selects
    one basis direction). Build it from the dense matrix with `operator_product_state` /
    `operator_local_sum_state` instead, as above. See [Operator Space](operator-space.md) for
    the full basis ordering and this caveat in more detail.

`spin(...; D = d)` composes multi-site operator strings with `kron` applied left-to-right,
first character = leftmost site — the same convention `operator_gate` and every
operator-space state builder in MPSToolkit use — so dense Hamiltonians from EDKit.jl (or any
other source following the same convention) plug in with no site relabeling.

## Relation to TEBD

DMT and TEBD share one scheduling abstraction and differ only in the truncation kernel:

- [`LocalGateEvolution`](@ref) and `DMTGateEvolution` both represent scheduled local updates: a
  gate (or gate provider), a schedule of bond positions, and a per-step truncation budget.
- `tebd_evolve!` applies each gate and truncates with an ordinary SVD that keeps the largest
  singular values — the correct choice for compressing a state by entanglement, and the
  general-purpose operator-space truncation.
- `dmt_step!` applies the gate the same way (it even reuses the TEBD gate-application routine
  internally for the raw step) but then truncates with the identity-preserving DMT rule.

Keeping the two on the same footing means switching from generic operator-space TEBD to a
transport simulation is a change of *truncation backend*, not a different evolution family. For
the TEBD side of this story — schedules, gate providers, Strang splitting, and the SVD
truncation budget — see [TEBD](tebd.md). For the Pauli-basis vectorization,
state builders, and gate builders shared by both, see [Operator Space](operator-space.md).

For *when* to prefer DMT, DAOE, or plain TEBD — and in particular why a dynamical exponent
``z`` and a diffusion constant ``D`` place different demands on the method — see
[Method Comparison](transport-methods.md).

## Tips and pitfalls

- **Use DMT only for the density-operator melt.** The whole point of the method is the
  identity-preserving bias, so it applies only to a near-infinite-temperature *density operator*
  whose transport you track (the domain-wall melt). If you are evolving a *traceless* operator —
  a two-point / autocorrelation function of a conserved density, generic operator growth,
  out-of-time-order correlators, or an arbitrary observable's Heisenberg evolution — that bias is
  wrong (there is no trace component to anchor it), and ordinary operator-space TEBD truncation
  ([`LocalGateEvolution`](@ref) with `tebd_evolve!`) is the right tool.
- **Tune `maxdim` and `gate_maxdim` together; `preserve_diameter` sets a floor, not a dial to
  balance against them.** A large `gate_maxdim` lets a gate inflate the bond a lot before DMT
  truncates it back to `maxdim`; the default `gate_maxdim = 0` (apply the gate exactly) is safe
  everywhere this package's own tests reach (`d <= 4`) and is a correctness fix at `d >= 5`, so
  raise it above `0` only to cap peak memory, never for accuracy. `preserve_diameter` is fixed
  by which diameter of observable needs preserving, not traded off against the others; `maxdim`
  must simply clear the floor `2 d^(2n) + 1` it implies (see the budget table above).
- **Pick `maxdim` by convergence, not by guesswork.** As with any MPS method, increase
  `maxdim` until the transport observable of interest (e.g. a diffusion constant or a density
  profile) stops moving. DMT typically reaches converged hydrodynamics at far smaller bond
  dimensions than plain SVD truncation would for the same operator, which is its main practical
  advantage — but this is not automatic, and near the budget floor it inverts. At a budget close
  to the protected block, DMT spends nearly everything on structure and is measurably *worse* than
  the plain SVD it replaces. No budget-only rule predicts where that changes: the `d = 3`
  measurement in
  [examples/dmt/spin1_semiexact_validation.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/spin1_semiexact_validation.jl)
  finds DMT behind at `chi' = maxdim - (2 d^(2n) - 1) = 3` in every model tested, and finds spin-1
  Heisenberg behind while ULS is ahead **at the same `d`, the same `chi' = 22` and the same time**.
  So treat the equal-`maxdim` DMT-vs-SVD comparison as a convergence check to run for your own
  model, budget `chi'` well into the tens, and do not read a crossover constant off any single
  model — it moved by at least 2.2x in `chi'` between two models here.
- **`cutoff` is a repair-SVD setting, not the primary control.** It governs only the final
  re-factorization of the already-truncated bond. The transport-relevant decision is made by
  the DMT rule itself, so `maxdim` and `preserve_diameter` are the dials that matter.
- **Higher-spin budgets grow fast.** The protected block has dimension `2 d^(2n)`
  (`n = (preserve_diameter - 1) / 2`), so at fixed `preserve_diameter` a spin-1 (`d = 3`) run
  needs a substantially larger `maxdim` than the equivalent spin-1/2 (`d = 2`) run for the same
  complement resolution — compare the `d = 2` and `d = 3` columns of the budget table above.
  arXiv:2205.02853 ran `d = 3` at `chi = 128-256`.

## API

```@docs
DMTOptions
DMTGateEvolution
dmt_step!
dmt_evolve!
constrained_dmt_evolve!
```

PXP model and constraint helpers used by the constrained-transport workflow:

```@docs
pxp_term_hamiltonian
pxp_term_support
pxp_constraint_mpo
pauli_pxp_constraint_state
pauli_pxp_constraint_projector
```

## Examples

- [examples/operator_space/dmt_scheduler.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/dmt_scheduler.ipynb)
- [examples/dmt/domain_wall_melting.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/domain_wall_melting.jl)
- [examples/dmt/pxp_energy_melting.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/pxp_energy_melting.jl) — the energy domain-wall melt (the DMT protocol)
- [examples/dmt/spin1_melt.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/spin1_melt.jl) — higher spin (`d = 3`): spin-1 magnetization melt, diffusive Heisenberg vs. KPZ at the integrable ULS/SU(3) point, with a `maxdim` convergence ladder
- [examples/dmt/spin1_semiexact_validation.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/spin1_semiexact_validation.jl) — higher spin (`d = 3`): DMT vs. plain SVD truncation at equal `maxdim` against an uncapped reference, with a `d = 2` control; measures what the truncation kernel buys instead of fitting an exponent
- [examples/operator_space/pxp_energy_correlator.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/pxp_energy_correlator.jl) — *off-label*: evolves a traceless operator with DMT (see the warning above); prefer TEBD for correlators
- [examples/open_systems/pauli_lindblad_tebd.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/open_systems/pauli_lindblad_tebd.ipynb)

## References

- C. David White, Michael Zaletel, Roger S. K. Mong, and Gil Refael, [Quantum dynamics of thermalizing systems](https://arxiv.org/abs/1707.01506) — the DMT algorithm; Sec. III for the exact-preservation guarantee.
- Bingtian Ye, Francisco Machado, Christopher David White, Roger S. K. Mong, and Norman Y. Yao, [Emergent hydrodynamics in nonequilibrium quantum systems](https://arxiv.org/abs/1902.01859) — the `maxdim = chi_preserve + chi_extra` total-bond-dimension convention used here.
- Bingtian Ye, Francisco Machado, Jack Kemp, Ross B. Hutson, and Norman Y. Yao, [Universal KPZ dynamics in integrable quantum systems](https://arxiv.org/abs/2205.02853) — precedent for DMT at `d = 3` (SU(3) spin-1) and `d = 4`, at `chi` up to 256-512.
- Stuart Yi-Thomas, Brayden Ware, Jay D. Sau, and Christopher David White, [Comparing numerical methods for hydrodynamics in a one-dimensional lattice spin model](https://arxiv.org/abs/2310.06886)
- En-Jui Kuo, Brayden Ware, Peter Lunts, Mohammad Hafezi, and Christopher David White, [Energy diffusion in weakly interacting chains with fermionic dissipation-assisted operator evolution](https://arxiv.org/abs/2311.17148)
- Marko Ljubotina, Jean-Yves Desaules, Maksym Serbyn, and Zlatko Papić, [Superdiffusive energy transport in kinetically constrained models](https://doi.org/10.1103/PhysRevX.13.011033), Phys. Rev. X 13, 011033 (2023) — PXP energy transport and the KPZ exponent.
- Paul Fendley, K. Sengupta, and Subir Sachdev, [Competing density-wave orders in a one-dimensional hard-boson model](https://arxiv.org/abs/cond-mat/0309438), Phys. Rev. B 69, 075106 (2004) — origin of the PXP / hard-square constraint.
- C. J. Turner, A. A. Michailidis, D. A. Abanin, M. Serbyn, and Z. Papić, [Weak ergodicity breaking from quantum many-body scars](https://arxiv.org/abs/1711.03528), Nature Physics 14, 745 (2018) — quantum many-body scars in the constrained PXP sector.
