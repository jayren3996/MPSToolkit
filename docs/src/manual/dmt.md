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
basis. MPSToolkit uses the normalized Pauli basis ``(I, X, Y, Z)`` for this; see the
[Operator Space](operator-space.md) page for the vectorization conventions, `pauli_siteinds`,
`pauli_basis_state`, and the gate builders.

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

### The identity-preserving truncation rule

DMT changes which information is protected. At each bond it forms the *reduced* description of
the two-site (left block, right block) operator as seen through the trace — concretely, it
contracts the left and right halves of the operator MPS against the local identity
(trace) environment, so that the matrix being truncated is expressed in a basis whose first
direction is the identity/trace connector and whose leading directions span the local reduced
operators on each side. Truncation then proceeds on the *remaining* block only:

```math
M \;=\; \underbrace{C}_{\text{protected connector}} \;+\; \underbrace{M'}_{\text{truncated by SVD}} ,
```

where ``C`` is the rank-one identity/trace connector that is subtracted off and reinstated
*exactly*, and the SVD truncation acts only on ``M'`` after a buffer of leading "connector"
directions has been set aside. Because the identity direction and a buffer of nearby low-weight
directions are carried through every bond untouched, the locally conserved densities and their
slow hydrodynamic tails survive truncation even though they correspond to small singular
values. What gets compressed is the high-weight, long-range connected correlations, the part
of the operator that does not feed back into transport.

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
  `connector_buffer`) so they can be passed to a single step without constructing a full
  evolution.
- `dmt_step!` is the lowest-level entry point: apply one local gate at a bond and then perform
  the associated DMT truncation. It exists both in a keyword form and in an overload that takes
  a `DMTOptions`.
- `dmt_evolve!` is the scheduled driver. It walks the forward schedule applying `dmt_step!`
  with `direction=:R`, then walks the reverse schedule with `direction=:L`, repeats that for
  `evo.nstep` sweeps, and finally normalizes the operator MPS.

`dmt_step!` first applies the raw gate with an inflated bond budget (`gate_maxdim`) and *no*
cutoff, then truncates the inflated bonds back down to `maxdim` using the DMT rule. In this
two-stage structure the gate is allowed to expand the bond and the transport-aware kernel then
compresses it, which is why `gate_maxdim`, `maxdim`, and `connector_buffer` should be chosen
together rather than independently.

The truncation budget has a few knobs worth understanding:

- `maxdim`: target bond dimension *after* DMT truncation.
- `cutoff`: cutoff used only in the final "repair" SVD that re-factorizes the truncated bond
  back into MPS form.
- `gate_maxdim`: temporary bond-dimension ceiling allowed while the raw gate is applied, before
  DMT compresses the bond. It defaults to `max(maxdim * 16, 64)`.
- `connector_buffer`: number of leading connector directions (starting with the identity/trace
  direction) that are carried through every bond untouched. It defaults to `8` and must satisfy
  `connector_buffer <= maxdim`.

DMT assumes the operator MPS lives on dimension-4 Pauli sites ordered `(I, X, Y, Z)`; the step
validates this and rejects non-Pauli sites. Periodic-boundary local gates are not supported.

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
  sites, so its induced Pauli-space gate spans two dimension-4 operator sites. `dmt_evolve!`
  infers this span automatically from the gate size and the site dimensions.
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
        maxdim, cutoff=1e-10, gate_maxdim=4 * maxdim, connector_buffer=4)
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
# 0.5  ->  z = 1.09    (ballistic,            z = 1)
# 1.0  ->  z = 1.45    (superdiffusive / KPZ, z = 3/2)
# 2.0  ->  z = 1.90    (diffusive,            z = 2)
```

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
transport is *superdiffusive*, with dynamical exponent ``z = 3/2`` in the KPZ universality
class (Ljubotina, Desaules, Serbyn, Papić, PRX **13**, 011033 (2023)).

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
   with ``P_G``, the preparation stays in the constrained sector up to truncation error. For
   a correlator run the initial state is instead a sector-projected local energy density,
   built by vectorizing the term's MPO and applying the constraint projector once.
3. **Constraint checkpoints during evolution.** Exact gates commute with ``P_G``, but DMT (or
   any) truncation leaks weight out of the sector. The operator-space projector MPO
   ``\rho \mapsto P_G \rho P_G`` — `pauli_pxp_constraint_projector(sites)`, bond dimension 4,
   built by the generic `pauli_superoperator_mpo` — removes the leaked weight, and
   `constrained_dmt_evolve!` interleaves it with the DMT sweeps.

Two ready-made protocols use these ingredients; both ship as full scripts under
[`examples/operator_space/`](https://github.com/jayren3996/MPSToolkit/tree/main/examples/operator_space).

### Protocol 1 — infinite-temperature energy correlator (the efficient exponent route)

Heisenberg-evolve a single **sector-projected local energy density**,
``O(0) \propto P_G h_c P_G``, and watch it spread: the profile
``C(x, t) = \mathrm{tr}(P_G h_x O(t))`` is the constrained infinite-temperature
energy-energy correlator, and its second moment grows as
``M_2(t) \sim t^{2/z}`` — asymptotic log-log slope ``4/3`` for ``z = 3/2``. This is the
protocol of PRX **13**, 011033, condensed from
[`examples/operator_space/pxp_energy_correlator.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/pxp_energy_correlator.jl):

```julia
using MPSToolkit, ITensors, ITensorMPS

nsites, c, chi, dt = 48, 24, 32, 0.1
sites = pauli_siteinds(nsites)
terms = [(first(pxp_term_support(nsites, j)), pxp_term_hamiltonian(nsites, j)) for j in 1:nsites]
projector = pauli_pxp_constraint_projector(sites)

phys = siteinds("S=1/2", nsites)                  # O(0) = P_G h_c P_G, vectorized
os = OpSum(); os += "ProjUp", c - 1, "X", c, "ProjUp", c + 1   # "Up" = our |0>
O = pauli_state_from_mpo(MPO(os, phys), sites)
O = apply(projector, O; maxdim=4chi, cutoff=1e-12); normalize!(O)

gates = [pauli_gate_from_hamiltonian(h, dt) for (_, h) in terms]
schedule = [start for (start, _) in terms]
evo = DMTGateEvolution(gates, dt; schedule, reverse_schedule=reverse(schedule),
    nstep=1, maxdim=chi, cutoff=1e-12, gate_maxdim=4chi, connector_buffer=8)

M2 = Float64[]
for k in 1:60                                     # one call advances t by 2*dt
    constrained_dmt_evolve!(O, evo, projector; project_every=1, normalize=false)
    p = real.(pauli_expectation_profile(O, terms; normalize=false))  # traceless: unnormalized
    push!(M2, sum((j - c)^2 * p[j] for j in 1:nsites) / sum(p))
end
```

Two conventions in this snippet are load-bearing for a **traceless** operator:
`normalize=false` in the evolution (truncation sheds Hilbert–Schmidt norm faster than the
DMT-protected components, so the default per-sweep renormalization silently inflates
absolute traces — and with it the conserved total ``\sum_x C(x,t) = \mathrm{tr}(P_G H\,
O(t))``, which with `normalize=false` becomes the run's truncation error bar), and
`normalize=false` in the profile measurement (the trace it would divide by vanishes).
Production-validated reference points for the crossover (``N = 64``), committed as a
reproducible reference in `examples/operator_space/data/pxp_energy_correlator_chi48_N64.csv`
(re-fit by `examples/operator_space/pxp_energy_correlator_reference.jl`, and pinned in
`test/test_pxp.jl`): the local ``M_2`` slope plateaus at **1.46–1.50 over** ``t = 8\!-\!14``
at ``\chi = 48`` (LSQ ``2/z = 1.48``), robustly above the diffusive ``1.0`` and **descending
toward** ``4/3`` **from above** — it has *not* reached the asymptote. The `normalize=false`
conserved-total drift is the run's truncation error bar: it stays within **1–6% only through**
``t \approx 10`` and grows to **~9–12% by** ``t = 12\!-\!14``, so the top of the window is the
least trustworthy part of the plateau, and reaching the asymptotic ``4/3`` (expected only at
``t \gtrsim 20\!-\!30``) requires ``\chi > 48`` — not merely longer times. By contrast
``\chi = 32`` (the shipped demo default) reads **1.30–1.38**: truncation suppresses ``M_2``
growth and pulls the slope *down*, so an undershooting slope near the target there is a
truncation artifact, not convergence. Raise `maxdim` until the slope stops moving up, then
extend in time.

### Protocol 2 — energy domain-wall melt

The quench protocol: prepare ``\rho(0) \propto e^{-\beta H_L} \otimes e^{+\beta H_R}`` and
track the transferred energy. From
[`examples/operator_space/pxp_energy_transport.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/pxp_energy_transport.jl):

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
    nstep=1, maxdim=32, cutoff=1e-12, gate_maxdim=128, connector_buffer=8)
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
a useful truncation gauge**. Both protocols are mirror-symmetric setups, so the profile
asymmetry should vanish in exact arithmetic; in practice it grows with truncation pressure
and *shrinks* with increasing ``\chi`` (it is truncation noise, not a sweep-direction bias).
Monitoring it alongside the conserved total gives a cheap convergence diagnostic.

In both protocols the measurement is `pauli_expectation_profile(rho, terms)`, which
evaluates every ``\mathrm{tr}(\rho h_j)`` (over ``\mathrm{tr}(\rho)`` when normalized) in
one ``O(N)`` pass with cumulative identity environments, mixed window sizes included. The
full scripts add the diagnostics discussed above — the local-slope crossover column, the
conserved-total drift, and the constraint-leakage probe (residual leakage after projection
vs. leakage accrued by one unprojected sweep):
[`pxp_energy_correlator.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/pxp_energy_correlator.jl)
and
[`pxp_energy_transport.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/pxp_energy_transport.jl).

!!! tip "The checkpoint pattern is generic"
    Nothing in `constrained_dmt_evolve!` is PXP-specific: it interleaves DMT sweeps with any
    operator-space MPO. Other kinetically constrained models need only their own constraint
    MPO — build it in physical space (cf. `pxp_constraint_mpo`) and lift it with
    `pauli_superoperator_mpo`. The same pattern also fits DAOE-style projectors
    ([DAOE](daoe.md)).

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

## Tips and pitfalls

- **Use DMT only for transport.** The whole point of the method is the identity-preserving
  bias. If you are computing something that is *not* a conserved-density transport quantity —
  generic operator growth, out-of-time-order correlators, an arbitrary observable's Heisenberg
  evolution — that bias is wrong and ordinary operator-space TEBD truncation
  ([`LocalGateEvolution`](@ref) with `tebd_evolve!`) is the right tool.
- **Tune `maxdim`, `gate_maxdim`, and `connector_buffer` together.** A large `gate_maxdim`
  lets a gate inflate the bond a lot before DMT truncates it back to `maxdim`; with a small
  `connector_buffer` that can reduce accuracy, because fewer low-weight directions are
  protected through the compression. Treat the three as a single budget. The constraint
  `connector_buffer <= maxdim` is enforced by the constructors.
- **Pick `maxdim` by convergence, not by guesswork.** As with any MPS method, increase
  `maxdim` until the transport observable of interest (e.g. a diffusion constant or a density
  profile) stops moving. DMT typically reaches converged hydrodynamics at far smaller bond
  dimensions than plain SVD truncation would for the same operator, which is its main practical
  advantage.
- **`cutoff` is a repair-SVD setting, not the primary control.** It governs only the final
  re-factorization of the already-truncated bond. The transport-relevant decision is made by
  the connector-preserving step, so `maxdim` and `connector_buffer` are the dials that matter.
- **Pauli sites only.** DMT requires dimension-4 ``(I, X, Y, Z)`` operator sites and rejects
  anything else; it is not a state-space or arbitrary-local-dimension method. Periodic-boundary
  local gates are not implemented.

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
- [examples/operator_space/pxp_energy_correlator.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/pxp_energy_correlator.jl)
- [examples/operator_space/pxp_energy_transport.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/operator_space/pxp_energy_transport.jl)
- [examples/open_systems/pauli_lindblad_tebd.ipynb](https://github.com/jayren3996/MPSToolkit/blob/main/examples/open_systems/pauli_lindblad_tebd.ipynb)

## References

- C. David White, Michael Zaletel, Roger S. K. Mong, and Gil Refael, [Quantum dynamics of thermalizing systems](https://arxiv.org/abs/1707.01506)
- Stuart Yi-Thomas, Brayden Ware, Jay D. Sau, and Christopher David White, [Comparing numerical methods for hydrodynamics in a one-dimensional lattice spin model](https://arxiv.org/abs/2310.06886)
- En-Jui Kuo, Brayden Ware, Peter Lunts, Mohammad Hafezi, and Christopher David White, [Energy diffusion in weakly interacting chains with fermionic dissipation-assisted operator evolution](https://arxiv.org/abs/2311.17148)
- Marko Ljubotina, Jean-Yves Desaules, Maksym Serbyn, and Zlatko Papić, [Superdiffusive energy transport in kinetically constrained models](https://doi.org/10.1103/PhysRevX.13.011033), Phys. Rev. X 13, 011033 (2023) — PXP energy transport and the KPZ exponent.
- Paul Fendley, K. Sengupta, and Subir Sachdev, [Competing density-wave orders in a one-dimensional hard-boson model](https://arxiv.org/abs/cond-mat/0309438), Phys. Rev. B 69, 075106 (2004) — origin of the PXP / hard-square constraint.
- C. J. Turner, A. A. Michailidis, D. A. Abanin, M. Serbyn, and Z. Papić, [Weak ergodicity breaking from quantum many-body scars](https://arxiv.org/abs/1711.03528), Nature Physics 14, 745 (2018) — quantum many-body scars in the constrained PXP sector.
