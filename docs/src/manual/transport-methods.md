# Choosing a transport method: TEBD, DMT, and DAOE

*When to reach for plain operator-space TEBD, for the identity-preserving DMT truncation, or for DAOE's size-dependent dissipation — and why the answer depends on whether you want a dynamical exponent or a diffusion constant.*

`MPSToolkit.jl` offers three ways to extract high-temperature transport from operator-space evolution. They are not really competitors; they are three points on one trade-off surface, and a single study often uses more than one. This page is both the decision guide and the theory behind it — the mechanics of each method live on the [TEBD](tebd.md), [DMT](dmt.md), and [DAOE](daoe.md) pages, and the shared vectorized-Pauli machinery on [Operator Space](operator-space.md).

Two orthogonal questions fix the choice:

1. **What are you evolving — a trace-ful density operator, or a traceless one?** This decides which truncation is *valid*.
2. **What do you want — the dynamical exponent ``z`` or the diffusion constant ``D``?** This decides which *bias* you can tolerate and how *long* you have to run.

## The shared bottleneck: operator-space entanglement

In the Heisenberg picture an operator evolves as ``O(t) = e^{iHt}\, O\, e^{-iHt}``, and in a chaotic chain it spreads into an exponentially growing cloud of long Pauli strings (the [DAOE](daoe.md) page draws this picture explicitly). Represented as a matrix product operator — equivalently a matrix product *state* in the doubled Pauli space — this proliferation appears as **operator-space entanglement** (the Hilbert–Schmidt entanglement of the vectorized operator), and in a thermalizing system it grows *linearly* in time. Faithful SVD truncation must then grow the bond dimension as ``\chi \sim e^{c t}``, so a naive evolution stalls after an ``O(1)`` time. This single fact is why DMT and DAOE exist.

The escape is physical, not numerical. Hydrodynamics — diffusion constants, conductivities, the slow tails of autocorrelators — is carried entirely by the **slow, local, conserved** part of the operator. The fast-spreading, high-weight cloud holds almost all of the entanglement but, under chaotic dynamics, behaves like a structureless bath that no longer feeds back coherently into few-body observables. The three methods differ in *how*, and by *what mechanism*, they exploit this:

- **TEBD** exploits it not at all — it truncates by entanglement (largest singular values), so it pays the full ``\chi \sim e^{c t}`` cost. Unbiased reference, exact up to ``\chi``, but bond-starved.
- **DMT** keeps the **local reduced-operator data** and discards long-range *connected* correlations — the part hydrodynamics ignores.
- **DAOE** keeps short operators and **physically damps** the long strings — removing the cloud that drives the entanglement growth.

DMT and DAOE thus reach far longer times than TEBD at the same ``\chi``, each by throwing away a different slice of the hydrodynamically-inert information.

## One backbone, three truncations

All three methods are the same algorithm — schedule local gates, apply them to a vectorized operator in the ``(I, X, Y, Z)`` Pauli basis, and compress. They differ *only* in the compression step.

- **TEBD** — ordinary SVD truncation (keep the largest singular values), via [`LocalGateEvolution`](@ref). The **unbiased baseline**: exact up to ``\chi``, valid for *any* operator. In the worked example below it leaks *several percent* of the conserved energy by ``t \sim 10`` at a moderate ``\chi`` — the truncation is already corrupting the hydrodynamics, exactly as the ``\chi \sim e^{ct}`` wall predicts.

- **DMT** — the identity-preserving truncation, via [`DMTGateEvolution`](@ref) / [`dmt_evolve!`](@ref). Instead of keeping the largest singular values it protects the **reduced operators on every window up to a fixed diameter** — the trace/identity component and its short-range Pauli dressing — and discards only the *connected* (cumulant) part of longer-range correlations. By linear response, that protected data is exactly what fixes the local densities and currents, i.e. the hydrodynamic content, so DMT is "hydrodynamically complete" by construction and converges in ``\chi`` far faster than SVD for transport. The flip side: it anchors on the trace, so it is valid **only on trace-ful, near-infinite-temperature density operators**. On a traceless operator the bias has nothing to anchor and is simply wrong (see the warning on the [DMT](dmt.md) page).

- **DAOE** — dissipation-assisted operator evolution, TEBD interleaved with [`pauli_daoe_projector`](@ref). A diagonal dissipator exponentially damps Pauli strings longer than a cutoff ``\ell^\ast`` (strength ``\gamma``). The induced bias has a clean interpretation through **operator backflow** (von Keyserlingk, Pollmann & Rakovszky, [arXiv:2111.09904](https://arxiv.org/abs/2111.09904)): the discarded long operators would otherwise feed current back into the conserved mode, so removing them *lowers* the transport coefficient — but that correction is **exponentially small in ``\ell^\ast``**, which is precisely why the ``\gamma \to 0`` / ``\ell^\ast \to \infty`` extrapolation converges quickly and controllably. The danger is the same mechanism near integrability: there the (quasi-)local conserved charges that protect ballistic or anomalous transport live on long strings, DAOE damps them, and it then **manufactures diffusion that is not there** (Yoo, White & Swingle, [arXiv:2210.06494](https://arxiv.org/abs/2210.06494)).

A compact way to hold them in mind: **TEBD is unbiased but bond-starved; DMT trades generality for reach on trace-ful densities; DAOE trades exactness for reach on operators, paying it back by extrapolation.**

## Trace structure and linear response: which method is valid

The melt and the autocorrelation are two faces of one object. A weakly biased domain-wall density ``\rho \approx 2^{-N} I + \varepsilon\, D`` (with ``D`` the traceless domain-wall operator) is, to ``O(\varepsilon)``, exactly the linear-response source whose evolution produces the infinite-temperature density autocorrelation. The **trace-ful** density (the melt) and the **traceless** correlator therefore differ only by the identity background — but that background is precisely what DMT protects. Hence:

- a near-``\infty``-``T`` density operator — a domain-wall melt, a thermal state, a sinusoidal energy mode — is the *right* input for DMT (and the trace it protects is physical);
- a traceless operator — a domain-wall *operator* evolved on its own, a current or energy-density autocorrelation, an OTOC — carries no identity for DMT to anchor, so it must use TEBD, optionally DAOE-assisted.

This is the guardrail enforced by [`dmt_evolve!`](@ref) / [`constrained_dmt_evolve!`](@ref) throughout the operator-space code: melt the density with DMT, spread the traceless density with TEBD(+DAOE).

## The exponent and the constant are different problems

This is the distinction that most often decides the method. A dynamical exponent ``z`` and a diffusion constant ``D`` are not two readouts of one quantity — they have different definitions and different numerical demands.

### The exponent ``z``

``z`` sets the *scaling*: a conserved density spreads as ``x \sim t^{1/z}``, read off as a log–log slope over a range of times. It classifies the transport into universality classes:

| ``z`` | transport | hallmark | finite ``D``? |
|---|---|---|---|
| ``1`` | ballistic | finite Drude weight | no — ``D`` diverges |
| ``3/2`` | superdiffusive (KPZ) | running ``D(t)\sim t^{1/3}`` | no |
| ``2`` | diffusive | finite diffusion constant | **yes** |
| ``>2`` (e.g. ``4``) | subdiffusive | dipole / multipole conservation | no — ``D = 0`` |

Two requirements follow:

- ``z`` only takes its asymptotic value once the hydrodynamic scaling regime is reached, which for near-integrable or kinetically constrained models can be a **long** time. (Bare PXP energy transport, for example, looks near-ballistic until ``t \sim 100`` before bending toward its asymptotic value.)
- ``z`` must be measured **without bias**. As above, DAOE's dissipation pushes transport toward ``z = 2`` and can manufacture spurious diffusion. **Never use DAOE to decide an unknown or anomalous exponent.**

So determining ``z`` calls for an **unbiased method run long**: DMT on the trace-ful melt (efficient and recommended), or plain TEBD (unbiased but bond-limited). This is exactly why anomalous-``z`` problems are hard — the exponent *is* the question, the cheap tool would lie, and the times required are large.

### The constant ``D``

``D`` is a transport *coefficient* that **exists only when ``z = 2``** (the table above: at any other ``z`` it diverges or vanishes). Its two standard definitions — the mean-square displacement (MSD) and the Green–Kubo current integral — are linked by an exact identity worth seeing, because it explains why either route works and what the numerics must achieve.

Write ``C(x,t) = \langle q_x(t)\, q_0(0)\rangle_c`` for the connected infinite-temperature autocorrelation of the conserved density ``q_x``, with lattice continuity ``\partial_t q_x = j_{x-1} - j_x`` and total current ``J = \sum_x j_x``. Conservation makes the *total* weight time-independent,

```math
\sum_x C(x,t) \;=\; \big\langle\, \textstyle\sum_x q_x(t)\; q_0(0)\big\rangle_c \;=\; \langle Q\, q_0\rangle_c \;=\; \chi,
```

the static susceptibility (the ``\infty``-``T`` variance of the density). This conserved total is exactly the ``\sum_b C_b(t) = \operatorname{tr}[h_c(t) H]`` quantity whose drift the worked example prints as a truncation error bar. The normalized second moment ``M_2(t) = \chi^{-1}\sum_x x^2\, C(x,t)`` then obeys the **running-diffusion-constant identity**

```math
D(t) \;\equiv\; \tfrac{1}{2}\,\frac{d M_2(t)}{dt}
\;=\; \frac{1}{\chi}\int_0^t \langle J(0)\, J(t')\rangle\, dt'
\qquad\xrightarrow{\;t\to\infty\;}\; D
```

(up to lattice normalization). The left side is *half the MSD growth rate*; the right side is the *running Green–Kubo integral* of the current autocorrelation. They are the same number: the MSD is the time-integral of the current autocorrelation, so a clean ``\mathrm{Var}(t)\sim 2Dt`` slope and a saturated Green–Kubo integral are one statement. The entire numerical task is to reach the **late-time plateau** ``D(t)\to D`` — precisely where operator entanglement is worst — and since ``z = 2`` is already assumed, a method that biases *toward* diffusion is harmless.

So extracting ``D`` calls for whatever reaches the plateau cheaply: **DAOE** with a ``\gamma \to 0`` extrapolation (how the canonical mixed-field-Ising value ``D \approx 1.40`` was obtained), or **DMT** via thermal-mode decay or a boundary-driven steady state (how the cross-method benchmark reached ``D \approx 1.446`` to ``0.23\%``). Plain TEBD reaches the plateau only at large ``\chi``.

**The logical order is therefore ``z`` first, then ``D``.** Establish the universality class with an unbiased method; only if it is diffusive does a finite ``D`` exist to measure.

## Decision matrix

| You want… | Evolve (trace) | Use | Avoid | Why |
|---|---|---|---|---|
| ``z``, anomalous (``3/2``, ``4``, …) | trace-ful melt | **DMT**, long ``t`` | **DAOE** | the exponent must be unbiased; DAOE forces ``z \to 2`` |
| ``z``, "is it even diffusive?" | melt or correlator | DMT or TEBD | DAOE for the verdict | same |
| ``D`` (diffusive), traceless probe | current / density | **DAOE**, ``\gamma \to 0`` | — | needs the late plateau; the diffusion bias is harmless |
| ``D`` (diffusive), trace-ful probe | thermal density | **DMT** (mode decay / NESS) | — | trace-ful ⇒ DMT's home turf |
| short-time ground truth / cross-check | any | **TEBD** | — | unbiased, exact up to ``\chi`` |

In one sentence: **TEBD is the unbiased but bond-starved reference; DMT is for trace-ful densities and is the right tool for the exponent ``z`` (and for ``D`` via mode decay); DAOE is for traceless operators and the right tool for the constant ``D`` — but never for an anomalous exponent.**

## Knowing when you've converged

Every method here is a controlled approximation, so a transport number is only as good as its error bars. Four checks, all visible in the worked example's output columns:

- **Truncation bar — a conserved quantity must not drift.** Total energy (the melt) or ``\sum_x C(x,t) = \chi`` (the spread) drifting from its ``t = 0`` value *is* the truncation error; keep it far below the signal. In the mixed-field-Ising example DAOE holds this to ``\sim 0.1\%``, where bare TEBD leaks several percent — a direct readout of which method is in control.
- **``\chi``-convergence.** Raise `maxdim` until the transport observable stops moving. DMT and DAOE converge at far smaller ``\chi`` than plain SVD, which is their whole point; if the answer still drifts with ``\chi``, it is not converged.
- **DAOE insensitivity window.** Genuine hydrodynamics is insensitive to ``\ell^\ast`` and ``\gamma`` over a window. The extrapolation must be smooth — ``D`` should rise monotonically and decelerate as ``\gamma \downarrow 0`` (and as ``\ell^\ast`` grows); if it is still moving steeply you are dissipating physics, not just the scrambled front.
- **Finite size and finite time.** Read the exponent or constant in the *plateau* between the fast-mode transient (``t \sim`` a few ``/\Omega``) and the moment the front reaches the chain edges — the example prints a front-at-edges guard for exactly this. For constrained models, also watch the sector-leakage diagnostic (see the [DMT](dmt.md) PXP example).

## Worked illustration: the mixed-field Ising chain

[`examples/dmt/mixed_field_ising_diffusion.jl`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/dmt/mixed_field_ising_diffusion.jl) puts all three on one model — the strongly chaotic mixed-field Ising chain

```math
H = \sum_i \left[\, J\, Z_i Z_{i+1} + g_x X_i + g_z Z_i \,\right],
\qquad J = 1,\; g_x = 1.4,\; g_z = 0.9045,
```

whose energy transport is diffusive with a well-established ``z = 2`` and ``D \approx 1.446``. The script extracts each quantity with the method matched to it:

- **melt → ``z``** uses **DMT** on the trace-ful energy domain wall: the energy transferred across the wall obeys ``\Delta E(t) \sim t^{1/z}``, and the running exponent drifts upward through ``z(t)\approx 1.5`` toward ``2`` while the energy-drift bar stays near ``1\%``.
- **spread → ``D``** uses **DAOE**-assisted TEBD on the *traceless* central energy density ``h_c(t)``: its spatial variance is exactly the ``M_2(t)`` above, so ``\tfrac12\,dM_2/dt`` is the running ``D(t)``. DAOE holds the conserved-energy bar at ``\approx 0.1\%`` (plain TEBD leaks several percent), ``D(t)`` plateaus near ``1.28`` at modest resources, and it climbs toward ``1.446`` as ``\gamma \to 0`` / ``\ell^\ast`` grows / ``t`` lengthens.
- **TEBD** is the unbiased short-time cross-check that both routes are validated against.

Same model, two requirements, two tools — the concrete embodiment of the decision matrix and of the running-``D`` identity above.

## References

- C. David White, M. Zaletel, R. S. K. Mong, G. Refael, [Quantum dynamics of thermalizing systems](https://arxiv.org/abs/1707.01506) — DMT.
- T. Rakovszky, C. W. von Keyserlingk, F. Pollmann, [Dissipation-assisted operator evolution method](https://arxiv.org/abs/2004.05177) — DAOE.
- C. W. von Keyserlingk, F. Pollmann, T. Rakovszky, [Operator backflow and the classical simulation of quantum transport](https://arxiv.org/abs/2111.09904) — why DAOE's bias is exponentially small in ``\ell^\ast``.
- C. Yoo, C. D. White, B. Swingle, [Open-system spin transport and operator weight dissipation in spin chains](https://arxiv.org/abs/2210.06494) — DAOE's controlled bias toward diffusion.
- D. Yi-Thomas, B. Ware, J. D. Sau, C. D. White, [Comparing numerical methods for hydrodynamics in a one-dimensional lattice spin model](https://arxiv.org/abs/2310.06886) — the cross-method ``D \approx 1.446`` benchmark.
