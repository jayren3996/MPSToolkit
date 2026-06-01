# Chebyshev (CheMPS)

*Chebyshev / kernel-polynomial expansion of spectral functions, evaluated with finite MPS, plus an energy-window cutoff for narrow target bands.*

## Background

A great many dynamical quantities reduce to a **spectral function** of the form

```math
A(\omega) = \sum_k |\langle k | \psi \rangle|^2 \, \delta(\omega - E_k),
```

where ``|\psi\rangle = O|0\rangle`` is a probe state obtained by acting with some operator ``O`` on a reference state ``|0\rangle`` (often the ground state), and ``|k\rangle``, ``E_k`` are eigenstates and eigenvalues of a Hamiltonian ``H``. Examples include local spectral densities, dynamical structure factors, and optical conductivities. The object ``A(\omega)`` is a weighted density of states: each eigenstate contributes a stick of weight ``|\langle k|\psi\rangle|^2`` at energy ``E_k``, and the total weight ``\int A(\omega)\,d\omega = \langle\psi|\psi\rangle = \|O|0\rangle\|^2``.

Computing ``A(\omega)`` by full diagonalization is exponentially expensive. The **kernel polynomial method** (KPM) sidesteps diagonalization by expanding ``A(\omega)`` in Chebyshev polynomials and computing only a finite number of expansion coefficients — the *moments* — each of which is a single inner product. For an MPS probe state and an MPO Hamiltonian, every moment is obtained from a matrix–product contraction, so the whole spectral function follows from one polynomial recursion on tensor networks. This is the Chebyshev MPS (CheMPS) approach of Holzner *et al.*

### Chebyshev polynomials

The Chebyshev polynomials of the first kind ``T_n(x)`` are defined on ``x \in [-1, 1]`` by ``T_n(\cos\theta) = \cos(n\theta)``. They satisfy the three-term recursion

```math
T_0(x) = 1, \qquad T_1(x) = x, \qquad T_{n+1}(x) = 2x\,T_n(x) - T_{n-1}(x),
```

and are orthogonal with respect to the weight ``w(x) = 1/(\pi\sqrt{1-x^2})``:

```math
\int_{-1}^{1} \frac{T_m(x)\,T_n(x)}{\pi\sqrt{1 - x^2}}\,dx =
\begin{cases}
1 & m = n = 0 \\
\tfrac{1}{2}\,\delta_{mn} & m = n > 0 .
\end{cases}
```

Expanding a function ``f(x)`` on ``[-1,1]`` in this basis gives

```math
f(x) = \frac{1}{\pi\sqrt{1 - x^2}}\left( \mu_0 + 2\sum_{n=1}^{\infty} \mu_n\,T_n(x) \right),
\qquad
\mu_n = \int_{-1}^{1} f(x)\,T_n(x)\,dx .
```

### Rescaling to ``[-1, 1]``

The expansion above only converges on ``[-1, 1]``, but a physical Hamiltonian has eigenvalues on some interval ``[E_{\min}, E_{\max}]``. We therefore introduce the **rescaled Hamiltonian**

```math
\tilde H = \frac{H - a}{b}, \qquad
a = \frac{E_{\max} + E_{\min}}{2}, \qquad
b = \frac{E_{\max} - E_{\min}}{2(1 - \epsilon)},
```

with center ``a`` and halfwidth ``b``. The small safety margin ``\epsilon`` pulls the extreme eigenvalues to ``\pm(1-\epsilon)`` so that no spectral weight reaches the band edge ``x = \pm 1``, where the Chebyshev weight ``1/(\pi\sqrt{1-x^2})`` diverges. A physical frequency ``\omega`` maps to the rescaled coordinate

```math
x = \frac{\omega - a}{b}.
```

Setting ``f(x) = A(a + b x)``, the spectral function in rescaled coordinates has the Chebyshev moments

```math
\mu_n = \int_{-1}^{1} A(a + b x)\,T_n(x)\,dx
      = \langle\psi|\,T_n(\tilde H)\,|\psi\rangle,
```

where the second equality follows from inserting the spectral resolution of ``\tilde H`` and using ``T_n(\tilde H)|k\rangle = T_n((E_k - a)/b)|k\rangle``. So **the moments are just expectation values of Chebyshev polynomials of the rescaled Hamiltonian in the probe state**.

### Computing moments by recursion

Applying the Chebyshev recursion to the *vector* ``|\psi\rangle`` rather than to scalars gives a sequence of states

```math
|t_0\rangle = |\psi\rangle, \qquad
|t_1\rangle = \tilde H\,|\psi\rangle, \qquad
|t_{n+1}\rangle = 2\tilde H\,|t_n\rangle - |t_{n-1}\rangle,
```

so that ``|t_n\rangle = T_n(\tilde H)\,|\psi\rangle`` and

```math
\mu_n = \langle\psi | t_n\rangle .
```

Each step applies the MPO ``2\tilde H`` to an MPS, subtracts the previous MPS, and truncates the result back to a bounded bond dimension. The cost per moment is therefore one MPO–MPS product plus one MPS addition, both controlled by `maxdim` and `cutoff`.

### Gibbs oscillations and the Jackson kernel

Truncating the series at ``N`` moments is equivalent to convolving ``A(\omega)`` with the Dirichlet kernel, which produces ringing (Gibbs oscillations) and can drive the reconstruction negative near sharp features. The KPM fix is to multiply each moment by a damping factor ``g_n`` before summing:

```math
A(x) \approx \frac{1}{\pi\sqrt{1 - x^2}}\left( g_0\mu_0 + 2\sum_{n=1}^{N-1} g_n\,\mu_n\,T_n(x) \right).
```

The **Jackson kernel** is the standard near-optimal choice. It corresponds to convolution with an almost-Gaussian of nearly minimal width that keeps the reconstruction non-negative:

```math
g_n = \frac{(N - n + 1)\cos\!\frac{n\pi}{N+1} + \sin\!\frac{n\pi}{N+1}\,\cot\!\frac{\pi}{N+1}}{N + 1}.
```

The effective resolution scales as ``\sim b\,\pi/N``, so resolving fine structure requires both a high moment order ``N`` and a halfwidth ``b`` that is not much larger than the relevant band.

### Reconstruction

Given the moments and a kernel, the reconstructed density at a rescaled coordinate ``x`` is evaluated with the same three-term recursion run on the scalars ``T_n(x)``, weighted by ``g_n\mu_n`` and divided by the Chebyshev measure. Mapping ``x`` back to physical frequency ``\omega = a + bx`` and dividing by ``b`` (the Jacobian of the rescaling) yields ``A(\omega)``.

## Rescaling

The affine map between physical frequency and the Chebyshev coordinate is stored in a [`ChebyshevRescaling`](@ref) value, which holds a `center` ``a`` and a positive `halfwidth` ``b``:

```julia
using MPSToolkit

rescaling = ChebyshevRescaling(0.0, 2.5)   # center = 0.0, halfwidth = 2.5
```

You rarely build one by hand. `chebyshev_rescaling(emin, emax; margin=0.05)` constructs it from spectral bounds, mapping the extreme eigenvalues to ``\pm(1-\text{margin})`` so that no weight lands on the band edge. It returns a `ChebyshevRescaling` with `center = (emax + emin)/2` and `halfwidth = (emax - emin) / (2(1 - margin))`.

To rescale an MPO Hamiltonian directly, use `rescale_hamiltonian`:

```julia
H_rescaled, rescaling = rescale_hamiltonian(H, emin, emax; margin=0.05, cutoff=1e-14)
```

This forms ``\tilde H = (H - a)/b`` by subtracting `a * MPO(firstsiteinds(H), "Id")` (truncated with `cutoff`) and dividing by `b`. It returns the rescaled MPO together with the `ChebyshevRescaling` you then reuse in [`spectral_function`](@ref). The bounds `emin`, `emax` may be exact or estimates; conservative estimates (e.g. from Gershgorin bounds or a sum of local operator norms) are fine as long as they bracket the true spectrum, since the `margin` keeps the mapped spectrum strictly inside ``(-1, 1)``.

!!! note "Input must be pre-rescaled"
    [`chebyshev_moments`](@ref) assumes its Hamiltonian argument **already** has its spectrum in ``[-1, 1]``. Pass the `H_rescaled` returned by `rescale_hamiltonian` (or your own ``(H - a)/b``), never the bare physical `H`. If you skip rescaling, the recursion ``|t_{n+1}\rangle = 2\tilde H|t_n\rangle - |t_{n-1}\rangle`` diverges because ``|T_n(x)|`` grows without bound for ``|x| > 1``.

If you would rather not subtract the shift into a single rescaled MPO, you can rescale by halfwidth only and keep `center = 0.0` when the band is already symmetric about zero, exactly as the two-peak example does: `h_rescaled = h_mpo / halfwidth`. The reconstruction still needs the same `center`/`halfwidth` to map back.

## Computing moments

[`chebyshev_moments`](@ref) runs the recursion and returns the moment vector ``\mu_0, \mu_1, \dots, \mu_{N-1}`` as a dense `Vector{Float64}`, stored in Julia's `1`-based convention so that `moments[n + 1]` holds ``\mu_n``:

```julia
moments = chebyshev_moments(
    H_rescaled,
    psi;
    order=80,        # number of moments N
    maxdim=64,       # bond-dimension cap during MPO·MPS and MPS addition
    cutoff=1e-12,    # truncation cutoff during the same operations
)
```

The two required pieces are `order` (how many moments, i.e. the truncation length ``N``) and the compression controls `maxdim` and `cutoff`, which bound the growth of the Chebyshev vectors ``|t_n\rangle``. Higher `order` sharpens the resolution; a larger `maxdim` reduces truncation error in the recursion at the cost of memory and time. The function checks that `H` and `psi` share site indices before starting.

### `normalize_initial` and absolute spectral weight

By default `normalize_initial=true`, so the probe MPS is normalized before the recursion begins. This is what you want when you only care about the **shape** of the spectrum: the reconstruction then integrates to one. The total spectral weight, however, is physically ``\int A(\omega)\,d\omega = \langle\psi|\psi\rangle = \|O|0\rangle\|^2``, which is lost under normalization. To recover **absolute intensity** you have two options:

- pass `normalize_initial=false` so the recursion seeds with the unnormalized probe and the moments carry the full weight (`moments[1]` then equals ``\langle\psi|\psi\rangle`` rather than `1.0`); or
- keep the default normalization and multiply the reconstructed ``A(\omega)`` by ``\|O|0\rangle\|^2`` afterwards.

Both give the same physical spectral function. The first keeps everything inside the moment vector; the second is convenient when the probe norm is known independently.

## Reconstruction

Once you have the moments, [`spectral_function`](@ref) packages them into a callable [`SpectralFunction`](@ref):

```julia
spectrum = spectral_function(moments; center=rescaling.center, halfwidth=rescaling.halfwidth, kernel=:jackson)

A_at_omega = spectrum(0.3)                 # evaluate at one physical frequency
A_on_grid  = spectrum(collect(-2.0:0.01:2.0))   # vectorized over a frequency grid
```

The keyword `kernel` accepts `:jackson` (the default Jackson kernel of the right length), `nothing` (the raw, undamped truncated series), or an explicit weight vector of length `order`. Calling the returned object maps each physical frequency ``\omega`` to ``x = (\omega - \text{center})/\text{halfwidth}``; frequencies that fall outside the open band (``|x| \ge 1``, including the edges ``x = \pm 1``) return `0.0`, so keep the spectrum strictly inside the band through the rescaling margin.

Under the hood, evaluation calls [`reconstruct_chebyshev`](@ref), which you can also use directly on the rescaled coordinate:

```julia
kernel = jackson_kernel(length(moments))
rho_x  = reconstruct_chebyshev(0.12, moments; kernel=kernel)   # x ∈ (-1, 1)
```

`reconstruct_chebyshev(x, moments; kernel=nothing)` already includes the Chebyshev weight factor ``1/(\pi\sqrt{1-x^2})``; with `kernel=nothing` it sums the raw moments without damping. Note that `reconstruct_chebyshev` works in the **rescaled** coordinate ``x`` and does not divide by the halfwidth, whereas a `SpectralFunction` works in physical ``\omega`` and applies the ``1/b`` Jacobian for you.

The damping factors come from [`jackson_kernel`](@ref) and [`jackson_damping`](@ref):

- `jackson_damping(n, order)` returns the single Jackson factor ``g_n`` for zero-based index ``n``.
- `jackson_kernel(order)` returns the full weight vector, with `kernel[n + 1]` holding ``g_n``.

You can build a custom kernel and pass it to either `spectral_function` or `reconstruct_chebyshev`; it only has to match `order` in length.

A [`SpectralFunction`](@ref) stores its `moments`, the resolved `kernel` (already a numeric vector, or `nothing`), and the `ChebyshevRescaling` used to convert frequencies, so it is fully self-describing and reusable.

## Energy-window cutoff (CheMPS)

When the interesting part of the spectrum occupies only a narrow slice of the full rescaled band — say a few low-lying excitations — most of the Chebyshev vector's bond dimension is spent representing high-energy weight that you do not care about. The energy-window cutoff is a CheMPS-style filter that projects each recursion vector back toward a target rescaled window ``[-\text{window}, \text{window}] \subset [-1, 1]``, freeing bond dimension for the part of the spectrum you keep.

[`energy_cutoff!`](@ref) is the standalone projector, mutating an MPS in place:

```julia
energy_cutoff!(psi, H_rescaled; sweeps=5, krylovdim=30, window=0.4, tol=1e-12, verbose=false)
```

It sweeps left–right then right–left, and at each site builds a small local Lanczos/Krylov approximation (dimension `krylovdim`) of the single-site effective Hamiltonian from the MPO environment, then removes the components whose local energies fall outside ``[-\text{window}, \text{window}]``. The `tol` argument is a stopping monitor on a heuristic per-sweep removed-weight estimate, not a calibrated error bound, and `verbose=true` prints per-sweep diagnostics.

The same projection is wired into the recursion through the `energy_cutoff=true` path of [`chebyshev_moments`](@ref), which applies `energy_cutoff!` after each Chebyshev step:

```julia
moments = chebyshev_moments(
    H_rescaled,
    psi;
    order=120,
    maxdim=32,
    cutoff=1e-12,
    energy_cutoff=true,
    energy_cutoff_sweeps=4,
    krylovdim=20,
    window=0.35,
    energy_cutoff_tol=1e-12,
    energy_cutoff_verbose=false,
)
```

The cutoff keyword arguments map onto `energy_cutoff!`: `energy_cutoff_sweeps` → `sweeps`, `krylovdim` → `krylovdim`, `window` → `window`, `energy_cutoff_tol` → `tol`, `energy_cutoff_verbose` → `verbose`.

This path helps most when **the target window is narrow relative to the full rescaled interval, the moment order is high, and the bond dimension is tightly capped** — exactly the regime where an unfiltered recursion would lose in-window accuracy to truncation. It is not free: each step now also pays for the cutoff sweeps.

!!! warning "Local, not global, filtering"
    `energy_cutoff!` filters each site against its *local* effective problem, not the global Hamiltonian spectrum. Out-of-window weight that is delocalized across entangled bonds is only weakly removed, so a strongly entangling state can retain such weight even after many sweeps. Treat the cutoff as a compression aid for narrow windows, and always cross-check against a higher-`maxdim` run on the same window before trusting it.

## Worked example

The following is fully self-contained: it builds a small transverse-field Ising (TFIM) MPO with `OpSum`/`MPO`, rescales it with the package's `rescale_hamiltonian` using a conservative analytic spectral bound, constructs a simple product-state superposition as the probe, runs the recursion, and reconstructs the spectrum. No external diagonalization is needed.

```julia
using MPSToolkit
using ITensors
using ITensorMPS

# --- model: transverse-field Ising on N sites, H = -J Σ Sz Sz - g Σ Sx ---
N = 6
J = 1.0
g = 0.8
sites = siteinds("S=1/2", N)

H = MPO(let opsum = OpSum()
    for j in 1:(N - 1)
        opsum += -J, "Sz", j, "Sz", j + 1
    end
    for j in 1:N
        opsum += -g, "Sx", j
    end
    opsum
end, sites)

# --- conservative spectral bound (sum of local operator norms, spin-1/2 ⇒ ‖Sα‖ = 1/2) ---
# |E| ≤ J·(N-1)·(1/2)² + g·N·(1/2)
ebound = J * (N - 1) * 0.25 + g * N * 0.5
H_rescaled, rescaling = rescale_hamiltonian(H, -ebound, ebound; margin=0.05)

# --- probe state: a superposition of two product states (a local-operator-like seed) ---
psi_a = productMPS(sites, ["Up", "Up", "Dn", "Dn", "Up", "Dn"])
psi_b = productMPS(sites, ["Dn", "Up", "Dn", "Dn", "Up", "Dn"])
psi = add(psi_a, psi_b; maxdim=8, cutoff=1e-14)
normalize!(psi)

# --- Chebyshev moments on the rescaled Hamiltonian ---
order = 80
moments = chebyshev_moments(H_rescaled, psi; order=order, maxdim=64, cutoff=1e-12)

@assert isapprox(moments[1], 1.0; atol=1e-10)   # μ_0 = ⟨ψ|ψ⟩ = 1 after normalization

# --- reconstruct the physical spectral function with Jackson damping ---
spectrum = spectral_function(
    moments;
    center=rescaling.center,
    halfwidth=rescaling.halfwidth,
    kernel=:jackson,
)

omegas = collect(range(-ebound, ebound; length=400))
A = spectrum(omegas)

println("μ_0 = ", moments[1])
println("max A(ω) = ", maximum(A), " at ω = ", omegas[argmax(A)])
```

Expected behavior: with `normalize_initial` left at its default, ``\mu_0 = \langle\psi|\psi\rangle = 1`` (the assertion holds). `A(ω)` is non-negative across the grid because the Jackson kernel suppresses Gibbs ringing, and it concentrates near the energies of the eigenstates that overlap the two-product-state probe. Raising `order` sharpens the peaks; the spectrum integrates (approximately) to ``\langle\psi|\psi\rangle = 1`` because the probe was normalized — multiply by ``\|O|0\rangle\|^2`` (or set `normalize_initial=false`) to recover absolute weight.

## API

```@docs
ChebyshevRescaling
chebyshev_rescaling
rescale_hamiltonian
SpectralFunction
chebyshev_moments
energy_cutoff!
jackson_damping
jackson_kernel
reconstruct_chebyshev
spectral_function
```

## Examples

- [`examples/chebyshev/two_peak_spectrum.ipynb`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/chebyshev/two_peak_spectrum.ipynb)
- [`examples/chebyshev/tfim_local_spectrum.ipynb`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/chebyshev/tfim_local_spectrum.ipynb)
- [`examples/chebyshev/energy_cutoff_comparison.ipynb`](https://github.com/jayren3996/MPSToolkit/blob/main/examples/chebyshev/energy_cutoff_comparison.ipynb)

## References

- Andreas Holzner, Andreas Weichselbaum, Ian P. McCulloch, Ulrich Schollwöck, and Jan von Delft, [Chebyshev matrix product state approach for spectral functions](https://arxiv.org/abs/1101.5895)
- Jad C. Halimeh, Fabian Kolley, and Ian P. McCulloch, [Chebyshev matrix product state approach for time evolution](https://arxiv.org/abs/1507.01226)
- Alexander Weiße, Gerhard Wellein, Andreas Alvermann, and Holger Fehske, [The kernel polynomial method](https://arxiv.org/abs/cond-mat/0504627)
