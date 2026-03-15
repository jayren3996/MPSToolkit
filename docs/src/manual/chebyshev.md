# Chebyshev

The Chebyshev module provides the low-level pieces needed to build spectral calculations:

- `ChebyshevRescaling`
- `SpectralFunction`
- `chebyshev_moments`
- `energy_cutoff!`
- `jackson_damping`
- `jackson_kernel`
- `reconstruct_chebyshev`
- `spectral_function`

## Rescaling

`chebyshev_moments` assumes the input Hamiltonian has already been rescaled into the Chebyshev window `[-1, 1]`.

In practice this means:

1. choose a physical center and halfwidth
2. shift and rescale the Hamiltonian
3. generate moments in the rescaled variable
4. reconstruct the physical spectrum using `spectral_function`

## Minimal Workflow

```julia
using MPSToolkit

moments = chebyshev_moments(
  h_rescaled,
  psi;
  order=80,
  maxdim=64,
  cutoff=1e-12,
)

spectrum = spectral_function(moments; center=center, halfwidth=halfwidth)
```

The `energy_cutoff=true` path is most useful when the relevant spectral window is narrow compared with the full rescaled interval and the recursion is compression-sensitive.

## Energy-Window Cutoff

`energy_cutoff!` provides a standalone CheMPS-style helper for projecting recursion vectors back into a target energy window. The integrated `chebyshev_moments(...; energy_cutoff=true, ...)` path is most relevant when:

- the target spectral window is narrower than the full rescaled interval
- the recursion order is high
- bond-dimension limits make the recursion compression-sensitive

For a worked example, see the Chebyshev notebook links on [Examples](../examples.md).

```julia
moments = chebyshev_moments(
  h_rescaled,
  psi;
  order=120,
  maxdim=32,
  cutoff=1e-12,
  energy_cutoff=true,
  energy_cutoff_sweeps=4,
  krylovdim=20,
  window=0.35,
)
```

## API

```@docs
ChebyshevRescaling
SpectralFunction
chebyshev_moments
energy_cutoff!
jackson_damping
jackson_kernel
reconstruct_chebyshev
spectral_function
```

## References

- Andreas Holzner, Andreas Weichselbaum, Ian P. McCulloch, Ulrich Schollwock, and Jan von Delft, [Chebyshev matrix product state approach for spectral functions](https://arxiv.org/abs/1101.5895)
- Jad C. Halimeh, Fabian Kolley, and Ian P. McCulloch, [Chebyshev matrix product state approach for time evolution](https://arxiv.org/abs/1507.01226)
