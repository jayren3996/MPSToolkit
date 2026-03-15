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

## Energy-Window Cutoff

`energy_cutoff!` provides a standalone CheMPS-style helper for projecting recursion vectors back into a target energy window. The integrated `chebyshev_moments(...; energy_cutoff=true, ...)` path is most relevant when:

- the target spectral window is narrower than the full rescaled interval
- the recursion order is high
- bond-dimension limits make the recursion compression-sensitive

For a worked example, see the Chebyshev notebook links on [Examples](../examples.md).
