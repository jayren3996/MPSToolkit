# DAOE

The DAOE utilities build diagonal projector MPOs that damp operator components according to a notion of operator size rather than bond dimension. They are intended for operator-space workflows where the goal is to preserve hydrodynamic transport while suppressing entanglement growth.

`MPSToolkit.jl` currently provides two related projectors:

- `pauli_daoe_projector` for the original DAOE construction, where size is measured by Pauli weight
- `fdaoe_projector` for the fermionic generalization, where size is measured by fermionic weight

## DAOE Vs FDAOE

The difference is the notion of locality that is being damped:

- DAOE is the natural choice for spin chains written in a local Pauli basis, where operator complexity is tracked by how many non-identity Pauli factors appear in a string.
- FDAOE is the natural choice for fermionic systems, where Jordan-Wigner strings make plain Pauli weight a poor proxy for physical operator size.

So the practical rule of thumb is:

- use `pauli_daoe_projector` for spin models and Pauli-basis transport calculations
- use `fdaoe_projector` when fermionic structure matters and Pauli-weight damping would distort the relevant operator hierarchy

## DAOE Example

```julia
using MPSToolkit

sites = pauli_siteinds(8)
projector = pauli_daoe_projector(sites; lstar=3, gamma=0.15)
```

## FDAOE Example

```julia
using MPSToolkit

sites = pauli_siteinds(8)
projector = fdaoe_projector(sites; wstar=3, gamma=0.15)
```

Because the projectors are returned as MPOs, you can combine them with your own operator-space pipelines instead of being forced into a single packaged workflow.

## API

```@docs
pauli_daoe_projector
fdaoe_projector
```

## Examples

For a simple script-level entry point, see [examples/daoe/projectors.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/daoe/projectors.jl).

## References

- Tibor Rakovszky, C. W. von Keyserlingk, and Frank Pollmann, [Dissipation-assisted operator evolution method for capturing hydrodynamic transport](https://arxiv.org/abs/2004.05177)
- Yongchan Yoo, Christopher David White, and Brian Swingle, [Open-system spin transport and operator weight dissipation in spin chains](https://arxiv.org/abs/2210.06494)
- En-Jui Kuo, Brayden Ware, Peter Lunts, Mohammad Hafezi, and Christopher David White, [Energy diffusion in weakly interacting chains with fermionic dissipation-assisted operator evolution](https://arxiv.org/abs/2311.17148)
