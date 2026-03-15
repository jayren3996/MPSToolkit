# DAOE

The DAOE utilities build diagonal projector MPOs in the Pauli basis. They are intended for operator-space workflows where Pauli strings should be damped according to a structured weight rule rather than by generic bond-dimension truncation.

`MPSToolkit.jl` currently provides two related projectors:

- `pauli_daoe_projector` for Pauli-weight damping beyond a cutoff `lstar`
- `fdaoe_projector` for the fermionic-weight variant used in FDAOE-style constructions

## Minimal Example

```julia
using MPSToolkit

sites = pauli_siteinds(8)
projector = pauli_daoe_projector(sites; lstar=3, gamma=0.15)
```

Because the projectors are returned as MPOs, you can combine them with your own operator-space pipelines instead of being forced into a single packaged workflow.

## API

```@docs
pauli_daoe_projector
fdaoe_projector
```

## Examples

For a simple script-level entry point, see [examples/daoe/projectors.jl](https://github.com/jayren3996/MPSToolkit/blob/main/examples/daoe/projectors.jl).
