"""
    project!(psi, truncation)

Backend-dispatched in-place projection or truncation entry point.

# Arguments
- `psi`: Mutable state to project back into the chosen variational manifold.
- `truncation`: Backend-specific truncation or projection settings.

# Returns
- The same `psi` object after in-place mutation.

# Notes
- ScarFinder uses this after each evolution step, which keeps the projection logic
  explicit instead of hiding it inside the evolution backend.
"""
function project!(psi, truncation)
  throw(MethodError(project!, (psi, truncation)))
end

"""
    project!(psi::MPS, trunc)

Truncate a finite `MPS` in place using standard ITensor/ITensorMPS compression.

# Arguments
- `psi`: Finite `MPS` to truncate.
- `trunc`: [`BondDimTruncation`](@ref) object providing `maxdim` and `cutoff`.

# Returns
- The mutated `psi`.

# Notes
- This is the default ScarFinder projection rule for finite-state workflows.
- Compression is delegated to `truncate!`, so orthogonality handling follows the
  underlying ITensor implementation.
"""
function project!(psi::MPS, trunc::BondDimTruncation)
  truncate!(psi; maxdim=trunc.maxdim, cutoff=trunc.cutoff)
  return psi
end
