"""
    project!(psi, truncation)

Backend-dispatched in-place projection or truncation entry point used after evolution.
"""
function project!(psi, truncation)
  throw(MethodError(project!, (psi, truncation)))
end

"""
    project!(psi::MPS, trunc)

Truncate a finite OBC MPS in place using standard MPS compression.
"""
function project!(psi::MPS, trunc::BondDimTruncation)
  truncate!(psi; maxdim=trunc.maxdim, cutoff=trunc.cutoff)
  return psi
end
