"""
    _dmt_complement_budget(maxdim, d, radius)

Return the number of complement singular directions `chi'` that fit inside a total bond budget
`maxdim`, given `radius` protected sites per side at local dimension `d`.

`maxdim = chi_preserve + chi_extra` with `chi_preserve = 2 d^(2 radius)`; this is the convention
of arXiv:1902.01859 and of the reference implementations.
"""
function _dmt_complement_budget(maxdim::Integer, d::Integer, radius::Integer)
  return Int(maxdim) - 2 * Int(d)^(2 * Int(radius))
end

"""
    _validate_dmt_budget(psi, maxdim, preserve_diameter)

Throw an `ArgumentError` unless `maxdim` leaves room for the protected block. Called once at the
entry to a DMT bond update, so the failure is immediate rather than mid-sweep.
"""
function _validate_dmt_budget(psi::MPS, maxdim::Integer, preserve_diameter::Integer)
  isodd(preserve_diameter) && preserve_diameter >= 1 ||
    throw(ArgumentError("DMT preserve_diameter must be a positive odd integer, got $(preserve_diameter)"))
  radius = (Int(preserve_diameter) - 1) ÷ 2
  d = local_dimension(siteind(psi, 1))
  floor_value = 2 * d^(2 * radius) + 1
  Int(maxdim) >= floor_value || throw(ArgumentError(
    "DMT requires maxdim >= 2 d^(preserve_diameter - 1) + 1 = $(floor_value) " *
    "for local dimension d = $(d) at preserve_diameter = $(preserve_diameter); got maxdim = $(maxdim). " *
    "maxdim is the total bond dimension, inclusive of the protected block."))
  return nothing
end

"""
    _dmt_protected_sites(bond, nsites, radius)

Return `(left_count, right_count)`: how many sites can actually be protected on each side of
`bond`, clamped by the chain edges.
"""
function _dmt_protected_sites(bond::Integer, nsites::Integer, radius::Integer)
  return (min(Int(radius), Int(bond)), min(Int(radius), Int(nsites) - Int(bond)))
end

"""
    _unit_direction(x)

Return `x / norm(x)`, or an exact zero vector when `x` vanishes.

# Notes
- The identity/trace direction of a *traceless* operator can be exactly zero. Plain `normalize`
  would return `NaN`s there and poison the whole bond; a zero vector instead makes
  [`_dmt_connector`](@ref) report no connector, which is the intended traceless behaviour.
"""
function _unit_direction(x::AbstractVector)
  scale = norm(x)
  scale == 0 && return zeros(ComplexF64, length(x))
  return ComplexF64.(x ./ scale)
end

"""
    _dmt_bond_truncate!(psi, bond; maxdim, cutoff, direction=:R, preserve_diameter=3,
                        truncation=:dense, orthogonalize=true, cache=nothing)

Perform one DMT-preserving bond truncation.

The bond matrix is expressed in a basis whose leading directions span the local-operator
subspace on the sites adjacent to the cut; those `d^(2 radius)` rows and columns, together with
the rank-one trace connector, are carried through **exactly**, and only the doubly-orthogonal
complement is truncated — to `chi' = maxdim - 2 d^(2 radius)` directions, so the reinstated
protected block always fits inside `maxdim` and never has to be clipped.

# Arguments
- `psi`: Operator-space `MPS` to mutate in place.
- `bond`: Bond index to truncate, in `1:length(psi)-1`.

# Keyword Arguments
- `maxdim`: Total post-truncation bond dimension, inclusive of the protected block.
- `cutoff`: Relative cutoff applied in the final refactorization.
- `direction`: `:R` leaves the orthogonality center at `bond + 1`, `:L` at `bond`.
- `preserve_diameter`: Odd; every observable of diameter at most this value is preserved
  exactly. `radius = (preserve_diameter - 1) / 2` sites are protected per side.
- `truncation`: `:dense` or `:random` complement truncation (see [`_truncated_svd`](@ref)).
- `orthogonalize`: Re-gauge so the orthogonality center is at `bond` first. Set `false` only
  when the center is already known to be there.
- `cache`: Optional [`_DMTEnvCache`](@ref) memoizing the identity/trace environments.

# Returns
- The mutated `psi`.

# Notes
- The left protected block is conjugated: the paper's pairing is `M = Q_L^T s Q_R`, a transpose
  rather than an adjoint. Omitting the conjugation is silently correct for a Hermitian operator
  and badly wrong otherwise.
- Both `direction` branches produce the same physical `psi[bond] * psi[bond + 1]`; they differ
  only in which tensor carries the singular values.
"""
function _dmt_bond_truncate!(
  psi::MPS,
  bond::Integer;
  maxdim::Integer,
  cutoff::Real,
  direction::Symbol=:R,
  preserve_diameter::Integer=3,
  truncation::Symbol=:dense,
  orthogonalize::Bool=true,
  cache::Union{Nothing,_DMTEnvCache}=nothing,
)
  direction === :R || direction === :L || throw(ArgumentError("DMT direction must be :R or :L"))
  1 <= bond < length(psi) || throw(ArgumentError("DMT bond must lie in 1:length(psi)-1"))
  _validate_dmt_budget(psi, maxdim, preserve_diameter)
  link = linkind(psi, bond)
  isnothing(link) && return psi
  dim(link) <= maxdim && return psi

  radius = (Int(preserve_diameter) - 1) ÷ 2
  nsites = length(psi)
  left_count, right_count = _dmt_protected_sites(bond, nsites, radius)

  if orthogonalize
    isnothing(cache) ? orthogonalize!(psi, bond) : _orthogonalize_env!(cache, psi, bond)
  end

  left_site = siteind(psi, bond)
  d = local_dimension(left_site)

  # The identity/trace environments depend only on the untouched tensors outside the protected
  # window, so a supplied `cache` memoizes them (see `_DMTEnvCache`); `_DMT_VERIFY_ENVS[]`
  # asserts the memoized value equals the from-scratch rebuild.
  if isnothing(cache)
    left_env = _left_identity_environment(psi, bond - left_count)
    right_env = _right_identity_environment(psi, bond + 1 + right_count)
  else
    left_env = _left_env_at!(cache, psi, bond - left_count)
    right_env = _right_env_at!(cache, psi, bond + 1 + right_count)
    if _DMT_VERIFY_ENVS[]
      _assert_env_matches("left b=$bond", left_env,
        _left_identity_environment(psi, bond - left_count))
      _assert_env_matches("right b=$bond", right_env,
        _right_identity_environment(psi, bond + 1 + right_count))
    end
  end

  previous_link = linkind(psi, bond - 1)
  left_inds = isnothing(previous_link) ? (left_site,) : (previous_link, left_site)
  u, s, v = svd(psi[bond], left_inds)
  left_link = commonind(u, s)
  right_link = commonind(v, s)
  right_block = v * psi[bond + 1]

  # Protected blocks: identity on every site except the `radius` sites adjacent to the cut.
  left_protected = left_env
  for site in (bond - left_count + 1):(bond - 1)
    left_protected *= psi[site]
  end
  left_protected *= u
  right_protected = right_block
  for site in (bond + 2):(bond + right_count)
    right_protected *= psi[site]
  end
  right_protected *= right_env

  left_sites = [siteind(psi, site) for site in (bond - left_count + 1):bond]
  right_sites = [siteind(psi, site) for site in (bond + 1):(bond + right_count)]
  left_combiner = combiner(left_sites...)
  right_combiner = combiner(right_sites...)
  # conj: the paper's pairing is M = Q_L^T s Q_R, a transpose on the left. Omitting this is
  # silently correct for a Hermitian operator and badly wrong otherwise.
  protected_left = conj(matrix(left_protected * left_combiner, left_link, combinedind(left_combiner)))
  protected_right = matrix(right_protected * right_combiner, right_link, combinedind(right_combiner))

  bond_matrix = Diagonal(real.(diag(matrix(s, left_link, right_link))))
  chi = size(bond_matrix, 1)
  ql = _protected_basis(protected_left)
  qr_basis = _protected_basis(protected_right)
  # Column 1 of each protected block is the all-identity multi-index, i.e. the trace direction.
  q0 = _unit_direction(protected_left[:, 1])
  r0 = _unit_direction(protected_right[:, 1])
  a, b, _ = _dmt_connector(bond_matrix, q0, r0)
  ops = _dmt_complement_ops(bond_matrix, a, b, ql, qr_basis)

  budget = max(_dmt_complement_budget(maxdim, d, radius), 1)
  uc, sc, vc = _truncated_svd(ops.mul, ops.adj, chi, budget; mode=truncation)

  # M' = C + QL (QL' B) + BQRc QR' + Uc Sc Vc'  in factored form.
  factor_left = hcat(a, ql, ops.BQRc, uc * Diagonal(sc))
  factor_right = hcat(conj(b), ops.QLtB', qr_basis, vc)
  new_u, new_s, new_v = _dmt_refactor(factor_left, factor_right, Int(maxdim), cutoff)

  # Absorb the singular values on the side the sweep is moving away from, so the orthogonality
  # centre ends up where the next step expects it: at bond + 1 for :R, at bond for :L.
  new_link = Index(length(new_s), "Link,l=$(bond)")
  if direction === :R
    psi[bond] = u * ITensor(new_u, left_link, new_link)
    psi[bond + 1] = ITensor(Diagonal(new_s) * new_v', dag(new_link), right_link) * right_block
  else
    psi[bond] = u * ITensor(new_u * Diagonal(new_s), left_link, new_link)
    psi[bond + 1] = ITensor(Matrix(new_v'), dag(new_link), right_link) * right_block
  end
  isnothing(cache) || _invalidate_env!(cache, bond - left_count, bond + 1 + right_count)
  return psi
end
