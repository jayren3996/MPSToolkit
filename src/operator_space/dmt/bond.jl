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
entry to a DMT step or sweep, so the failure is immediate rather than mid-sweep.

# Notes
- `_validate_dmt_step` runs this *before* the gate is applied, so a rejected budget never
  leaves a partially updated `psi` behind. [`_dmt_bond_truncate!`](@ref) repeats the check as a
  backstop for callers that reach the kernel directly.
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
    _dmt_bond_factorize(psi, bond, left_inds; factorize=:qr)

Split the bond tensor `psi[bond]` into an orthonormal left basis, the bond matrix expressed in
that basis, and the right block, returning
`(left_isometry, bond_matrix, right_block, left_link, right_link)`.

# Arguments
- `psi`: Operator-space `MPS`, already gauged so that `psi[bond]` is the orthogonality centre.
- `bond`: Bond index being truncated.
- `left_inds`: Indices of `psi[bond]` that belong to the left of the cut.

# Keyword Arguments
- `factorize`: `:qr` for `psi[bond] = Q R`, `:svd` for `psi[bond] = U S V'`.

# Returns
- `left_isometry`: `ITensor` carrying `left_inds` and `left_link`, orthonormal in `left_link`.
- `bond_matrix`: `UpperTriangular` (`:qr`) or `Diagonal` (`:svd`) matrix over
  `(left_link, right_link)`.
- `right_block`: `psi[bond + 1]` in the right basis; `right_link` is one of its indices.

# Notes
- DMT needs an orthonormal basis on each side of the cut and the bond matrix expressed in them;
  it never needs the Schmidt form. The two factorizations therefore yield the *same* truncated
  operator: they differ by a unitary on each side, every step of the rule (protected projectors,
  trace connector, truncation of the doubly-projected complement) is covariant under those, and
  [`_dmt_refactor`](@ref) recovers the Schmidt form at the end either way.
- `:qr` is the default because this factorization dominates a bond step: at `chi = 1600, d = 3`
  it is 7.24 s of a 10.16 s step, and QR does it in 5.19 s. The price is that the bond matrix is
  triangular rather than diagonal, so complement products cost `O(chi^2)` instead of `O(chi)`;
  `dev/bench_dmt.jl` measures the net, and `:qr` came out ahead in every cell of a
  `d in 2:4` by `chi in (400, 800, 1600)` sweep, by 1.2x to 1.9x on the whole step.
- `:qr` needs a *square* `R`. A bond tensor with fewer rows than columns is rank deficient by
  construction, and there the rank-revealing SVD is both required and cheap (its cost is set by
  the small row count), so this falls back to `:svd`.
- The `UpperTriangular` wrapper is only applied when `R` really is upper triangular, because the
  wrapper *discards* whatever lies below the diagonal. Wrapping is a pure speed optimization —
  the dense matrix takes the same code path in [`_bond_mul`](@ref) — so a future `qr` that
  returned a permuted factor would slow this down rather than corrupt it.
"""
function _dmt_bond_factorize(psi::MPS, bond::Integer, left_inds; factorize::Symbol=:qr)
  factorize === :qr || factorize === :svd ||
    throw(ArgumentError("DMT factorize must be :qr or :svd, got $(factorize)"))
  if factorize === :qr
    q, r = qr(psi[bond], left_inds)
    left_link = commonind(q, r)
    right_link = commonind(r, psi[bond + 1])
    # `ComplexF64` even for a real-coefficient operator: everything the bond matrix multiplies
    # (`_protected_basis`, the connector, the randomized probes) is complex by construction, and
    # a real triangular factor against complex blocks falls off BLAS onto generic matmul.
    bond_matrix = convert(Matrix{ComplexF64}, matrix(r, left_link, right_link))
    if size(bond_matrix, 1) == size(bond_matrix, 2)
      istriu(bond_matrix) && (bond_matrix = UpperTriangular(bond_matrix))
      return (q, bond_matrix, psi[bond + 1], left_link, right_link)
    end
  end
  u, s, v = svd(psi[bond], left_inds)
  left_link = commonind(u, s)
  right_link = commonind(v, s)
  bond_matrix = Diagonal(real.(diag(matrix(s, left_link, right_link))))
  return (u, bond_matrix, v * psi[bond + 1], left_link, right_link)
end

"""
    _dmt_bond_truncate!(psi, bond; maxdim, cutoff, direction=:R, preserve_diameter=3,
                        truncation=:dense, factorize=:qr, orthogonalize=true, cache=nothing)

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
- `factorize`: `:qr` or `:svd` bond factorization (see [`_dmt_bond_factorize`](@ref)). Both give
  the same truncated operator; `:qr` is the faster, by 1.2x to 1.9x on a whole bond step.
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
  only in which tensor carries the singular values. Under `truncation = :random` this holds only
  to randomized-SVD accuracy, because the two calls draw independent sketches — one of the
  reasons `:dense` is the default (see [`_truncated_svd`](@ref)).
"""
function _dmt_bond_truncate!(
  psi::MPS,
  bond::Integer;
  maxdim::Integer,
  cutoff::Real,
  direction::Symbol=:R,
  preserve_diameter::Integer=3,
  truncation::Symbol=:dense,
  factorize::Symbol=:qr,
  orthogonalize::Bool=true,
  cache::Union{Nothing,_DMTEnvCache}=nothing,
)
  direction === :R || direction === :L || throw(ArgumentError("DMT direction must be :R or :L"))
  # Checked up front as well as inside `_dmt_bond_factorize`, so a typo is rejected even on a
  # bond that is already under budget and short-circuits below.
  factorize === :qr || factorize === :svd ||
    throw(ArgumentError("DMT factorize must be :qr or :svd, got $(factorize)"))
  1 <= bond < length(psi) || throw(ArgumentError("DMT bond must lie in 1:length(psi)-1"))
  # Backstop: `_validate_dmt_step` already rejected an under-floor budget before the gate ran,
  # so this only fires for callers that enter the kernel directly.
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
  left_isometry, bond_matrix, right_block, left_link, right_link =
    _dmt_bond_factorize(psi, bond, left_inds; factorize=factorize)

  # Protected blocks: identity on every site except the `radius` sites adjacent to the cut.
  left_protected = left_env
  for site in (bond - left_count + 1):(bond - 1)
    left_protected *= psi[site]
  end
  left_protected *= left_isometry
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
  # The fused width is the `d^(2 radius)` the budget arithmetic assumes; assert it rather than
  # trust that `left_sites`/`right_sites` still enumerate the protected window. `combiner` fuses
  # in an unspecified order, which is harmless for the span, but column 1 must still be the
  # all-identity multi-index for `q0`/`r0` below -- index 1 of every factor maps to index 1 of
  # any product ordering, and `test_dmt_higher_spin.jl` pins that against an identity cap.
  @assert size(protected_left, 2) == d^(2 * left_count)
  @assert size(protected_right, 2) == d^(2 * right_count)

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
    psi[bond] = left_isometry * ITensor(new_u, left_link, new_link)
    psi[bond + 1] = ITensor(Diagonal(new_s) * new_v', dag(new_link), right_link) * right_block
  else
    psi[bond] = left_isometry * ITensor(new_u * Diagonal(new_s), left_link, new_link)
    psi[bond + 1] = ITensor(Matrix(new_v'), dag(new_link), right_link) * right_block
  end
  isnothing(cache) || _invalidate_env!(cache, bond - left_count, bond + 1 + right_count)
  return psi
end
