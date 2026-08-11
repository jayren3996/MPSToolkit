"""
    _dmt_complement_budget(maxdim, left_count, right_count)

Return the number of complement singular directions `chi'` that fit inside a total bond budget
`maxdim`, given `left_count` and `right_count` protected directions on the two sides of the cut.

`maxdim = chi_preserve + chi_extra` with `chi_preserve = left_count + right_count`; this is the
convention of arXiv:1902.01859 and of the reference implementations. For the default full local
block both counts are `d^(2 radius)`; with `preserve_operators` they are the number of selected
operators per side (see [`_preserved_operator_tensors`](@ref)).

# Notes
- `left_count + right_count` is **not** one too many, and tightening it by one is a bug
  that the suite catches only via the traceless testset. On a trace-carrying bond the reinstated
  rank really is `2 d^(2 radius) - 1`, because [`_dmt_connector`](@ref) leaves
  `B = S - a b^T` annihilating the identity directions (`B r0 = 0`, `q0' B = 0`), costing
  `QL' B` a row and `B QR` a column while the rank-one connector puts one back. But on a
  *traceless* operator -- a transport current -- the connector is disabled, `B = S`, neither
  identity is available, and the reinstated rank saturates `maxdim`. Reserving one fewer then
  makes [`_dmt_refactor`](@ref) clip a protected direction: measured, the trace of a traceless
  operator drifts to 2e-3. The budget is sized for that worst case.
- Jack Kemp's reference implementation does reserve `sdimL + sdimR - 1` (`DMT.h:719-720`), but
  it has no traceless branch to be correct for -- its connector divides by `D[1,1]` unguarded,
  which is ~0 exactly where this kernel switches the connector off. The one-direction
  difference buys support for operators that implementation cannot evolve.
"""
function _dmt_complement_budget(maxdim::Integer, left_count::Integer, right_count::Integer)
  return Int(maxdim) - (Int(left_count) + Int(right_count))
end

"""
    _preserved_count(d, radius, specs)

Return how many directions [`_preserved_operator_tensors`](@ref) protects per side in the bulk,
without building them. Used by [`_validate_dmt_budget`](@ref), which runs before a state is in
hand.

# Notes
- The bulk value is the worst case: a window clipped by a chain edge is narrower and protects
  fewer directions, so a budget that clears the bulk clears every bond.
"""
function _preserved_count(d::Integer, radius::Integer, specs)
  isnothing(specs) && return Int(d)^(2 * Int(radius))
  total = 1                                    # the identity is always protected
  for op in specs
    support = _operator_span(size(op, 1), Int(d))
    support <= Int(radius) || throw(ArgumentError(
      "preserve_operators entry has support $(support) but only $(radius) sites are protected " *
      "per side at preserve_diameter = $(2 * Int(radius) + 1)"))
    total += Int(radius) - support + 1         # every offset that fits in the window
  end
  return total
end

"""
    _validate_dmt_budget(psi, maxdim, preserve_diameter, preserve_operators)

Throw an `ArgumentError` unless `maxdim` leaves room for the protected block. Called once at the
entry to a DMT step or sweep, so the failure is immediate rather than mid-sweep.

# Notes
- `_validate_dmt_step` runs this *before* the gate is applied, so a rejected budget never
  leaves a partially updated `psi` behind. [`_dmt_bond_truncate!`](@ref) repeats the check as a
  backstop for callers that reach the kernel directly.
"""
function _validate_dmt_budget(psi::MPS, maxdim::Integer, preserve_diameter::Integer,
                              preserve_operators=nothing)
  isodd(preserve_diameter) && preserve_diameter >= 1 ||
    throw(ArgumentError("DMT preserve_diameter must be a positive odd integer, got $(preserve_diameter)"))
  radius = (Int(preserve_diameter) - 1) ÷ 2
  d = local_dimension(siteind(psi, 1))
  per_side = _preserved_count(d, radius, preserve_operators)
  floor_value = 2 * per_side + 1
  if Int(maxdim) < floor_value
    detail = isnothing(preserve_operators) ?
      "2 d^(preserve_diameter - 1) + 1 = $(floor_value) for local dimension d = $(d) at " *
      "preserve_diameter = $(preserve_diameter)" :
      "2 * $(per_side) + 1 = $(floor_value) for the $(length(preserve_operators)) selected " *
      "preserve_operators at preserve_diameter = $(preserve_diameter) (each is protected at " *
      "every offset that fits the window, plus the identity)"
    throw(ArgumentError(
      "DMT requires maxdim >= $(detail); got maxdim = $(maxdim). " *
      "maxdim is the total bond dimension, inclusive of the protected block."))
  end
  return nothing
end

"""
    _preserved_operator_tensors(window_sites, d, specs)

Return the operator-space coefficient tensors whose span DMT should protect on one side of a
cut, as `ITensor`s over `window_sites`. The identity is always first.

# Arguments
- `window_sites`: Operator-space site indices of the protected window, in chain order.
- `d`: Local Hilbert space dimension.
- `specs`: `nothing` for the full local-operator block (`d^(2 length(window_sites))` directions),
  or a collection of dense physical matrices, each acting on `1 <= support <= length(window_sites)`
  contiguous sites. Each is placed at **every** offset that fits inside the window.

# Returns
- A `Vector{ITensor}`, element 1 the identity. Its length is what the budget arithmetic must
  reserve, so callers take it from `length(...)` rather than recomputing it.

# Notes
- Selecting a proper subset is a **weaker** guarantee than the full block: only observables in
  the selected span survive truncation exactly. It exists because the full block costs
  `d^(2 radius)` per side, which at `d = 4, preserve_diameter = 5` is 256 per side and a floor of
  513 — protecting energy density and a charge instead costs of order ten. This is the
  `presAll_ = false` path of Jack Kemp's reference implementation (`DMT.h:159-206`), which is
  marked experimental there; the construction here is independent.
- Each operator is padded to the full window with identity and expanded on
  [`_operator_basis_operators`](@ref), rather than tracked per site. The window is `radius` sites
  wide, so the padded matrix is `d^radius x d^radius` — 16 x 16 at the largest configuration
  this unlocks — and the cost is irrelevant next to the bond step.
- The flat expansion enumerates the site-1 label slowest (it is built by `kron`), while `reshape`
  makes the first dimension fastest, so the site indices are attached in reverse. Getting this
  backwards would silently mirror a multi-site operator across the window; the full-block
  equivalence test is what pins it.
"""
function _preserved_operator_tensors(window_sites, d::Integer, specs)
  width = length(window_sites)
  local_dim = Int(d)
  basis = _operator_basis_operators(width, local_dim)
  full_dim = local_dim^width
  coefficient_tensor(matrix) = ITensor(
    reshape(ComplexF64[tr(adjoint(element) * matrix) for element in basis],
            ntuple(_ -> local_dim^2, width)),
    reverse(collect(window_sites))...)

  identity_matrix = Matrix{ComplexF64}(I, full_dim, full_dim)
  isnothing(specs) && return [coefficient_tensor(element) for element in basis]

  tensors = ITensor[coefficient_tensor(identity_matrix)]
  for op in specs
    size(op, 1) == size(op, 2) ||
      throw(ArgumentError("preserve_operators entries must be square, got $(size(op))"))
    support = _operator_span(size(op, 1), local_dim)
    support <= width || throw(ArgumentError(
      "preserve_operators entry has support $(support) but only $(width) sites are protected " *
      "per side at this preserve_diameter; widen preserve_diameter or shrink the operator"))
    dense = Matrix{ComplexF64}(op)
    for offset in 0:(width - support)
      before = Matrix{ComplexF64}(I, local_dim^offset, local_dim^offset)
      after_dim = local_dim^(width - support - offset)
      after = Matrix{ComplexF64}(I, after_dim, after_dim)
      push!(tensors, coefficient_tensor(kron(before, dense, after)))
    end
  end
  return tensors
end

"""
    _protected_columns(protected_tensor, link, operator_tensors)

Contract the protected half-chain tensor against each preserved operator, returning the
`chi x length(operator_tensors)` protected block whose column 1 is the trace direction.

# Notes
- Contracting operator by operator rather than fusing the window with a `combiner` keeps this
  independent of the order `combiner` happens to fuse in, and is what lets the selected and
  full-block paths share one code path.
"""
function _protected_columns(protected_tensor::ITensor, link::Index, operator_tensors)
  columns = [vector(protected_tensor * op, link) for op in operator_tensors]
  return reduce(hcat, columns)
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

Return `x / norm(x)`, or an exact zero vector when `x` vanishes, in `x`'s own (floating) element
type.

# Notes
- The identity/trace direction of a *traceless* operator can be exactly zero. Plain `normalize`
  would return `NaN`s there and poison the whole bond; a zero vector instead makes
  [`_dmt_connector`](@ref) report no connector, which is the intended traceless behaviour.
- No element type is imposed here. The result only ever feeds [`_dmt_connector`](@ref), which
  widens to the bond step's own `T`, so coercing twice would be redundant — and the zero branch
  has to agree with the division branch, which is what `float` is for.
"""
function _unit_direction(x::AbstractVector)
  scale = norm(x)
  scale == 0 && return zeros(float(eltype(x)), length(x))
  return x ./ scale
end

"""
    _mps_eltype(psi)

Return the single element type one DMT bond step should run in: the promotion over the element
types of **every** tensor in `psi`.

# Notes
- The scan covers the whole chain rather than just `psi[bond]` because a bond step consumes
  rather more than the bond tensor: `psi[bond + 1]`, the interior protected sites, and both
  identity/trace environments, which contract over every tensor in the chain. A single complex
  tensor anywhere therefore makes `protected_left`/`protected_right` complex, and the `hcat`s
  that assemble the refactorization would re-promote a partially converted factor with no error
  (`hcat(::Matrix{Float64}, ::Matrix{ComplexF64})` is `Matrix{ComplexF64}`). The scan costs
  ~100 ns against the `O(chi^3)` work of the step it types.
- The result is a *promotion*, so every conversion inside the kernel widens and an imaginary
  component can never be dropped. It is also gauge invariant: `orthogonalize!`, `qr` and `svd`
  preserve element types, so scanning before the re-gauge below gives the same answer as after.
- `float` maps an `Int`-valued tensor (legal in ITensors) onto `Float64`; the DMT algebra divides.
- `NDTensors.EmptyNumber`, the element type of a tensor that has never been assigned, has to be
  screened out **by name**: it subtypes `Real` (verified on ITensors 0.9.30 --
  `supertype(EmptyNumber) === Real`), so no subtype test against `Number`, `Real` or `Complex`
  excludes it, and `float(EmptyNumber)` is `Float64` -- which would type the whole step off a
  chain carrying no data at all. Nothing can be inferred there, so it falls back to the
  `ComplexF64` this kernel used unconditionally before the type was threaded.
"""
function _mps_eltype(psi::MPS)
  T = eltype(psi[1])
  for n in 2:length(psi)
    T = promote_type(T, eltype(psi[n]))
  end
  T === ITensors.NDTensors.EmptyNumber && return ComplexF64
  T <: Union{Real,Complex} || return ComplexF64
  return float(T)
end

"""
    _dmt_bond_factorize(psi, bond, left_inds, T=ComplexF64; factorize=:qr)

Split the bond tensor `psi[bond]` into an orthonormal left basis, the bond matrix expressed in
that basis, and the right block, returning
`(left_isometry, bond_matrix, right_block, left_link, right_link)`.

# Arguments
- `psi`: Operator-space `MPS`, already gauged so that `psi[bond]` is the orthogonality centre.
- `bond`: Bond index being truncated.
- `left_inds`: Indices of `psi[bond]` that belong to the left of the cut.
- `T`: Element type for the `:qr` bond matrix, normally [`_mps_eltype`](@ref) of the chain.
  Defaults to `ComplexF64`, the type this kernel used unconditionally.

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
- The `:svd` bond matrix is a **real** `Diagonal` whatever `T` is, and deliberately does not
  follow it. Singular values are real, and `_bond_mul(::Diagonal, x)` is elementwise rather than
  BLAS, so a real `Diagonal` against complex blocks costs nothing. It is also why `T` must be
  threaded in rather than read back off `bond_matrix`: `eltype` of this factor is `Float64` even
  for a complex state, and typing the step off it would throw on the first protected block.
"""
function _dmt_bond_factorize(psi::MPS, bond::Integer, left_inds, ::Type{T}=ComplexF64;
                             factorize::Symbol=:qr) where {T}
  factorize === :qr || factorize === :svd ||
    throw(ArgumentError("DMT factorize must be :qr or :svd, got $(factorize)"))
  if factorize === :qr
    q, r = qr(psi[bond], left_inds)
    left_link = commonind(q, r)
    right_link = commonind(r, psi[bond + 1])
    # One element type for everything the bond matrix meets. The conversion has to be
    # all-or-nothing -- `T` is a single promotion over the whole chain, so it widens here and
    # nowhere narrows -- because `NDTensors` answers a mismatched element type by materializing a
    # promoted copy of the larger operand (184 MB per bond at chi = 800), and a mixed contraction
    # falls off BLAS onto generic matmul. Both costs are paid silently.
    bond_matrix = convert(Matrix{T}, matrix(r, left_link, right_link))
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
    _dmt_complement_keep(complement_values, cutoff)

Return how many of the complement singular values `complement_values` to keep under a relative
`cutoff`, measured against the **complement's own** leading value.

# Notes
- The complement `D = P_L^perp B P_R^perp` is the only subspace DMT is allowed to discard, so it
  is the only one a cutoff may be applied to. Applying the cutoff to the reassembled bond matrix
  instead — which is what [`_dmt_refactor`](@ref) does when handed a physics cutoff — voids the
  preservation guarantee outright, because the top-`k` singular directions of the reassembled
  matrix need not contain the protected row and column spaces. That failure is not hypothetical
  and not graceful: in operator space the identity direction dominates `sigma_1`, so on a state
  carrying a signal `eps` above an infinite-temperature background *every* cutoff larger than
  `eps` discards the whole protected block. Measured at `d = 2, eps = 1e-6, maxdim = 14`, a
  `cutoff = 1e-6` collapsed the bond to `chi = 1` — the pure identity — for a diameter-3
  preservation error of `0.98`, against `1.5e-7` at `cutoff = 0`. The `d = 3` cell of the same
  sweep gave `chi = 1` and `0.76`, and even at `eps = 1e-3` the guarantee was gone by
  `cutoff = 1e-4`.
- Measuring against `complement_values[1]` rather than the bond's own scale is what makes the
  cutoff scale-free inside the subspace it governs. It cannot then interact with the identity
  background at all, which is what the manual's "invariant to the cutoff itself" claims and what
  was previously true only in the regime where the signal was a few percent of the norm.
"""
function _dmt_complement_keep(complement_values::AbstractVector, cutoff::Real)
  (cutoff <= 0 || isempty(complement_values)) && return length(complement_values)
  return count(>(cutoff * complement_values[1]), complement_values)
end

"""
    _dmt_refactor_tolerance(factor)

Return the relative tolerance for the final repair SVD: a *numerical rank* threshold, not a
physics cutoff.

# Notes
- `LinearAlgebra.rank`'s convention, `min(size(A)...) * eps`. Directions below it are null to the
  accuracy of the QR/SVD that just reassembled the bond, so dropping them changes nothing
  measurable — while dropping anything above it could clip a protected direction (see
  [`_dmt_complement_keep`](@ref)). The traceless bond needs this: there the connector is disabled
  and `hcat`'s first column is exactly zero, so without a rank cut the bond would carry one
  guaranteed-null direction.
"""
_dmt_refactor_tolerance(factor::AbstractMatrix) = minimum(size(factor)) * eps(Float64)

"""
    _dmt_bond_solve(T, bond_matrix, protected_left, protected_right, maxdim, cutoff,
                    truncation)

Apply the DMT rule to a bond matrix already expressed in orthonormal bases on both sides of the
cut, returning the truncated bond matrix in SVD form `(U, S, V)`.

The protected row and column spaces and the rank-one trace connector are reinstated **exactly**;
only the doubly-orthogonal complement `D = P_L^perp B P_R^perp` is truncated, to
`maxdim` minus the protected widths.

# Arguments
- `T`: Element type to work in, a promotion over the whole chain (see [`_mps_eltype`](@ref)).
- `bond_matrix`: Bond matrix over `(left_link, right_link)` from [`_dmt_bond_factorize`](@ref).
- `protected_left`, `protected_right`: `chi x d^(2 n)` protected blocks, conjugated on the left,
  with column 1 the all-identity multi-index (the trace direction).
- `maxdim`: Total post-truncation bond dimension, inclusive of the protected block.
- `cutoff`: Relative cutoff on the **complement** singular values only (see
  [`_dmt_complement_keep`](@ref)). The final repair SVD runs at a numerical-rank tolerance
  instead, so no cutoff can clip a protected direction.
- `truncation`: `:dense` or `:random` (see [`_truncated_svd`](@ref)).

# Returns
- `(U, S, V)` with `U * Diagonal(S) * V'` the truncated bond matrix.

# Notes
- Contains **no ITensors**, only matrices. That buys three things: the DMT algebra is unit-testable
  without building a state, `T` specializes the whole rule rather than each helper separately, and
  the [`_dmt_complement_ops`](@ref) closures are created *inside* a body where their captures are
  concretely typed — so the `mul` that `_truncated_svd`'s power iteration calls `2 power + 1`
  times per step is a static call rather than a dynamic dispatch.
- `_protected_basis`/`_dmt_connector` are handed `T` explicitly rather than left to infer it,
  because their inputs can be narrower than `T`: the `:svd` bond matrix is a real `Diagonal` even
  for a complex state, and `protected_left` contracts only sites `1:bond` (so it stays real when
  the sole complex tensor sits to the right of the cut). Convert them all and the two `hcat`s
  below are homogeneous; convert only some and `hcat` re-promotes the rest silently, which is the
  failure mode `_DMT_VERIFY_ELTYPE` exists to catch.
"""
function _dmt_bond_solve(::Type{T}, bond_matrix::AbstractMatrix, protected_left::AbstractMatrix,
                         protected_right::AbstractMatrix,
                         maxdim::Integer, cutoff::Real, truncation::Symbol) where {T}
  chi = size(bond_matrix, 1)
  ql = _protected_basis(protected_left, T)
  qr_basis = _protected_basis(protected_right, T)
  # Column 1 of each protected block is the all-identity multi-index, i.e. the trace direction.
  q0 = _unit_direction(protected_left[:, 1])
  r0 = _unit_direction(protected_right[:, 1])
  a, b, _ = _dmt_connector(bond_matrix, q0, r0, T)
  ops = _dmt_complement_ops(bond_matrix, a, b, ql, qr_basis)

  budget = max(_dmt_complement_budget(maxdim, size(protected_left, 2),
                                      size(protected_right, 2)), 1)
  uc, sc, vc = _truncated_svd(ops.mul, ops.adj, chi, budget, T; mode=truncation)
  # `cutoff` binds *here*, on the complement, and nowhere else. See `_dmt_complement_keep`: the
  # protected block and the trace connector are reinstated exactly, so a cutoff applied after
  # reassembly would be free to clip them.
  complement = _dmt_complement_keep(sc, cutoff)
  uc, sc, vc = uc[:, 1:complement], sc[1:complement], vc[:, 1:complement]

  # M' = C + QL (QL' B) + BQRc QR' + Uc Sc Vc'  in factored form.
  factor_left = hcat(a, ql, ops.BQRc, uc * Diagonal(sc))
  factor_right = hcat(conj(b), ops.QLtB', qr_basis, vc)
  if _DMT_VERIFY_ELTYPE[]
    eltype(factor_left) === T && eltype(factor_right) === T || error(
      "DMT eltype verify: the step is threaded at T = $(T) but the refactorization factors are " *
      "$(eltype(factor_left)) (left) / $(eltype(factor_right)) (right). `hcat` re-promotes a " *
      "partially converted factor with no error, so this is a silently widened step.")
  end
  return _dmt_refactor(factor_left, factor_right, Int(maxdim),
                       _dmt_refactor_tolerance(factor_left))
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
- `cutoff`: Relative cutoff on the discarded complement, measured against the leading complement
  singular value. It cannot reach the protected block or the trace connector, so the preservation
  guarantee holds at every `cutoff` (see [`_dmt_complement_keep`](@ref)).
- `direction`: `:R` leaves the orthogonality center at `bond + 1`, `:L` at `bond`.
- `preserve_diameter`: Odd; every observable of diameter at most this value is preserved
  exactly. `radius = (preserve_diameter - 1) / 2` sites are protected per side.
- `preserve_operators`: `nothing` (default) protects the full `d^(2 radius)` local block per
  side. A collection of dense physical matrices instead protects only their span, at every
  offset that fits the window, plus the identity — a **weaker** guarantee, and the only way to
  reach `preserve_diameter = 5` at `d = 4`, whose full-block floor is 513. See
  [`_preserved_operator_tensors`](@ref).
- `truncation`: `:dense` or `:random` complement truncation (see [`_truncated_svd`](@ref)).
- `factorize`: `:qr` or `:svd` bond factorization (see [`_dmt_bond_factorize`](@ref)). Both give
  the same truncated operator; `:qr` is the faster, by 1.2x to 1.9x on a whole bond step.
- `orthogonalize`: Re-gauge so the orthogonality center is at `bond` first. Set `false` only
  when the center is already known to be there.
- `cache`: Optional [`_DMTEnvCache`](@ref) memoizing the identity/trace environments.

# Returns
- The mutated `psi`.

# Notes
- The step runs in the element type of the state, [`_mps_eltype`](@ref) of the whole chain, not
  in an unconditional `ComplexF64`. Every operator-space basis element is Hermitian, so a
  Hermitian operator — any density matrix — has real coefficients and a real chain stays real all
  the way through, at roughly half the memory. The conversion is all-or-nothing on purpose: a
  mixed chain makes `NDTensors` materialize a promoted copy of the larger operand at every
  contraction across the seam, and `orthogonalize!` then spreads the mixedness.
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
  preserve_operators=nothing,
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
  _validate_dmt_budget(psi, maxdim, preserve_diameter, preserve_operators)
  link = linkind(psi, bond)
  isnothing(link) && return psi
  dim(link) <= maxdim && return psi

  # One element type for the whole step, taken from the state rather than forced to `ComplexF64`
  # (see `_mps_eltype`). Scanned here, before the re-gauge below, because gauging preserves
  # element types and the promotion over the chain is therefore the same either side of it.
  elt = _mps_eltype(psi)
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
    _dmt_bond_factorize(psi, bond, left_inds, elt; factorize=factorize)

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
  # conj: the paper's pairing is M = Q_L^T s Q_R, a transpose on the left. Omitting this is
  # silently correct for a Hermitian operator and badly wrong otherwise.
  if isnothing(preserve_operators)
    # Default: the whole local block, fused in one contraction. Going operator by operator here
    # would cost `d^(2 radius)` contractions per side -- 256 at `d = 4, radius = 2` -- for a
    # basis the combiner spans in one.
    left_combiner = combiner(left_sites...)
    right_combiner = combiner(right_sites...)
    protected_left = conj(matrix(left_protected * left_combiner, left_link,
                                 combinedind(left_combiner)))
    protected_right = matrix(right_protected * right_combiner, right_link,
                             combinedind(right_combiner))
    # The fused width is the `d^(2 radius)` the budget arithmetic assumes; assert it rather than
    # trust that `left_sites`/`right_sites` still enumerate the protected window. `combiner` fuses
    # in an unspecified order, which is harmless for the span, but column 1 must still be the
    # all-identity multi-index for `q0`/`r0` below -- index 1 of every factor maps to index 1 of
    # any product ordering, and `test_dmt_higher_spin.jl` pins that against an identity cap.
    @assert size(protected_left, 2) == d^(2 * left_count)
    @assert size(protected_right, 2) == d^(2 * right_count)
  else
    protected_left = conj(_protected_columns(left_protected, left_link,
      _preserved_operator_tensors(left_sites, d, preserve_operators)))
    protected_right = _protected_columns(right_protected, right_link,
      _preserved_operator_tensors(right_sites, d, preserve_operators))
  end

  new_u, new_s, new_v = _dmt_bond_solve(elt, bond_matrix, protected_left, protected_right,
                                        maxdim, cutoff, truncation)

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
