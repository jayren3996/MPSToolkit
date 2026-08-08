# Higher-Spin DMT Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make DMT correct and efficient for local Hilbert space dimension `d > 2`, by replacing the truncation kernel with one that actually delivers DMT's preservation guarantee and generalizing the operator-space layer away from spin-1/2.

**Architecture:** A generic onsite operator basis (generalized Gell-Mann, identity first) supplied through a provider function; an `operator_*` API layer with the existing `pauli_*` names kept as `d = 2` wrappers; and a basis-free DMT bond truncation that protects the `d^(2n)` local-operator subspace on each side exactly, sizing the complement budget as `chi' = maxdim - 2 d^(2n)` so nothing protected is ever clipped.

**Tech Stack:** Julia 1.10+, ITensors 0.9.25, ITensorMPS 0.3.44/0.4, LinearAlgebra, KrylovKit. Tests with `Test`. Interop target: EDKit.jl `spin(...; D = d)`.

**Spec:** `docs/superpowers/specs/2026-08-08-higher-spin-dmt-design.md` — read it before starting. It contains the measured evidence behind every design choice here.

**Deviation from the spec's file layout, deliberate:** the spec sketches a five-way split of
`src/operator_space/dmt.jl`. This plan creates only `dmt/lowrank.jl` and `dmt/bond.jl` — the two
units with genuinely separate responsibilities and their own tests — and leaves the options type,
the environment cache and the sweep driver in `dmt.jl`, which drops to roughly 400 lines once the
superseded kernel is deleted. Splitting those three further would be churn without a testing or
comprehension benefit.

## Global Constraints

- Branch: `feat/higher-spin-dmt`. Baseline on this branch is fully green; keep it green after every task.
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'` from the repo root. Single-file runs: `julia --project=. test/test_dmt.jl` (test files are self-contained and `using MPSToolkit`).
- Julia compat floor is `1.10`; do not use syntax newer than that.
- Every `operator_*` function must have a `pauli_*` wrapper that is **exactly** the `d = 2` case, and existing `pauli_*` call sites must keep working unchanged.
- Local basis element 1 is always `I / sqrt(d)`. The DMT kernel depends on this and on Hilbert-Schmidt orthonormality, and on nothing else about the basis.
- `maxdim` is the **total** post-truncation bond dimension, inclusive of the protected block: `chi' = maxdim - 2 d^(2n)`, floor `maxdim >= 2 d^(2n) + 1`, with `n = (preserve_diameter - 1) / 2`.
- The left protected block uses `conj(...)` — the paper's pairing is `M = Q_L^T s Q_R`, a transpose. Omitting it is silently correct for Hermitian operators and ~50% wrong otherwise.
- Docstrings on every exported symbol. `test/test_docstrings.jl` enforces this; run it if you add exports.
- Commit after every task with the message given in the task's final step.

---

### Task 1: Generic onsite operator basis

**Files:**
- Create: `src/bases/operator_basis.jl`
- Modify: `src/MPSToolkit.jl` (include + export)
- Create: `test/test_operator_basis.jl`
- Modify: `test/runtests.jl` (add the include)

**Interfaces:**
- Consumes: nothing.
- Produces: `operator_basis_matrices(d::Integer) -> Vector{Matrix{ComplexF64}}` (length `d^2`, element 1 is `I/sqrt(d)`); `local_dimension(x) -> Int` accepting an `Int` site dimension or an `Index`.

- [ ] **Step 1: Write the failing test**

Create `test/test_operator_basis.jl`:

```julia
using ITensors
using LinearAlgebra
using MPSToolkit
using Test

@testset "generic onsite operator basis" begin
  @testset "d = 2 reproduces the normalized Pauli basis exactly" begin
    basis = operator_basis_matrices(2)
    expected = [m / sqrt(2) for m in values(pauli_matrices())]   # (I, X, Y, Z) / sqrt(2)
    @test length(basis) == 4
    for k in 1:4
      @test basis[k] ≈ expected[k] atol = 1e-15
    end
  end

  @testset "orthonormal, Hermitian, identity first" begin
    for d in 2:5
      basis = operator_basis_matrices(d)
      @test length(basis) == d^2
      @test basis[1] ≈ Matrix{ComplexF64}(I, d, d) / sqrt(d) atol = 1e-14
      gram = [tr(basis[i]' * basis[j]) for i in eachindex(basis), j in eachindex(basis)]
      @test gram ≈ Matrix{ComplexF64}(I, d^2, d^2) atol = 1e-12
      for m in basis
        @test m ≈ m' atol = 1e-14
      end
      # every element other than the identity is traceless
      for k in 2:length(basis)
        @test abs(tr(basis[k])) < 1e-14
      end
    end
  end

  @testset "the cache returns the same matrices, not fresh copies" begin
    @test operator_basis_matrices(3) === operator_basis_matrices(3)
  end

  @testset "local_dimension inverts the site dimension" begin
    @test local_dimension(4) == 2
    @test local_dimension(9) == 3
    @test local_dimension(Index(16, "OperatorSpace,n=1")) == 4
    @test_throws ArgumentError local_dimension(8)     # not a perfect square
    @test_throws ArgumentError local_dimension(1)     # d = 1 is not a spin space
    @test_throws ArgumentError operator_basis_matrices(1)
  end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_operator_basis.jl`
Expected: FAIL with `UndefVarError: operator_basis_matrices not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/bases/operator_basis.jl`:

```julia
const _OPERATOR_BASIS_CACHE = Dict{Int,Vector{Matrix{ComplexF64}}}()

"""
    operator_basis_matrices(d)

Return the normalized generalized Gell-Mann basis for a local Hilbert space of dimension `d`.

# Arguments
- `d`: Local Hilbert space dimension (`2S + 1` for spin `S`).

# Returns
- A cached `Vector` of `d^2` dense `d x d` matrices that is Hermitian, traceless apart from
  element 1, and orthonormal under the Hilbert-Schmidt inner product `tr(A' * B)`. Element 1 is
  the normalized identity `I / sqrt(d)`.

# Notes
- The ordering is: the identity; then, for each pair `j < k`, the real and imaginary
  off-diagonal generators; then the `d - 1` diagonal generators. For `d = 2` this reproduces
  `(I, X, Y, Z) / sqrt(2)` exactly, so the `pauli_*` helpers are the `d = 2` case of the
  `operator_*` ones.
- The DMT kernel relies on exactly two properties of this basis: orthonormality, and the
  identity sitting at index 1. Any basis with those properties can be substituted; Hermiticity
  is not required by the algorithm.

# Examples
```jldoctest
julia> length(operator_basis_matrices(3))
9
```
"""
function operator_basis_matrices(d::Integer)
  local_dim = Int(d)
  local_dim >= 2 || throw(ArgumentError("local Hilbert space dimension must be at least 2, got $(local_dim)"))
  return get!(_OPERATOR_BASIS_CACHE, local_dim) do
    _build_operator_basis(local_dim)
  end
end

function _build_operator_basis(d::Int)
  basis = Vector{Matrix{ComplexF64}}()
  push!(basis, Matrix{ComplexF64}(I, d, d) / sqrt(d))
  for j in 1:d, k in (j + 1):d
    symmetric = zeros(ComplexF64, d, d)
    symmetric[j, k] = 1 / sqrt(2)
    symmetric[k, j] = 1 / sqrt(2)
    push!(basis, symmetric)
    antisymmetric = zeros(ComplexF64, d, d)
    antisymmetric[j, k] = -im / sqrt(2)
    antisymmetric[k, j] = im / sqrt(2)
    push!(basis, antisymmetric)
  end
  for level in 1:(d - 1)
    entries = zeros(ComplexF64, d)
    for j in 1:level
      entries[j] = 1
    end
    entries[level + 1] = -level
    push!(basis, Matrix(Diagonal(entries)) / sqrt(level * (level + 1)))
  end
  return basis
end

"""
    local_dimension(site_dimension)
    local_dimension(site::Index)

Return the physical local dimension `d` encoded by an operator-space site of dimension `d^2`.

# Arguments
- `site_dimension`: Operator-space site dimension, or an `Index` carrying it.

# Returns
- The integer `d` with `d^2 == site_dimension`.

# Notes
- Operator-space sites always have a perfect-square dimension because they carry a basis for
  the `d x d` onsite operator space. This is how every `operator_*` helper recovers `d` without
  being told.
"""
function local_dimension(site_dimension::Integer)
  squared = Int(site_dimension)
  d = isqrt(squared)
  d * d == squared || throw(ArgumentError("operator-space site dimension $(squared) is not a perfect square"))
  d >= 2 || throw(ArgumentError("operator-space site dimension $(squared) implies local dimension $(d) < 2"))
  return d
end

local_dimension(site::Index) = local_dimension(dim(site))
```

- [ ] **Step 4: Wire into the module**

In `src/MPSToolkit.jl`, add the include immediately after `include("bases/pauli.jl")` (line 22):

```julia
include("bases/operator_basis.jl")
```

In the `Bases` submodule (line 94-97), extend both the `using` and `export` lists with `operator_basis_matrices, local_dimension`. In the top-level exports, extend the `export pauli_matrices, pauli_basis, pauli_components` line to also export `operator_basis_matrices, local_dimension`.

In `test/runtests.jl`, add `include("test_operator_basis.jl")` alongside the other includes.

- [ ] **Step 5: Run tests to verify they pass**

Run: `julia --project=. test/test_operator_basis.jl`
Expected: PASS, all testsets.

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS — including `test_docstrings.jl`, which requires docstrings on the two new exports.

- [ ] **Step 6: Commit**

```bash
git add src/bases/operator_basis.jl src/MPSToolkit.jl test/test_operator_basis.jl test/runtests.jl
git commit -m "feat(bases): generalized Gell-Mann onsite operator basis for any local dimension"
```

---

### Task 2: Generic site indices and state builders

**Files:**
- Create: `src/operator_space/states.jl`
- Modify: `src/operator_space/daoe.jl` (remove `pauli_siteinds`, lines 1-19)
- Modify: `src/operator_space/helpers.jl` (replace `pauli_basis_state`, `pauli_total_sz_state`, `pauli_domain_wall_state`, `_pauli_basis_label` with wrappers)
- Modify: `src/MPSToolkit.jl`
- Modify: `test/test_operator_space.jl`

**Interfaces:**
- Consumes: `operator_basis_matrices(d)`, `local_dimension(x)` from Task 1.
- Produces:
  - `operator_siteinds(nsites; d=2, tagprefix="OperatorSpace") -> Vector{Index}`
  - `operator_basis_state(sites, labels; coefficient=1.0) -> MPS`
  - `operator_product_state(sites, ops; coefficient=1.0) -> MPS`
  - `operator_local_sum_state(sites, op, coeffs) -> MPS` (bond dimension 2)
  - `_operator_basis_label(label, d) -> Int`
  - `_operator_coefficients(op::AbstractMatrix, d) -> Vector{ComplexF64}` (length `d^2`, entry `mu` is `tr(P_mu' * op)`)

- [ ] **Step 1: Write the failing test**

Append to `test/test_operator_space.jl`:

```julia
@testset "generic operator-space state builders" begin
  @testset "pauli_* wrappers are exactly the d = 2 case" begin
    @test dim(first(operator_siteinds(3; d = 2))) == 4
    @test dim(first(operator_siteinds(3; d = 3))) == 9

    sites = pauli_siteinds(4)
    a = pauli_basis_state(sites, ["I", "Z", "X", "I"])
    b = operator_basis_state(sites, [1, 4, 2, 1])
    @test inner(a, b) ≈ 1.0 atol = 1e-14

    total_via_wrapper = pauli_total_sz_state(sites)
    total_via_generic = operator_local_sum_state(
      sites, pauli_matrices().Z / sqrt(2), fill(2.0^(4 / 2 - 1), 4))
    @test inner(total_via_wrapper, total_via_generic) ≈ inner(total_via_wrapper, total_via_wrapper) atol = 1e-12

    dw_wrapper = pauli_domain_wall_state(sites; kink = 2)
    dw_generic = operator_local_sum_state(
      sites, pauli_matrices().Z / sqrt(2), [j <= 2 ? -1.0 : 1.0 for j in 1:4])
    @test inner(dw_wrapper, dw_generic) ≈ inner(dw_wrapper, dw_wrapper) atol = 1e-12
  end

  @testset "operator_product_state decomposes dense local matrices" begin
    # A spin-1 S^z on site 2 of a 3-site chain, identity elsewhere.
    d = 3
    sites = operator_siteinds(3; d = d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    identity_matrix = Matrix{ComplexF64}(I, d, d)
    state = operator_product_state(sites, [identity_matrix, sz, identity_matrix])
    # Amplitude on basis label mu of site 2, with label 1 (the normalized identity) on the
    # others: tr(P_mu' * S^z) times tr(P_1' * I) = sqrt(d) on each of the two spectator sites.
    basis = operator_basis_matrices(d)
    for mu in 1:d^2
      probe = operator_basis_state(sites, [1, mu, 1])
      @test inner(probe, state) ≈ tr(basis[mu]' * sz) * d atol = 1e-12
    end
  end

  @testset "operator_local_sum_state builds a bond-dimension-2 sum" begin
    d = 3
    sites = operator_siteinds(4; d = d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    coeffs = [1.0, -2.0, 3.0, -4.0]
    state = operator_local_sum_state(sites, sz, coeffs)
    @test all(dim(linkind(state, b)) <= 2 for b in 1:3)
    basis = operator_basis_matrices(d)
    # Coefficient of "S^z on site j, identity elsewhere" must be coeffs[j] * tr(P_mu' * S^z).
    for j in 1:4, mu in 2:d^2
      labels = [k == j ? mu : 1 for k in 1:4]
      probe = operator_basis_state(sites, labels)
      @test inner(probe, state) ≈ coeffs[j] * tr(basis[mu]' * sz) atol = 1e-12
    end
  end

  @testset "label validation is dimension aware" begin
    @test_throws ArgumentError operator_basis_state(operator_siteinds(2; d = 3), [10, 1])
    # Pauli letter labels are only meaningful at d = 2.
    @test_throws ArgumentError operator_basis_state(operator_siteinds(2; d = 3), ["X", "I"])
  end
end
```

If the amplitude convention in the second testset turns out to differ by a power of `sqrt(d)`
from what `operator_product_state` produces, fix the **test** to match `operator_trace`, not the
other way round: `operator_trace` is pinned against dense linear algebra in Task 4 and is the
authority on the normalization.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_operator_space.jl`
Expected: FAIL with `UndefVarError: operator_siteinds not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/operator_space/states.jl`:

```julia
"""
    operator_siteinds(nsites; d=2, tagprefix="OperatorSpace")

Construct site indices for a vectorized operator-space `MPS` on `nsites` sites of local
Hilbert space dimension `d`.

# Arguments
- `nsites`: Number of operator-space sites.

# Keyword Arguments
- `d`: Local Hilbert space dimension (`2S + 1` for spin `S`). The generated indices have
  dimension `d^2`.
- `tagprefix`: Prefix used when naming the generated `Index` tags.

# Returns
- A vector of length `nsites` containing dimension-`d^2` `Index` objects.

# Notes
- The basis ordering on each site is [`operator_basis_matrices`](@ref)`(d)`, whose first
  element is the normalized identity. Downstream helpers recover `d` from the site dimension
  via [`local_dimension`](@ref), so it never has to be passed again.
"""
function operator_siteinds(nsites::Integer; d::Integer=2, tagprefix::AbstractString="OperatorSpace")
  nsites >= 1 || throw(ArgumentError("number of operator-space sites must be positive"))
  local_dim = Int(d)
  local_dim >= 2 || throw(ArgumentError("local Hilbert space dimension must be at least 2"))
  return [Index(local_dim^2, "$(tagprefix),n=$(n)") for n in 1:Int(nsites)]
end

"""
    _operator_basis_label(label, d)

Normalize one local operator-basis label into its integer index in the
[`operator_basis_matrices`](@ref) ordering. Pauli letters are accepted only at `d == 2`.
"""
function _operator_basis_label(label::Integer, d::Integer)
  1 <= label <= Int(d)^2 || throw(ArgumentError("operator-basis labels must lie in 1:$(Int(d)^2) for local dimension $(d)"))
  return Int(label)
end

function _operator_basis_label(label::Symbol, d::Integer)
  return _operator_basis_label(String(label), d)
end

function _operator_basis_label(label::AbstractString, d::Integer)
  normalized = uppercase(strip(label))
  normalized == "I" && return 1
  Int(d) == 2 || throw(ArgumentError("letter operator-basis labels are only defined for local dimension 2; use an integer label in 1:$(Int(d)^2)"))
  normalized == "X" && return 2
  normalized == "Y" && return 3
  normalized == "Z" && return 4
  throw(ArgumentError("unsupported operator-basis label: $(label)"))
end

"""
    _operator_coefficients(op, d)

Return the coefficient vector of a dense `d x d` operator in the normalized basis, i.e. the
entries `tr(P_mu' * op)` for `P_mu = operator_basis_matrices(d)[mu]`.
"""
function _operator_coefficients(op::AbstractMatrix, d::Integer)
  local_dim = Int(d)
  size(op) == (local_dim, local_dim) ||
    throw(ArgumentError("local operator must be $(local_dim) x $(local_dim), got $(size(op))"))
  basis = operator_basis_matrices(local_dim)
  dense = Matrix{ComplexF64}(op)
  return ComplexF64[tr(adjoint(matrix) * dense) for matrix in basis]
end

"""
    operator_basis_state(sites, labels; coefficient=1.0)

Build a product `MPS` in operator space selecting one basis element per site.

# Arguments
- `sites`: Operator-space site indices, typically from [`operator_siteinds`](@ref).
- `labels`: One local basis label per site: an integer in `1:d^2`, or `"I"` (any `d`), or
  `"X"`, `"Y"`, `"Z"` when `d == 2`.

# Keyword Arguments
- `coefficient`: Overall scalar prefactor stored on the first tensor.

# Returns
- A product-state `MPS` in operator space.
"""
function operator_basis_state(sites, labels::AbstractVector; coefficient::Number=1.0)
  length(sites) == length(labels) || throw(ArgumentError("operator-basis labels must have one entry per site"))
  tensors = ITensor[]
  for (index, (site, label)) in enumerate(zip(sites, labels))
    d = local_dimension(site)
    tensor = ITensor(site)
    tensor[site => _operator_basis_label(label, d)] = index == 1 ? coefficient : 1.0
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

"""
    operator_product_state(sites, ops; coefficient=1.0)

Build the operator-space `MPS` for the tensor product `⊗_j O_j` of dense local operators.

# Arguments
- `sites`: Operator-space site indices.
- `ops`: One dense `d x d` matrix per site.

# Keyword Arguments
- `coefficient`: Overall scalar prefactor stored on the first tensor.

# Returns
- A bond-dimension-1 `MPS` whose amplitude on the basis product `alpha` is
  `prod_j tr(P_{alpha_j}' * O_j)`.

# Notes
- Use this rather than [`operator_basis_state`](@ref) whenever the operator you want is not a
  single basis element. At `d = 3` the physical `S^z = diag(1, 0, -1)` is **not** proportional
  to any single Gell-Mann matrix, so a spin-1 `S^z` string must be built this way.
"""
function operator_product_state(sites, ops; coefficient::Number=1.0)
  length(sites) == length(ops) || throw(ArgumentError("operator_product_state needs one local operator per site"))
  tensors = ITensor[]
  for (index, (site, op)) in enumerate(zip(sites, ops))
    d = local_dimension(site)
    coefficients = _operator_coefficients(op, d)
    tensor = ITensor(ComplexF64, site)
    scale = index == 1 ? ComplexF64(coefficient) : one(ComplexF64)
    for mu in eachindex(coefficients)
      iszero(coefficients[mu]) && continue
      tensor[site => mu] = scale * coefficients[mu]
    end
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

"""
    operator_local_sum_state(sites, op, coeffs)

Build the operator-space `MPS` representing the local-density sum `sum_j c_j O_j`, where `O_j`
is the dense local operator `op` placed on site `j` with identity elsewhere.

# Arguments
- `sites`: Operator-space site indices.
- `op`: Dense `d x d` local operator.
- `coeffs`: One coefficient `c_j` per site.

# Returns
- An `MPS` of bond dimension 2 (link state 1 = "operator not yet placed", state 2 = "placed,
  identity henceforth").

# Notes
- This is the transport source builder: uniform `coeffs` give a total charge such as
  `sum_j S^z_j`, and a sign-split profile gives the domain-wall operator whose melting sets the
  dynamical exponent. `pauli_total_sz_state` and `pauli_domain_wall_state` are the `d = 2` cases.
- The identity coefficient inserted on unoccupied sites is `1 / sqrt(d)` times `sqrt(d)`, i.e.
  the amplitude `1` on basis element 1, matching the normalized-basis convention used by
  [`operator_trace`](@ref).
"""
function operator_local_sum_state(sites, op::AbstractMatrix, coeffs::AbstractVector)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("operator_local_sum_state requires at least one site"))
  length(coeffs) == nsites || throw(ArgumentError("operator_local_sum_state needs one coefficient per site"))
  d = local_dimension(first(sites))
  all(local_dimension(site) == d for site in sites) ||
    throw(ArgumentError("operator_local_sum_state requires a uniform local dimension"))
  weights = _operator_coefficients(op, d)

  function placed!(tensor, indices, site, scale)
    for mu in eachindex(weights)
      iszero(weights[mu]) && continue
      tensor[(indices..., site => mu)...] = scale * weights[mu]
    end
    return tensor
  end

  if nsites == 1
    tensor = ITensor(ComplexF64, sites[1])
    placed!(tensor, (), sites[1], ComplexF64(coeffs[1]))
    return MPS([tensor])
  end

  left_link = Index(2, "OperatorStateLink,n=1")
  first_tensor = ITensor(ComplexF64, sites[1], left_link)
  first_tensor[sites[1] => 1, left_link => 1] = 1.0
  placed!(first_tensor, (left_link => 2,), sites[1], ComplexF64(coeffs[1]))
  tensors = ITensor[first_tensor]

  for j in 2:(nsites - 1)
    right_link = Index(2, "OperatorStateLink,n=$(j)")
    tensor = ITensor(ComplexF64, dag(left_link), sites[j], right_link)
    tensor[dag(left_link) => 1, sites[j] => 1, right_link => 1] = 1.0
    tensor[dag(left_link) => 2, sites[j] => 1, right_link => 2] = 1.0
    placed!(tensor, (dag(left_link) => 1, right_link => 2), sites[j], ComplexF64(coeffs[j]))
    push!(tensors, tensor)
    left_link = right_link
  end

  last_tensor = ITensor(ComplexF64, dag(left_link), sites[nsites])
  last_tensor[dag(left_link) => 2, sites[nsites] => 1] = 1.0
  placed!(last_tensor, (dag(left_link) => 1,), sites[nsites], ComplexF64(coeffs[nsites]))
  push!(tensors, last_tensor)
  return MPS(tensors)
end
```

- [ ] **Step 4: Replace the spin-1/2 originals with wrappers**

Delete `pauli_siteinds` from `src/operator_space/daoe.jl` (lines 1-19) and delete
`pauli_basis_state`, `pauli_total_sz_state`, `pauli_domain_wall_state`, and the three
`_pauli_basis_label` methods from `src/operator_space/helpers.jl`. Add to
`src/operator_space/states.jl`:

```julia
"""
    pauli_siteinds(nsites; tagprefix="PauliSpace")

Spin-1/2 case of [`operator_siteinds`](@ref): dimension-4 sites in the `(I, X, Y, Z)` ordering.
"""
function pauli_siteinds(nsites::Integer; tagprefix::AbstractString="PauliSpace")
  return operator_siteinds(nsites; d=2, tagprefix=tagprefix)
end

"""
    pauli_basis_state(sites, labels; coefficient=1.0)

Spin-1/2 case of [`operator_basis_state`](@ref), with labels in the `(I, X, Y, Z)` ordering.
"""
function pauli_basis_state(sites, labels::AbstractVector; coefficient::Number=1.0)
  return operator_basis_state(sites, labels; coefficient=coefficient)
end

"""
    pauli_total_sz_state(sites; coefficient=nothing)

Spin-1/2 case of [`operator_local_sum_state`](@ref) for the uniform total-`S^z` operator.
The default per-site coefficient is `2^(N / 2 - 1)`, matching the normalized Pauli convention.
"""
function pauli_total_sz_state(sites; coefficient::Union{Nothing,Number}=nothing)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("pauli_total_sz_state requires at least one site"))
  weight = isnothing(coefficient) ? 2.0^(nsites / 2 - 1) : coefficient
  return operator_local_sum_state(sites, pauli_matrices().Z / sqrt(2), fill(weight, nsites))
end

"""
    pauli_domain_wall_state(sites; kink=length(sites) ÷ 2, coefficient=1.0)

Spin-1/2 case of [`operator_local_sum_state`](@ref) for the signed domain-wall operator
`D = sum_j sign(j - kink) * sigma^z_j`, the infinite-temperature transport source.
`D` is traceless, so measure its profile with `normalize=false`.
"""
function pauli_domain_wall_state(sites; kink::Integer=length(sites) ÷ 2, coefficient::Number=1.0)
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("pauli_domain_wall_state requires at least one site"))
  0 <= kink <= nsites || throw(ArgumentError("kink must lie in 0:length(sites)"))
  weights = [j <= kink ? -coefficient : coefficient for j in 1:nsites]
  return operator_local_sum_state(sites, pauli_matrices().Z / sqrt(2), weights)
end
```

In `src/MPSToolkit.jl` add `include("operator_space/states.jl")` immediately before
`include("operator_space/helpers.jl")`, and add `operator_siteinds`, `operator_basis_state`,
`operator_product_state`, `operator_local_sum_state` to the `OperatorSpace` submodule `using` and
`export` lists and to the top-level `export` line that currently lists `pauli_siteinds`.

- [ ] **Step 5: Run tests to verify they pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS. `test_operator_space.jl`, `test_dmt.jl`, `test_pxp.jl` and the doctests all
exercise `pauli_siteinds` / `pauli_basis_state` heavily, so this step is the real regression
check that the wrappers are faithful.

- [ ] **Step 6: Commit**

```bash
git add src/operator_space/states.jl src/operator_space/helpers.jl src/operator_space/daoe.jl src/MPSToolkit.jl test/test_operator_space.jl
git commit -m "feat(operator-space): generic site indices and state builders for any local dimension"
```

---

### Task 3: Generic gate builders

**Files:**
- Create: `src/operator_space/gates.jl`
- Modify: `src/operator_space/helpers.jl` (remove `pauli_gate`, `pauli_gate_from_hamiltonian`, `pauli_lindblad_generator`, `pauli_gate_from_lindbladian`, `_spinhalf_span`, `_pauli_basis_operators`)
- Modify: `src/operator_space/thermal.jl` (`pauli_gate_from_imaginary_time` becomes a wrapper)
- Modify: `src/MPSToolkit.jl`
- Modify: `test/test_operator_space.jl`

**Interfaces:**
- Consumes: `operator_basis_matrices(d)`, `local_dimension(x)`.
- Produces:
  - `_operator_span(matrix_dim::Integer, d::Integer) -> Int`
  - `_operator_basis_operators(nsites::Integer, d::Integer) -> Vector{Matrix{ComplexF64}}`
  - `operator_gate(op::AbstractMatrix; d=2) -> Matrix{ComplexF64}`
  - `operator_gate_from_hamiltonian(h, dt; d=2)`, `operator_gate_from_imaginary_time(h, dbeta; d=2)`
  - `operator_lindblad_generator(h, jumps; d=2)`, `operator_gate_from_lindbladian(h, jumps, dt; d=2)`

- [ ] **Step 1: Write the failing test**

Append to `test/test_operator_space.jl`:

```julia
@testset "generic operator-space gates" begin
  @testset "pauli_gate is the d = 2 case, bit-for-bit" begin
    x = ComplexF64[0 1; 1 0]
    z = ComplexF64[1 0; 0 -1]
    u = exp(-0.3im * kron(x, z))
    @test pauli_gate(u) ≈ operator_gate(u; d = 2) atol = 1e-14
    h = spinhalf_xyz_bond_hamiltonian(; Jx = 1.0, Jy = 0.7, Jz = 0.3)
    @test pauli_gate_from_hamiltonian(h, 0.1) ≈ operator_gate_from_hamiltonian(h, 0.1; d = 2) atol = 1e-14
  end

  @testset "d = 3 gate conjugates correctly" begin
    d = 3
    basis = operator_basis_matrices(d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    sx = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    h = kron(sx, sx) + kron(sz, sz)
    dt = 0.13
    gate = operator_gate_from_hamiltonian(h, dt; d = d)
    @test size(gate) == (d^4, d^4)
    # Compare against the definition on a few random two-site operators.
    u = exp(-im * dt * Matrix{ComplexF64}(h))
    two_site = [kron(a, b) for a in basis for b in basis]
    for column in (1, 7, 40, d^4)
      evolved = u * two_site[column] * u'
      expected = ComplexF64[tr(two_site[row]' * evolved) for row in eachindex(two_site)]
      @test gate[:, column] ≈ expected atol = 1e-10
    end
  end

  @testset "identity Hamiltonian gives the identity superoperator" begin
    for d in (2, 3, 4)
      gate = operator_gate(Matrix{ComplexF64}(I, d, d); d = d)
      @test gate ≈ Matrix{ComplexF64}(I, d^2, d^2) atol = 1e-12
    end
  end

  @testset "sparse Hamiltonians from EDKit-style builders are accepted" begin
    using SparseArrays
    d = 3
    sz = sparse(ComplexF64[1 0 0; 0 0 0; 0 0 -1])
    h = kron(sz, sz)
    @test h isa AbstractSparseMatrix
    dense_gate = operator_gate_from_hamiltonian(Matrix(h), 0.1; d = d)
    @test operator_gate_from_hamiltonian(h, 0.1; d = d) ≈ dense_gate atol = 1e-12
  end

  @testset "span inference rejects incompatible sizes" begin
    @test MPSToolkit._operator_span(9, 3) == 1
    @test MPSToolkit._operator_span(81, 3) == 2
    @test_throws ArgumentError MPSToolkit._operator_span(8, 3)
  end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_operator_space.jl`
Expected: FAIL with `UndefVarError: operator_gate not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/operator_space/gates.jl`:

```julia
"""
    _operator_span(matrix_dim, d)

Infer how many sites of local dimension `d` a dense physical matrix of size
`matrix_dim x matrix_dim` acts on.
"""
function _operator_span(matrix_dim::Integer, d::Integer)
  local_dim = Int(d)
  local_dim >= 2 || throw(ArgumentError("local Hilbert space dimension must be at least 2"))
  size_value = Int(matrix_dim)
  size_value >= local_dim || throw(ArgumentError("operator dimension $(size_value) is smaller than the local dimension $(local_dim)"))
  span = 1
  product = local_dim
  while product < size_value
    product *= local_dim
    span += 1
  end
  product == size_value ||
    throw(ArgumentError("matrix dimension $(size_value) is not a power of the local dimension $(local_dim)"))
  return span
end

"""
    _operator_basis_operators(nsites, d)

Enumerate the normalized multi-site operator-basis matrices for `nsites` sites of local
dimension `d`, ordered so the site-`1` label varies slowest (matching `kron`).
"""
function _operator_basis_operators(nsites::Integer, d::Integer)
  local_basis = operator_basis_matrices(d)
  operators = local_basis
  for _ in 2:Int(nsites)
    operators = [kron(left, right) for left in operators for right in local_basis]
  end
  return operators
end

"""
    operator_gate(op; d=2)

Convert a dense physical map `A` on one or more sites of local dimension `d` into the
corresponding dense superoperator `rho -> A rho A'` in the normalized operator basis.

# Arguments
- `op`: Dense (or sparse) physical matrix of size `d^span x d^span`. Unitary for real-time
  evolution; Hermitian positive for an imaginary-time slice.

# Keyword Arguments
- `d`: Local Hilbert space dimension.

# Returns
- A dense `d^(2 span) x d^(2 span)` matrix `G` with `G[a, b] = tr(P_a' * A * P_b * A')`,
  where `P` are the multi-site operators of [`operator_basis_matrices`](@ref).

# Notes
- Built as `G = W' * kron(conj(A), A) * W`, where `W[:, b] = vec(P_b)`. This costs
  `O(d^(6 span))` rather than the `O(d^(7 span))` of a per-column loop, which matters once
  `d^span` exceeds a few: a three-site spin-1 gate is `19683 x 19683` in operator space, so
  keep `span` small at higher spin.
- Sparse input (for example from `EDKit.spin(...; D = d)`) is densified internally.
"""
function operator_gate(op::AbstractMatrix; d::Integer=2)
  size(op, 1) == size(op, 2) || throw(ArgumentError("operator-space gate requires a square matrix"))
  local_dim = Int(d)
  span = _operator_span(size(op, 1), local_dim)
  dense = Matrix{ComplexF64}(op)
  basis = _operator_basis_operators(span, local_dim)
  nbasis = length(basis)
  physical_dim = size(dense, 1)
  # W[:, b] = vec(P_b) with the same column-major vec used by kron(conj(A), A).
  w = Matrix{ComplexF64}(undef, physical_dim^2, nbasis)
  for b in 1:nbasis
    w[:, b] = vec(basis[b])
  end
  return adjoint(w) * (kron(conj(dense), dense) * w)
end

"""
    operator_gate_from_hamiltonian(h, dt; d=2)

Build the operator-space superoperator induced by `exp(-im * dt * h)` for a dense local
Hamiltonian `h` on sites of local dimension `d`.
"""
function operator_gate_from_hamiltonian(h::AbstractMatrix, dt::Real; d::Integer=2)
  size(h, 1) == size(h, 2) || throw(ArgumentError("operator-space Hamiltonian must be square"))
  return operator_gate(exp(-im * dt * Matrix{ComplexF64}(h)); d=d)
end

"""
    operator_gate_from_imaginary_time(h, dbeta; d=2)

Build the operator-space superoperator for one two-sided imaginary-time slice
`rho -> e^{-(dbeta/2) h} rho e^{-(dbeta/2) h}` for a dense local Hermitian `h`.
"""
function operator_gate_from_imaginary_time(h::AbstractMatrix, dbeta::Real; d::Integer=2)
  size(h, 1) == size(h, 2) || throw(ArgumentError("imaginary-time Hamiltonian must be square"))
  dense = Matrix{ComplexF64}(h)
  norm(dense - dense') <= sqrt(eps(Float64)) * max(norm(dense), one(Float64)) ||
    throw(ArgumentError("imaginary-time Hamiltonian must be Hermitian"))
  return operator_gate(exp(-(dbeta / 2) * dense); d=d)
end

"""
    operator_lindblad_generator(h, jumps; d=2)

Build the dense operator-space Lindbladian generator for a local Hamiltonian `h` and local
jump operators `jumps`, implementing `-im[H, rho] + sum_j (L_j rho L_j' - {L_j' L_j, rho}/2)`.
"""
function operator_lindblad_generator(h::AbstractMatrix, jumps; d::Integer=2)
  size(h, 1) == size(h, 2) || throw(ArgumentError("operator-space Hamiltonian must be square"))
  local_dim = Int(d)
  span = _operator_span(size(h, 1), local_dim)
  basis = _operator_basis_operators(span, local_dim)
  dense_h = Matrix{ComplexF64}(h)
  jump_list = _lindblad_jump_list(jumps)
  for jump in jump_list
    size(jump) == size(dense_h) || throw(ArgumentError("jump operators must match the Hamiltonian dimension"))
  end
  jump_data = [(Matrix{ComplexF64}(jump), Matrix{ComplexF64}(jump)' * Matrix{ComplexF64}(jump)) for jump in jump_list]

  generator = Matrix{ComplexF64}(undef, length(basis), length(basis))
  for column in eachindex(basis)
    evolved = -im * (dense_h * basis[column] - basis[column] * dense_h)
    for (jump, jump_dag_jump) in jump_data
      evolved += jump * basis[column] * jump'
      evolved -= 0.5 * (jump_dag_jump * basis[column] + basis[column] * jump_dag_jump)
    end
    for row in eachindex(basis)
      generator[row, column] = tr(basis[row]' * evolved)
    end
  end
  return generator
end

"""
    operator_gate_from_lindbladian(h, jumps, dt; d=2)

Build the dense operator-space TEBD gate `exp(dt * operator_lindblad_generator(h, jumps; d))`.
"""
function operator_gate_from_lindbladian(h::AbstractMatrix, jumps, dt::Real; d::Integer=2)
  return exp(dt * operator_lindblad_generator(h, jumps; d=d))
end

"""
    pauli_gate(unitary)

Spin-1/2 case of [`operator_gate`](@ref), in the `(I, X, Y, Z)` ordering.
"""
pauli_gate(unitary::AbstractMatrix) = operator_gate(unitary; d=2)

"""
    pauli_gate_from_hamiltonian(h, dt)

Spin-1/2 case of [`operator_gate_from_hamiltonian`](@ref).
"""
pauli_gate_from_hamiltonian(h::AbstractMatrix, dt::Real) = operator_gate_from_hamiltonian(h, dt; d=2)

"""
    pauli_lindblad_generator(h, jumps)

Spin-1/2 case of [`operator_lindblad_generator`](@ref).
"""
pauli_lindblad_generator(h::AbstractMatrix, jumps) = operator_lindblad_generator(h, jumps; d=2)

"""
    pauli_gate_from_lindbladian(h, jumps, dt)

Spin-1/2 case of [`operator_gate_from_lindbladian`](@ref).
"""
pauli_gate_from_lindbladian(h::AbstractMatrix, jumps, dt::Real) = operator_gate_from_lindbladian(h, jumps, dt; d=2)
```

Delete the corresponding originals from `src/operator_space/helpers.jl` (`pauli_gate`,
`pauli_gate_from_hamiltonian`, `pauli_lindblad_generator`, `pauli_gate_from_lindbladian`,
`_spinhalf_span`, `_pauli_basis_operators`), keeping `_lindblad_jump_list` there. In
`src/operator_space/thermal.jl` replace the body of `pauli_gate_from_imaginary_time` with
`return operator_gate_from_imaginary_time(h, dbeta; d=2)`.

`src/operator_space/expectations.jl` calls `_spinhalf_span(size(op, 1))` at line 48; change it to
`_operator_span(size(op, 1), local_dimension(siteind(rho, Int(start))))`. Task 4 rewrites that
file properly; this keeps the package compiling in the meantime.

Add `include("operator_space/gates.jl")` to `src/MPSToolkit.jl` after the `states.jl` include, and
add the five `operator_*` gate names to the `OperatorSpace` submodule and the top-level exports.

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS. Watch `test_operator_space.jl` and the open-systems Lindblad tests in particular.

- [ ] **Step 5: Commit**

```bash
git add src/operator_space/gates.jl src/operator_space/helpers.jl src/operator_space/thermal.jl src/operator_space/expectations.jl src/MPSToolkit.jl test/test_operator_space.jl
git commit -m "feat(operator-space): generic gate and Lindbladian builders for any local dimension"
```

---

### Task 4: Generic trace, expectations, vectorization, and Gibbs preparation

**Files:**
- Modify: `src/operator_space/expectations.jl` (whole file)
- Modify: `src/operator_space/vectorize.jl` (whole file)
- Modify: `src/operator_space/thermal.jl` (`pauli_gibbs_state` -> `operator_gibbs_state` + wrapper)
- Modify: `src/MPSToolkit.jl`
- Modify: `test/test_operator_space.jl`

**Interfaces:**
- Consumes: everything from Tasks 1-3.
- Produces:
  - `operator_trace(rho) -> ComplexF64`
  - `operator_expectation_profile(rho, terms; normalize=true) -> Vector{ComplexF64}`
  - `operator_expectation(rho, op, start; normalize=true) -> ComplexF64`
  - `operator_state_from_mpo(op, sites) -> MPS`, `operator_superoperator_mpo(op, sites) -> MPO`
  - `operator_gibbs_state(sites, terms, weights; nsteps=4, maxdim=64, cutoff=1e-12, initial_state=nothing) -> MPS`
  - `_identity_env(site) -> ITensor` (basis element 1 cap; replaces `_pauli_identity_env`)
  - `_validate_operator_space(psi, start, span) -> Int` returning the uniform local dimension

- [ ] **Step 1: Write the failing test**

Append to `test/test_operator_space.jl`:

```julia
@testset "generic operator-space measurement" begin
  @testset "trace and expectations at d = 3 match dense linear algebra" begin
    d = 3
    nsites = 3
    sites = operator_siteinds(nsites; d = d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    identity_matrix = Matrix{ComplexF64}(I, d, d)
    # rho = I x I x I  +  0.4 * (I x S^z x I) : a genuine near-infinite-temperature operator
    rho = add(
      operator_product_state(sites, [identity_matrix, identity_matrix, identity_matrix]),
      operator_product_state(sites, [identity_matrix, sz, identity_matrix]; coefficient = 0.4);
      maxdim = 4, cutoff = 0.0)

    dense = kron(kron(identity_matrix, identity_matrix), identity_matrix) +
            0.4 * kron(kron(identity_matrix, sz), identity_matrix)
    @test operator_trace(rho) ≈ tr(dense) atol = 1e-10

    for (start, op) in ((2, sz), (1, kron(sz, sz)))
      span = start == 1 ? 2 : 1
      padded = start == 1 ? kron(op, identity_matrix) : kron(kron(identity_matrix, op), identity_matrix)
      @test operator_expectation(rho, op, start; normalize = false) ≈ tr(dense * padded) atol = 1e-10
      @test operator_expectation(rho, op, start) ≈ tr(dense * padded) / tr(dense) atol = 1e-10
    end
  end

  @testset "pauli_* measurement wrappers are unchanged" begin
    sites = pauli_siteinds(4)
    rho = add(pauli_basis_state(sites, fill(1, 4)),
              pauli_domain_wall_state(sites; kink = 2); maxdim = 4, cutoff = 0.0)
    z = pauli_matrices().Z
    terms = [(x, ComplexF64.(z)) for x in 1:4]
    @test pauli_trace(rho) ≈ operator_trace(rho) atol = 1e-12
    @test pauli_expectation_profile(rho, terms; normalize = false) ≈
          operator_expectation_profile(rho, terms; normalize = false) atol = 1e-12
  end

  @testset "operator_state_from_mpo round-trips at d = 3" begin
    d = 3
    nsites = 3
    phys = siteinds("S=1", nsites)
    sites = operator_siteinds(nsites; d = d)
    mpo = MPO(phys, "Id")
    vectorized = operator_state_from_mpo(mpo, sites)
    @test operator_trace(vectorized) ≈ ComplexF64(d^nsites) atol = 1e-10
  end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_operator_space.jl`
Expected: FAIL with `UndefVarError: operator_trace not defined`.

- [ ] **Step 3: Rewrite `expectations.jl` generically**

Replace the spin-1/2 assumptions:

- `_pauli_identity_env(site)` becomes `_identity_env(site)` (unchanged body — it sets basis
  element 1, which is `I / sqrt(d)` for every `d`). Keep a `const _pauli_identity_env = _identity_env`
  alias so `dmt.jl` keeps compiling until Task 6 renames its uses.
- `_validate_pauli_operator_space(psi, start, span)` becomes `_validate_operator_space(psi, start, span)`,
  which computes `d = local_dimension(siteind(psi, start))`, asserts every site in the window has
  the same dimension, and returns `d`.
- `pauli_trace` becomes:

```julia
function operator_trace(rho::MPS)
  d = _validate_operator_space(rho, 1, length(rho))
  return Float64(d)^(length(rho) / 2) * scalar(_right_identity_environment(rho, 1))
end

pauli_trace(rho::MPS) = operator_trace(rho)
```

- `_pauli_window_cap(window_sites, op; local_basis)` becomes `_operator_window_cap(window_sites, op, d)`,
  iterating `Iterators.product(ntuple(_ -> 1:d^2, span)...)` and using
  `operator_basis_matrices(d)` for `local_basis`. Guard `size(op) == (d^span, d^span)`.
- `_validated_pauli_windows(rho, terms)` becomes `_validated_operator_windows(rho, terms, d)` using
  `_operator_span(size(op, 1), d)`.
- In `pauli_expectation_profile`, replace the two `2.0^(...)` factors with `Float64(d)^(...)`:
  `results[k] = normalize ? raw / (Float64(d)^(span / 2) * denominator) : Float64(d)^((nsites - span) / 2) * raw`.
  Rename to `operator_expectation_profile` / `operator_expectation` and add wrappers
  `pauli_expectation_profile(rho, terms; normalize=true) = operator_expectation_profile(rho, terms; normalize=normalize)`
  and the analogous one for `pauli_expectation`.
- Update the docstring note: the `(sqrt d)^N` prefactor overflows `Float64` for
  `N ≳ 2048 / log2(d)`.

- [ ] **Step 4: Rewrite `vectorize.jl` generically**

- `_validate_vectorization_sites(op, sites)` returns `d = local_dimension(first(sites))` after
  checking lengths and uniformity, instead of requiring dimension 4.
- In `operator_state_from_mpo`, replace `_pauli_basis_operators(1)` with
  `operator_basis_matrices(d)`, replace the `dim(phys) == 2` check with `dim(phys) == d`, and
  replace the `for row in 1:2, column in 1:2` loops with `for row in 1:d, column in 1:d`.
- Same three substitutions in `operator_superoperator_mpo`.
- Add wrappers `pauli_state_from_mpo(op, sites) = operator_state_from_mpo(op, sites)` and
  `pauli_superoperator_mpo(op, sites) = operator_superoperator_mpo(op, sites)`.

- [ ] **Step 5: Generalize `thermal.jl`**

Rename `pauli_gibbs_state` to `operator_gibbs_state`. Inside it:
- derive `d = local_dimension(first(sites))`;
- the default `initial_state` becomes `operator_basis_state(sites, fill(1, length(sites)))`
  (unchanged code, but now dimension-generic);
- the gate build becomes `operator_gate_from_imaginary_time(h, weight / (2 * nsteps); d=d)`.

Add `pauli_gibbs_state(args...; kwargs...) = operator_gibbs_state(args...; kwargs...)`.

- [ ] **Step 6: Guard the paths that stay spin-1/2**

DAOE and the PXP constraint helpers are genuinely spin-1/2-specific — `daoe.jl:69` keys off
`local_state == 4` being the `Z` label, and the PXP constraint is a two-level blockade. They are
out of scope, but they must now *say so* rather than silently producing nonsense when handed
`d = 3` sites. At the top of `pauli_daoe_projector`, `pauli_fdaoe_projector`,
`pauli_pxp_constraint_state` and `pauli_pxp_constraint_projector`, add:

```julia
  all(local_dimension(site) == 2 for site in sites) ||
    throw(ArgumentError("this helper is defined for spin-1/2 (local dimension 2) operator space only"))
```

Add a test asserting the error fires:

```julia
  @testset "spin-1/2-only helpers reject higher local dimension" begin
    sites = operator_siteinds(4; d = 3)
    @test_throws ArgumentError pauli_daoe_projector(sites, 2, 0.5)
    @test_throws ArgumentError pauli_pxp_constraint_state(sites)
  end
```

Check the exact positional signatures of those four functions before writing the test; adjust the
call arguments to match.

- [ ] **Step 7: Update the module and run tests**

Add `operator_trace`, `operator_expectation`, `operator_expectation_profile`,
`operator_state_from_mpo`, `operator_superoperator_mpo`, `operator_gibbs_state` to the
`OperatorSpace` submodule and top-level exports in `src/MPSToolkit.jl`.

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS, including `test_transport_reference.jl` and `test_oracles.jl`, which pin
numerical values through these functions.

- [ ] **Step 8: Commit**

```bash
git add src/operator_space/expectations.jl src/operator_space/vectorize.jl src/operator_space/thermal.jl src/operator_space/daoe.jl src/operator_space/pxp.jl src/MPSToolkit.jl test/test_operator_space.jl
git commit -m "feat(operator-space): generic trace, expectations, vectorization and Gibbs preparation"
```

---

### Task 5: DMT truncation primitives (pure matrix level)

Nothing is wired into the MPS path in this task; it builds and tests the linear algebra alone, so
a failure here is unambiguous.

**Files:**
- Create: `src/operator_space/dmt/lowrank.jl`
- Modify: `src/MPSToolkit.jl` (include)
- Create: `test/test_dmt_kernel.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `_protected_basis(protected::AbstractMatrix) -> Matrix{ComplexF64}` — thin orthonormal basis (QR), `chi x k`
  - `_bond_mul(bond_matrix, x)`, `_bond_adj(bond_matrix, x)` — bond-matrix products with a `Diagonal` fast path
  - `_dmt_connector(bond_matrix::AbstractMatrix, q0, r0) -> (a, b, has_connector)` where the connector is `a * transpose(b)`
  - `_dmt_complement_ops(bond_matrix, a, b, QL, QR) -> NamedTuple` with fields `mul`, `adj`, `QLtB`, `BQR`, `BQRc`
  - `_truncated_svd(mul, adj, chi, rank; mode=:dense, oversample=10, power=2) -> (U, S, V)`
  - `_dmt_refactor(factors_left, factors_right, maxdim, cutoff) -> (U, S, V)` — QR-based low-rank SVD of `F * G'`

The bond matrix is taken as a **matrix** from the start, not a vector of Schmidt values, so that
Task 9 can swap the SVD bond factorization for a QR (which yields a triangular bond matrix)
without changing any of these signatures.

- [ ] **Step 1: Write the failing test**

Create `test/test_dmt_kernel.jl`:

```julia
using LinearAlgebra
using MPSToolkit
using Random
using Test

@testset "DMT truncation primitives" begin
  Random.seed!(20260808)

  @testset "_protected_basis is orthonormal and spans the input" begin
    for (chi, k) in ((40, 4), (40, 9), (12, 16))
      protected = randn(ComplexF64, chi, k)
      basis = MPSToolkit._protected_basis(protected)
      @test size(basis, 1) == chi
      @test size(basis, 2) == min(chi, k)
      @test basis' * basis ≈ Matrix{ComplexF64}(I, size(basis, 2), size(basis, 2)) atol = 1e-12
      if k <= chi
        @test norm(protected - basis * (basis' * protected)) < 1e-10
      end
    end
  end

  @testset "the connector is rank one and annihilated by B" begin
    for bond_matrix in (Diagonal(ComplexF64.(sort(rand(30); rev = true))),
                        UpperTriangular(randn(ComplexF64, 30, 30)))
      chi = size(bond_matrix, 1)
      q0 = normalize(randn(ComplexF64, chi))
      r0 = normalize(randn(ComplexF64, chi))
      a, b, has_connector = MPSToolkit._dmt_connector(bond_matrix, q0, r0)
      @test has_connector
      connector = a * transpose(b)
      @test rank(connector; atol = 1e-10) == 1
      bmat = Matrix(bond_matrix) - connector
      # B annihilates the identity directions: this is why the connector needs no extra budget.
      @test norm(bmat * r0) < 1e-10
      @test norm(q0' * bmat) < 1e-10
    end
  end

  @testset "a negligible trace overlap disables the connector rather than blowing up" begin
    chi = 20
    bond_matrix = Diagonal(ComplexF64.(sort(rand(chi); rev = true)))
    q0 = zeros(ComplexF64, chi); q0[1] = 1
    r0 = zeros(ComplexF64, chi); r0[2] = 1     # orthogonal directions => q0' S r0 == 0
    a, b, has_connector = MPSToolkit._dmt_connector(bond_matrix, q0, r0)
    @test !has_connector
    @test all(iszero, a)
  end

  @testset "the complement operator is doubly orthogonal to the protected subspaces" begin
    chi, k = 40, 9
    bond_matrix = Diagonal(ComplexF64.(sort(rand(chi); rev = true)))
    QL = MPSToolkit._protected_basis(randn(ComplexF64, chi, k))
    QR = MPSToolkit._protected_basis(randn(ComplexF64, chi, k))
    q0 = QL[:, 1]
    r0 = QR[:, 1]
    a, b, _ = MPSToolkit._dmt_connector(bond_matrix, q0, r0)
    ops = MPSToolkit._dmt_complement_ops(bond_matrix, a, b, QL, QR)
    dense = ops.mul(Matrix{ComplexF64}(I, chi, chi))
    @test norm(QL' * dense) < 1e-10        # rows orthogonal to the left protected space
    @test norm(dense * QR) < 1e-10         # columns orthogonal to the right protected space
    # adj really is the adjoint
    x = randn(ComplexF64, chi, 3)
    y = randn(ComplexF64, chi, 3)
    @test tr(y' * ops.mul(x)) ≈ tr((ops.adj(y))' * x) atol = 1e-10
  end

  @testset "randomized truncated SVD matches the dense one" begin
    chi, rank_target = 200, 20
    base = randn(ComplexF64, chi, chi)
    mul(x) = base * x
    adj(x) = base' * x
    u_dense, s_dense, v_dense = MPSToolkit._truncated_svd(mul, adj, chi, rank_target; mode = :dense)
    u_rand, s_rand, v_rand = MPSToolkit._truncated_svd(mul, adj, chi, rank_target; mode = :random)
    @test length(s_dense) == rank_target
    @test length(s_rand) == rank_target
    @test s_rand ≈ s_dense rtol = 0.05
    dense_error = norm(base - u_dense * Diagonal(s_dense) * v_dense')
    rand_error = norm(base - u_rand * Diagonal(s_rand) * v_rand')
    @test rand_error <= 1.05 * dense_error
  end

  @testset "_dmt_refactor reproduces the dense SVD of a low-rank product" begin
    chi, r = 300, 25
    f = randn(ComplexF64, chi, r)
    g = randn(ComplexF64, chi, r)
    u, s, v = MPSToolkit._dmt_refactor(f, g, 1000, 1e-14)
    @test size(u, 2) == length(s) == size(v, 2)
    @test length(s) <= r
    @test u * Diagonal(s) * v' ≈ f * g' atol = 1e-9
    @test u' * u ≈ Matrix{ComplexF64}(I, length(s), length(s)) atol = 1e-10
    # maxdim clips
    u2, s2, v2 = MPSToolkit._dmt_refactor(f, g, 10, 0.0)
    @test length(s2) == 10
  end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_dmt_kernel.jl`
Expected: FAIL with `UndefVarError: _protected_basis not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/operator_space/dmt/lowrank.jl`:

```julia
"""
    _protected_basis(protected)

Return a thin orthonormal basis whose columns span the columns of `protected`.

# Notes
- Uses a QR factorization rather than an SVD: only the span matters, the input is
  `chi x d^(2n)` with `d^(2n)` tiny compared with `chi`, and a thin QR is orders of magnitude
  cheaper than `svd(...; full=true)` (measured 0.0013 s versus 0.505 s at `chi = 1800`).
- No rank detection is performed. A rank-deficient `protected` yields a few arbitrary extra
  orthonormal columns, which over-protects rather than under-protects, and is therefore safe.
"""
function _protected_basis(protected::AbstractMatrix)
  size(protected, 2) == 0 && return zeros(ComplexF64, size(protected, 1), 0)
  ncols = min(size(protected, 1), size(protected, 2))
  factorization = qr(Matrix{ComplexF64}(protected))
  return Matrix(factorization.Q)[:, 1:ncols]
end

"""
    _bond_mul(bond_matrix, x)
    _bond_adj(bond_matrix, x)

Apply the bond matrix (or its adjoint) to `x`, with an `O(chi)` fast path when the bond matrix
is `Diagonal`. The bond matrix is diagonal when the bond tensor is factorized by SVD and
triangular when it is factorized by QR; the DMT rule is identical either way.
"""
_bond_mul(bond_matrix::Diagonal, x) = diag(bond_matrix) .* x
_bond_mul(bond_matrix::AbstractMatrix, x) = bond_matrix * x
_bond_adj(bond_matrix::Diagonal, x) = conj(diag(bond_matrix)) .* x
_bond_adj(bond_matrix::AbstractMatrix, x) = bond_matrix' * x

"""
    _dmt_connector(bond_matrix, q0, r0)

Return `(a, b, has_connector)` describing the rank-one trace connector
`C = a * transpose(b) = (S r0)(q0' S) / (q0' S r0)` for the bond matrix `S`.

# Notes
- When the identity overlap `q0' S r0` is negligible relative to the operator scale — a
  traceless operator has no trace component to anchor on — the connector is disabled and `a`
  is returned as zeros, so the caller still truncates and still enforces `maxdim`.
- `B = S - C` annihilates the identity directions (`B r0 = 0`, `q0' B = 0`), so the connector
  costs no extra bond dimension: its range already lies inside the protected block.
"""
function _dmt_connector(bond_matrix::AbstractMatrix, q0::AbstractVector, r0::AbstractVector)
  scaled_right = _bond_mul(bond_matrix, r0)
  denominator = dot(q0, scaled_right)
  has_connector = abs(denominator) > sqrt(eps(Float64)) * norm(bond_matrix)
  a = has_connector ? ComplexF64.(scaled_right ./ denominator) : zeros(ComplexF64, length(r0))
  b = ComplexF64.(_bond_adj(bond_matrix, q0))    # (q0' S)^T = conj(S' q0)... see note
  return a, b, has_connector
end
```

The `b` line deserves care: the connector is `C[i, j] = (S r0)_i * (q0' S)_j / den`, so
`b_j = (q0' S)_j = sum_i conj(q0_i) S[i, j]`, i.e. `b = transpose(S) * conj(q0) = conj(S' * q0)`.
Write it as `b = ComplexF64.(conj(_bond_adj(bond_matrix, q0)))` and pin it with this test, which
must be added to the connector testset above:

```julia
      # b is the row vector q0' S, not S' q0
      @test transpose(b) ≈ q0' * Matrix(bond_matrix) atol = 1e-10
```

Continuing `lowrank.jl`:

```julia
"""
    _dmt_complement_ops(bond_matrix, a, b, QL, QR)

Return matrix-free operators for `D = P_L^perp B P_R^perp` with `B = S - a b^T` and
`P^perp = I - Q Q'`.

# Returns
- A named tuple with `mul(X) = D * X`, `adj(X) = D' * X`, and the reusable intermediates
  `QLtB = QL' * B`, `BQR = B * QR`, `BQRc = B QR - QL (QL' B QR)`.

# Notes
- Nothing here is `chi x chi`. With a diagonal bond matrix every product costs
  `O(chi * ncols * d^(2n))`, which is what makes the randomized truncation in
  [`_truncated_svd`](@ref) cheap; with a triangular one it is `O(chi^2 * ncols)`.
"""
function _dmt_complement_ops(bond_matrix::AbstractMatrix, a::AbstractVector, b::AbstractVector,
                             QL::AbstractMatrix, QR::AbstractMatrix)
  bmul(x) = _bond_mul(bond_matrix, x) .- a * (transpose(b) * x)
  badj(x) = _bond_adj(bond_matrix, x) .- conj(b) * (a' * x)
  qltb = badj(QL)'
  bqr = bmul(QR)
  bqrc = bqr - QL * (qltb * QR)
  mul(x) = bmul(x) - QL * (qltb * x) - bqrc * (QR' * x)
  function adj(x)
    w = badj(x) - qltb' * (QL' * x)
    return w - QR * (QR' * w)
  end
  return (mul=mul, adj=adj, QLtB=qltb, BQR=bqr, BQRc=bqrc)
end

"""
    _truncated_svd(mul, adj, chi, rank; mode=:dense, oversample=10, power=2)

Return the leading `rank` singular triplet `(U, S, V)` of the `chi x chi` operator defined by
the matrix-free products `mul` and `adj`.

# Keyword Arguments
- `mode`: `:dense` materializes the operator and calls LAPACK; `:random` uses a block
  randomized range finder with `power` power iterations and `oversample` extra probe vectors.

# Notes
- The randomized branch cannot break DMT's preservation guarantee: its range lies inside the
  range of the operator, which is already doubly orthogonal to the protected subspaces, so the
  protected border and connector are still reinstated exactly. It only affects how optimal the
  discarded weight is (measured overlap 0.998-0.999 with the dense result).
"""
function _truncated_svd(mul, adj, chi::Integer, rank::Integer; mode::Symbol=:dense,
                        oversample::Integer=10, power::Integer=2)
  keep = max(min(Int(rank), Int(chi)), 0)
  keep == 0 && return (zeros(ComplexF64, chi, 0), Float64[], zeros(ComplexF64, chi, 0))
  if mode === :dense
    factorization = svd(mul(Matrix{ComplexF64}(I, chi, chi)))
    return (factorization.U[:, 1:keep], factorization.S[1:keep], factorization.V[:, 1:keep])
  end
  mode === :random || throw(ArgumentError("DMT truncated SVD mode must be :dense or :random"))
  probes = min(Int(chi), keep + Int(oversample))
  sketch = mul(randn(ComplexF64, Int(chi), probes))
  for _ in 1:Int(power)
    sketch = mul(adj(Matrix(qr(sketch).Q)))
  end
  range_basis = Matrix(qr(sketch).Q)
  projected = adj(range_basis)'                   # probes x chi
  factorization = svd(projected)
  keep = min(keep, length(factorization.S))
  return (range_basis * factorization.U[:, 1:keep], factorization.S[1:keep], factorization.V[:, 1:keep])
end

"""
    _dmt_refactor(left, right, maxdim, cutoff)

Return the SVD `(U, S, V)` of the low-rank product `left * right'` without forming it.

# Notes
- Costs `O(chi r^2 + r^3)` for `chi x r` inputs, versus `O(chi^3)` for a dense SVD (measured
  0.60 s versus 14.2 s at `chi = 3600, r = 400`). The DMT bond matrix is always available in
  this factored form, so the dense matrix is never needed.
"""
function _dmt_refactor(left::AbstractMatrix, right::AbstractMatrix, maxdim::Integer, cutoff::Real)
  size(left, 1) == size(right, 1) || throw(ArgumentError("DMT refactorization needs matching row counts"))
  size(left, 2) == size(right, 2) || throw(ArgumentError("DMT refactorization needs matching column counts"))
  left_qr = qr(Matrix{ComplexF64}(left))
  right_qr = qr(Matrix{ComplexF64}(right))
  core = svd(Matrix(left_qr.R) * Matrix(right_qr.R)')
  keep = length(core.S)
  if keep > 0 && cutoff > 0
    keep = max(count(>(cutoff * core.S[1]), core.S), 1)
  end
  keep = min(keep, Int(maxdim))
  return (Matrix(left_qr.Q) * core.U[:, 1:keep], core.S[1:keep], Matrix(right_qr.Q) * core.V[:, 1:keep])
end
```

Add `include("operator_space/dmt/lowrank.jl")` to `src/MPSToolkit.jl` immediately before the
existing `include("operator_space/dmt.jl")`, and `include("test_dmt_kernel.jl")` to
`test/runtests.jl`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=. test/test_dmt_kernel.jl`
Expected: PASS.

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS (nothing else uses these yet).

- [ ] **Step 5: Commit**

```bash
git add src/operator_space/dmt/lowrank.jl src/MPSToolkit.jl test/test_dmt_kernel.jl test/runtests.jl
git commit -m "feat(dmt): matrix-level primitives for the faithful truncation kernel"
```

---

### Task 6: Faithful bond truncation wired in at d = 2

This is the load-bearing task. It replaces the kernel and the option plumbing in one atomic
change, because the old `connector_buffer` semantics and the new budget semantics cannot both be
true at once.

**Files:**
- Create: `src/operator_space/dmt/bond.jl`
- Modify: `src/operator_space/dmt.jl` (delete `_mat_trunc!`, `_complete_orthonormal_basis`, `_dmt_bond_truncate!`; keep the env cache, window, step and driver code)
- Modify: `src/evolution/types.jl` (`DMTGateEvolution`: `connector_buffer` -> `preserve_diameter`)
- Modify: `src/operator_space/constrained.jl:74`, `src/scarfinder/algorithm.jl:102`
- Modify: `test/test_dmt.jl`
- Create: `test/dmt_test_helpers.jl`, `test/test_dmt_preservation.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: `_protected_basis`, `_dmt_connector`, `_dmt_complement_ops`, `_truncated_svd`, `_dmt_refactor` (Task 5); `local_dimension` (Task 1); `_identity_env`, `_left_identity_environment`, `_right_identity_environment` (existing, generalized in Task 4).
- Produces:
  - `_dmt_protected_sites(bond, nsites, radius) -> (left_count, right_count)`
  - `_dmt_bond_truncate!(psi, bond; maxdim, cutoff, direction=:R, preserve_diameter=3, truncation=:dense, orthogonalize=true, cache=nothing) -> MPS`
  - `_dmt_complement_budget(maxdim, d, radius) -> Int`
  - `_validate_dmt_budget(psi, maxdim, preserve_diameter) -> Nothing`
  - `DMTOptions(; maxdim=30, cutoff=1e-12, gate_maxdim=..., preserve_diameter=3, truncation=:dense)`
  - `DMTGateEvolution(...; preserve_diameter=3, truncation=:dense, ...)` with a `preserve_diameter::Int` and `truncation::Symbol` field replacing `connector_buffer::Int`

- [ ] **Step 1: Write the failing preservation test**

First create `test/dmt_test_helpers.jl`, shared by all three DMT test files:

```julia
using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit

"""
All dense probes of diameter <= `diameter` on an `nsites` chain of local dimension `d`,
as `(start, matrix)` pairs suitable for `operator_expectation_profile`.
"""
function diameter_probes(nsites::Int, d::Int, diameter::Int)
  basis = operator_basis_matrices(d)
  # locals[1] is the identity; the rest are a few non-identity basis elements. Every
  # combination over the window is probed, so a middle site is not silently left as identity
  # -- the guarantee covers ALL operators of this diameter, so the test must too.
  locals = [basis[k] * sqrt(d) for k in 1:min(length(basis), 4)]
  probes = Tuple{Int,Matrix{ComplexF64}}[]
  for width in 1:diameter, start in 1:(nsites - width + 1)
    for labels in Iterators.product(ntuple(_ -> eachindex(locals), width)...)
      all(isequal(1), labels) && continue          # skip the pure-identity probe (that is the trace)
      op = locals[labels[1]]
      for offset in 2:width
        op = kron(op, locals[labels[offset]])
      end
      push!(probes, (start, ComplexF64.(op)))
    end
  end
  return probes
end

function preservation_error(before, after)
  scale = maximum(abs, before)
  scale == 0 && return maximum(abs, after)
  return maximum(abs, after .- before) / scale
end
```

Then create `test/test_dmt_preservation.jl`. This is the test the suite has never had.

```julia
using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Test

include("dmt_test_helpers.jl")

@testset "DMT preserves local observables exactly" begin
  Random.seed!(20260808)

  @testset "diameter <= 3 preserved to machine precision (d = 2)" begin
    nsites, chi, maxdim = 7, 40, 24
    sites = operator_siteinds(nsites; d = 2)
    probes = diameter_probes(nsites, 2, 3)
    for (label, noise) in (("Hermitian (real coefficients)", random_mps(sites; linkdims = chi)),
                           ("non-Hermitian (complex)", random_mps(ComplexF64, sites; linkdims = chi)))
      rho = add(operator_basis_state(sites, fill(1, nsites)), 0.3 * noise;
                maxdim = chi + 1, cutoff = 0.0)
      before = operator_expectation_profile(rho, probes; normalize = false)
      trace_before = operator_trace(rho)
      MPSToolkit._dmt_bond_truncate!(rho, 4; maxdim = maxdim, cutoff = 0.0)
      after = operator_expectation_profile(rho, probes; normalize = false)
      @test dim(linkind(rho, 4)) <= maxdim
      @test preservation_error(before, after) < 1e-11
      @test abs(operator_trace(rho) - trace_before) <= 1e-11 * abs(trace_before)
    end
  end

  @testset "the guarantee survives a full sweep" begin
    nsites, maxdim, dt = 10, 24, 0.1
    sites = operator_siteinds(nsites; d = 2)
    gate = operator_gate_from_hamiltonian(
      spinhalf_xyz_bond_hamiltonian(; Jx = 1.0, Jy = 1.0, Jz = 1.0), dt; d = 2)
    rho = add(operator_basis_state(sites, fill(1, nsites)),
              0.25 * pauli_domain_wall_state(sites; kink = nsites ÷ 2); maxdim = 8, cutoff = 0.0)
    z = ComplexF64[1 0; 0 -1]
    charge(state) = sum(real.(operator_expectation_profile(
      state, [(x, z) for x in 1:nsites]; normalize = false)))
    initial_charge = charge(rho)
    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule = schedule, reverse_schedule = reverse(schedule),
                           nstep = 1, maxdim = maxdim, cutoff = 1e-14, normalize = false)
    for _ in 1:15
      dmt_evolve!(rho, evo)
      @test abs(charge(rho) - initial_charge) < 1e-10
      @test maximum(dim(linkind(rho, b)) for b in 1:(nsites - 1)) <= maxdim
    end
  end

  @testset "purity monotone: DMT keeps Z >= 1 where plain SVD truncation does not" begin
    nsites, chi, maxdim = 8, 40, 20
    sites = operator_siteinds(nsites; d = 2)
    rho = add(operator_basis_state(sites, fill(1, nsites)),
              0.3 * random_mps(sites; linkdims = chi); maxdim = chi + 1, cutoff = 0.0)
    # Z = tr(rho) / sqrt(tr(rho^2)); in an orthonormal operator basis tr(rho^2) = <psi|psi>.
    purity_ratio(state) = real(operator_trace(state)) / sqrt(real(inner(state, state)))
    dmt_state = copy(rho)
    MPSToolkit._dmt_bond_truncate!(dmt_state, 4; maxdim = maxdim, cutoff = 0.0)
    svd_state = copy(rho)
    orthogonalize!(svd_state, 4)
    u, s, v = svd(svd_state[4], (linkind(svd_state, 3), siteind(svd_state, 4));
                  maxdim = maxdim, cutoff = 0.0)
    svd_state[4] = u
    svd_state[5] = s * v * svd_state[5]
    @test purity_ratio(dmt_state) >= purity_ratio(svd_state)
  end

  @testset "budget validation names the local dimension" begin
    sites = operator_siteinds(6; d = 2)
    rho = random_mps(ComplexF64, sites; linkdims = 20)
    # 2 d^2 + 1 = 9 for d = 2
    @test_throws ArgumentError MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim = 8, cutoff = 0.0)
    @test MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim = 9, cutoff = 0.0) === rho
    @test_throws ArgumentError DMTOptions(maxdim = 30, preserve_diameter = 4)   # must be odd
  end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_dmt_preservation.jl`
Expected: FAIL — the first testset reports `preservation_error` around `0.2-0.5` rather than
`< 1e-11`, since the current kernel clips the protected block. Record the number you see; it is
the defect being fixed.

- [ ] **Step 3: Write the faithful bond truncation**

Create `src/operator_space/dmt/bond.jl`:

```julia
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
    _dmt_bond_truncate!(psi, bond; maxdim, cutoff, direction=:R, preserve_diameter=3,
                        truncation=:dense, orthogonalize=true, cache=nothing)

Perform one DMT-preserving bond truncation.

The bond matrix is expressed in a basis whose leading directions span the local-operator
subspace on the sites adjacent to the cut; those `d^(2 radius)` rows and columns, together with
the rank-one trace connector, are carried through **exactly**, and only the doubly-orthogonal
complement is truncated — to `chi' = maxdim - 2 d^(2 radius)` directions, so the reinstated
protected block always fits inside `maxdim` and never has to be clipped.

# Keyword Arguments
- `maxdim`: Total post-truncation bond dimension, inclusive of the protected block.
- `cutoff`: Relative cutoff applied in the final refactorization.
- `direction`: `:R` leaves the orthogonality center at `bond + 1`, `:L` at `bond`.
- `preserve_diameter`: Odd; every observable of diameter at most this value is preserved
  exactly. `radius = (preserve_diameter - 1) / 2` sites are protected per side.
- `truncation`: `:dense` or `:random` complement truncation (see [`_truncated_svd`](@ref)).
- `orthogonalize`: Re-gauge so the orthogonality center is at `bond` first. Set `false` only
  when the center is already known to be there.

# Returns
- The mutated `psi`.
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
  right_site = siteind(psi, bond + 1)
  d = local_dimension(left_site)

  left_env = isnothing(cache) ? _left_identity_environment(psi, bond - left_count) :
             _left_env_at!(cache, psi, bond - left_count)
  right_env = isnothing(cache) ? _right_identity_environment(psi, bond + 1 + right_count) :
              _right_env_at!(cache, psi, bond + 1 + right_count)

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

  bond_values = real.(diag(matrix(s, left_link, right_link)))
  chi = length(bond_values)
  ql = _protected_basis(protected_left)
  qr_basis = _protected_basis(protected_right)
  q0 = normalize(protected_left[:, 1])
  r0 = normalize(protected_right[:, 1])
  a, b, _ = _dmt_connector(bond_values, q0, r0)
  ops = _dmt_complement_ops(bond_values, a, b, ql, qr_basis)

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
```

In both branches `psi[bond] * psi[bond + 1]` equals `u * (new_u * Diagonal(new_s) * new_v') * right_block`;
they differ only in which tensor carries the singular values. Assert that in the test suite by
checking `abs(inner(a, b)) ≈ 1` between a `:R` and an `:L` truncation of the same input.

- [ ] **Step 4: Delete the superseded code and re-plumb the options**

From `src/operator_space/dmt.jl` delete `_mat_trunc!` (lines 238-283), `_complete_orthonormal_basis`
(285-313), the old `_dmt_bond_truncate!` (315-410), and `_validate_pauli_operator_space` (211-216).
In the remaining code:

- `_validate_dmt_step`: drop the `connector_buffer` argument and the dimension-4 check; call
  `_validate_operator_space(psi, start, span)` instead, and keep the span / boundary checks.
- `_dmt_window_truncate!` and `dmt_step!`: replace the `connector_buffer` keyword with
  `preserve_diameter::Integer=3, truncation::Symbol=:dense` and forward both.
- `DMTOptions`: replace the `connector_buffer::Int` field with `preserve_diameter::Int` and
  `truncation::Symbol`; validate `isodd(preserve_diameter) && preserve_diameter >= 1` and
  `truncation in (:dense, :random)`. Delete the `connector_buffer <= maxdim` check.
- `dmt_evolve!`: forward `evo.preserve_diameter` and `evo.truncation`.
- Add `include("operator_space/dmt/bond.jl")` to `src/MPSToolkit.jl` before `dmt.jl`.

In `src/evolution/types.jl`, replace the `connector_buffer::Int` field of `DMTGateEvolution` with
`preserve_diameter::Int` and `truncation::Symbol`, update the constructor keyword defaults to
`preserve_diameter=3, truncation=:dense`, and update the docstrings. Add an explicit migration
error at the top of the constructor:

```julia
function DMTGateEvolution(gate, dt; schedule, reverse_schedule=reverse(schedule), nstep=1,
                          maxdim=30, cutoff=1e-12, gate_maxdim=max(Int(maxdim) * 16, 64),
                          preserve_diameter=3, truncation=:dense, normalize=true,
                          connector_buffer=nothing)
  isnothing(connector_buffer) || throw(ArgumentError(
    "connector_buffer was removed: DMT now protects the d^(preserve_diameter - 1) local-operator " *
    "subspace on each side structurally. Use preserve_diameter (odd, default 3) instead, and note " *
    "that maxdim is now the total bond dimension including the protected block " *
    "(floor 2 d^(preserve_diameter - 1) + 1)."))
  ...
```

Add the same `connector_buffer=nothing` guard to `DMTOptions` and `dmt_step!`.

Update `src/operator_space/constrained.jl:74` and `src/scarfinder/algorithm.jl:102` to forward
`preserve_diameter=evo.preserve_diameter` and `truncation=evo.truncation` instead of
`connector_buffer=evo.connector_buffer`.

- [ ] **Step 5: Update the existing DMT tests to the new contract**

In `test/test_dmt.jl`:
- Replace every `connector_buffer = N` keyword with `preserve_diameter = 3` (the structural
  protection replaces the tuned buffer), and raise any `maxdim` below `9` to at least `9`.
- Delete the `"DMT validates connector buffer budget"` testset and the
  `"reduced-matrix truncation enforces maxdim for traceless operators"` testset (they test
  `_mat_trunc!`, which no longer exists). Replace them with:

```julia
  @testset "connector_buffer raises a migration error" begin
    @test_throws ArgumentError DMTOptions(maxdim = 30, connector_buffer = 8)
    @test_throws ArgumentError DMTGateEvolution(_identity_gate(2), 0.1; schedule = [1],
                                                maxdim = 30, connector_buffer = 8)
  end

  @testset "maxdim is the total budget" begin
    sites = operator_siteinds(6; d = 2)
    psi = random_mps(ComplexF64, sites; linkdims = 40)
    normalize!(psi)
    MPSToolkit._dmt_bond_truncate!(psi, 3; maxdim = 20, cutoff = 1e-14)
    @test dim(linkind(psi, 3)) <= 20
    # 2 d^2 = 8 of those 20 are the protected block
    @test dim(linkind(psi, 3)) >= 8
  end
```

- Replace the `"DMT preserves total S^z ..."` testset tolerance `1e-2` with `1e-10`: the faithful
  kernel conserves exactly, and a loose tolerance here is what hid the defect.
- Delete the `"_complete_orthonormal_basis edge cases"`, `"complex DMT projection uses adjoint
  orthonormal bases"`, and `"orthonormal basis completion ..."` testsets; `_protected_basis` is
  covered by `test_dmt_kernel.jl`.
- Keep the environment-cache, schedule, and reverse-sweep testsets unchanged apart from the
  keyword rename.
- Change the `"DMT validates Pauli dimension for single-site gates"` testset: dimension-3 sites
  are no longer invalid *as operator sites* (3 is not a perfect square, so it still throws) —
  keep it, but rename it to `"DMT rejects sites whose dimension is not a perfect square"`.

Add `include("test_dmt_preservation.jl")` to `test/runtests.jl`.

- [ ] **Step 6: Run tests to verify they pass**

Run: `julia --project=. test/test_dmt_preservation.jl`
Expected: PASS — `preservation_error` now below `1e-11` where Step 2 measured `0.2-0.5`.

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS. `test_transport_reference.jl` pins a committed chi=48 correlator reference
produced by the old kernel; if it fails, regenerate the reference with
`julia --project=. examples/operator_space/pxp_energy_correlator_reference.jl` and commit the new
CSV together with a note in the commit message that the values moved because the kernel was fixed.

- [ ] **Step 7: Commit**

```bash
git add src/operator_space/dmt/bond.jl src/operator_space/dmt.jl src/evolution/types.jl \
        src/operator_space/constrained.jl src/scarfinder/algorithm.jl src/MPSToolkit.jl \
        test/test_dmt.jl test/test_dmt_preservation.jl test/dmt_test_helpers.jl test/runtests.jl
git commit -m "fix(dmt): faithful truncation that preserves diameter-3 observables exactly

The protected block was reinstated and then clipped away by the final repair
SVD, so no setting of connector_buffer delivered DMT's guarantee: diameter-<=3
observables came out 26% wrong from a single bond truncation. Size the
complement budget as chi' = maxdim - 2 d^2 instead, protect the local-operator
subspace structurally rather than by a tuned direction count, and fix the
left-basis conjugation to match the paper's M = Q_L^T s Q_R pairing."
```

---

### Task 7: Generic local dimension and wider preservation

**Files:**
- Modify: `test/test_dmt_preservation.jl`
- Create: `test/test_dmt_higher_spin.jl`
- Modify: `test/runtests.jl`
- Modify: `src/operator_space/dmt/bond.jl` only if a test exposes a `d`-specific bug

The kernel written in Task 6 is already generic in `d` and in `radius`. This task proves it, and
fixes whatever the proof breaks.

**Interfaces:**
- Consumes: everything from Tasks 1-6.
- Produces: no new API.

- [ ] **Step 1: Write the failing higher-spin test**

Create `test/test_dmt_higher_spin.jl`:

```julia
using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Random
using Test

include("dmt_test_helpers.jl")   # diameter_probes / preservation_error

@testset "DMT at higher local dimension" begin
  Random.seed!(20260808)

  @testset "diameter <= 3 preserved at d = 3 and d = 4" begin
    for d in (3, 4)
      nsites, chi = 6, 60
      maxdim = 2 * d^2 + 12
      sites = operator_siteinds(nsites; d = d)
      probes = diameter_probes(nsites, d, 3)
      for noise in (random_mps(sites; linkdims = chi),
                    random_mps(ComplexF64, sites; linkdims = chi))
        rho = add(operator_basis_state(sites, fill(1, nsites)), 0.3 * noise;
                  maxdim = chi + 1, cutoff = 0.0)
        before = operator_expectation_profile(rho, probes; normalize = false)
        MPSToolkit._dmt_bond_truncate!(rho, 3; maxdim = maxdim, cutoff = 0.0)
        @test dim(linkind(rho, 3)) <= maxdim
        @test preservation_error(before, operator_expectation_profile(rho, probes; normalize = false)) < 1e-11
      end
    end
  end

  @testset "preserve_diameter = 5 preserves diameter 5 and stops at diameter 6" begin
    d, nsites, chi = 2, 9, 60
    maxdim = 2 * d^4 + 12      # 2 * 16 + 12 = 44
    sites = operator_siteinds(nsites; d = d)
    kept = diameter_probes(nsites, d, 5)
    rho = add(operator_basis_state(sites, fill(1, nsites)),
              0.3 * random_mps(ComplexF64, sites; linkdims = chi); maxdim = chi + 1, cutoff = 0.0)
    reference = copy(rho)
    before = operator_expectation_profile(rho, kept; normalize = false)
    MPSToolkit._dmt_bond_truncate!(rho, 5; maxdim = maxdim, cutoff = 0.0)
    @test preservation_error(before, operator_expectation_profile(rho, kept; normalize = false)) < 1e-11

    # At preserve_diameter = 3 the same diameter-5 probes are NOT preserved: this pins that the
    # parameter does something rather than being cosmetic.
    narrow = copy(reference)
    MPSToolkit._dmt_bond_truncate!(narrow, 5; maxdim = maxdim, preserve_diameter = 3, cutoff = 0.0)
    wide_only = [p for p in kept if size(p[2], 1) == d^5]
    @test preservation_error(
      operator_expectation_profile(reference, wide_only; normalize = false),
      operator_expectation_profile(narrow, wide_only; normalize = false)) > 1e-8
  end

  @testset "spin-1 Heisenberg melt conserves total S^z end to end" begin
    d, nsites, dt = 3, 8, 0.1
    maxdim = 2 * d^2 + 22      # 40
    sites = operator_siteinds(nsites; d = d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    sx = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    sy = ComplexF64[0 -im 0; im 0 -im; 0 im 0] / sqrt(2)
    h = kron(sx, sx) + kron(sy, sy) + kron(sz, sz)
    gate = operator_gate_from_hamiltonian(h, dt; d = d)
    identity_matrix = Matrix{ComplexF64}(I, d, d)
    rho = add(operator_product_state(sites, fill(identity_matrix, nsites)),
              operator_local_sum_state(sites, sz, [j <= nsites ÷ 2 ? -0.25 : 0.25 for j in 1:nsites]);
              maxdim = 8, cutoff = 0.0)
    charge(state) = sum(real.(operator_expectation_profile(
      state, [(x, sz) for x in 1:nsites]; normalize = false)))
    initial_charge = charge(rho)
    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule = schedule, reverse_schedule = reverse(schedule),
                           nstep = 1, maxdim = maxdim, cutoff = 1e-14, normalize = false)
    for _ in 1:10
      dmt_evolve!(rho, evo)
      @test abs(charge(rho) - initial_charge) < 1e-9
      @test maximum(dim(linkind(rho, b)) for b in 1:(nsites - 1)) <= maxdim
    end
  end

  @testset "exact-diagonalization oracle at d = 3" begin
    # With a budget large enough that no truncation fires, DMT evolution must equal dense
    # Liouvillian evolution of the same Trotter circuit.
    d, nsites, dt, nstep = 3, 4, 0.05, 3
    sites = operator_siteinds(nsites; d = d)
    sz = ComplexF64[1 0 0; 0 0 0; 0 0 -1]
    sx = ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
    h = kron(sx, sx) + kron(sz, sz)
    gate = operator_gate_from_hamiltonian(h, dt; d = d)
    identity_matrix = Matrix{ComplexF64}(I, d, d)
    rho0 = add(operator_product_state(sites, fill(identity_matrix, nsites)),
               operator_product_state(sites, [j == 2 ? sz : identity_matrix for j in 1:nsites]);
               maxdim = 4, cutoff = 0.0)

    dense = kron(kron(kron(identity_matrix, identity_matrix), identity_matrix), identity_matrix) +
            kron(kron(kron(identity_matrix, sz), identity_matrix), identity_matrix)
    u2 = exp(-im * dt * Matrix{ComplexF64}(h))
    function embed(u, start)
      left = start == 1 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^(start - 1), d^(start - 1))
      right_sites = nsites - start - 1
      right = right_sites == 0 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^right_sites, d^right_sites)
      return kron(kron(left, u), right)
    end

    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, dt; schedule = schedule, reverse_schedule = reverse(schedule),
                           nstep = 1, maxdim = 400, cutoff = 1e-14, normalize = false)
    rho = copy(rho0)
    for _ in 1:nstep
      dmt_evolve!(rho, evo)
      for start in schedule
        big = embed(u2, start)
        dense = big * dense * big'
      end
      for start in reverse(schedule)
        big = embed(u2, start)
        dense = big * dense * big'
      end
    end

    for (start, op) in ((1, sz), (2, sz), (3, sz), (1, kron(sz, sz)))
      span = size(op, 1) == d ? 1 : 2
      left = start == 1 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^(start - 1), d^(start - 1))
      right_sites = nsites - start - span + 1
      right = right_sites == 0 ? ComplexF64[1] : Matrix{ComplexF64}(I, d^right_sites, d^right_sites)
      padded = kron(kron(left, ComplexF64.(op)), right)
      @test operator_expectation(rho, op, start; normalize = false) ≈ tr(dense * padded) atol = 1e-8
    end
  end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_dmt_higher_spin.jl`
Expected: FAIL or ERROR. Likely causes, in order of probability: the `combiner`-based protected
block in `_dmt_bond_truncate!` orders indices differently than `matrix(...)` expects when
`radius > 1`; `_left_env_at!` / `_right_env_at!` watermark arithmetic assumes `radius == 1`.

- [ ] **Step 3: Fix what the tests expose**

Work through the failures one at a time. The two known-risky spots:

1. **Environment cache watermarks.** `_invalidate_env!(cache, lo, hi)` uses a one-site margin
   (`lo - 2`, `hi + 2`) justified for `radius = 1`. With `radius` sites folded into the protected
   block, the environments stop at `bond - radius` and `bond + 1 + radius`, so the invalidation
   range must widen correspondingly: call
   `_invalidate_env!(cache, bond - left_count, bond + 1 + right_count)` (already in the Task 6
   code) and confirm with the existing verification flag:

```julia
MPSToolkit._DMT_VERIFY_ENVS[] = true
# rerun the higher-spin sweep testset; it throws if any cached env differs from the rebuild
```

2. **Combiner index order.** `combinedind(left_combiner)` fuses the site indices in an
   unspecified order. That is fine — the protected *span* is all that matters, and `_protected_basis`
   only needs the column space — but the identity direction `protected_left[:, 1]` must still be
   the all-identity column. Assert it explicitly rather than assuming:

```julia
# in the kernel, immediately after building protected_left / protected_right
@assert size(protected_left, 2) == d^(2 * left_count)
```

   and pick `q0` by contracting the identity cap directly rather than by column index if the
   combiner reorders. The safe construction is to build `q0` from a separate contraction with the
   identity cap on all `left_count` sites, which is exactly `_left_identity_environment(psi, bond)`
   applied to `u`.

- [ ] **Step 4: Run the full suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/operator_space/dmt/bond.jl test/test_dmt_higher_spin.jl test/runtests.jl
git commit -m "feat(dmt): validate the faithful kernel at d = 3, 4 and preserve_diameter = 5"
```

---

### Task 8: Migrate examples and remaining call sites

**Files:**
- Modify: `examples/dmt/domain_wall_melting.jl`, `examples/dmt/mixed_field_ising_diffusion.jl`, `examples/dmt/pxp_energy_melting.jl`
- Modify: `examples/operator_space/pxp_energy_correlator.jl`
- Modify: `test/test_core.jl`, `test/test_pxp.jl` (wherever they pass `connector_buffer`)

**Interfaces:**
- Consumes: the Task 6 API.
- Produces: no new API.

- [ ] **Step 1: Find every remaining call site**

Run: `grep -rn "connector_buffer" --include=*.jl --include=*.ipynb .`
Expected: hits only in `examples/`, `test/test_core.jl`, `test/test_pxp.jl`, and the migration
guards themselves.

- [ ] **Step 2: Update each one**

For each hit, delete the `connector_buffer=...` keyword. Then check the `maxdim` at that call
site against the new floor `2 d^2 + 1 = 9` (all examples are `d = 2`) and raise it if needed.
Because `maxdim` now includes the 8 protected directions, an example that ran at `maxdim = 32`
with `connector_buffer = 8` now has 24 complement directions rather than 24 — unchanged in
practice — but an example at `maxdim = 16` now has only 8. Where the spec's
"re-measure, don't assume" rule applies, run the example and record the new numbers.

- [ ] **Step 3: Re-run the examples and record the new outputs**

Run each of the three DMT examples and capture stdout:

```bash
for f in examples/dmt/domain_wall_melting.jl examples/dmt/mixed_field_ising_diffusion.jl examples/dmt/pxp_energy_melting.jl; do
  echo "=== $f ==="; julia --project=. "$f" 2>&1 | tail -25
done
```

Compare the reported transport exponents against the values quoted in the file headers and in
`docs/src/manual/dmt.md`. Update every quoted number to what the run actually produces. If an
exponent moves by more than its own quoted uncertainty, say so explicitly in the commit message
rather than quietly editing the number.

- [ ] **Step 4: Run the full suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add examples test
git commit -m "chore(examples): migrate to preserve_diameter and re-measure DMT outputs"
```

---

### Task 9: Stage A performance — QR bond factorization

**Files:**
- Modify: `src/operator_space/dmt/bond.jl`
- Create: `dev/bench_dmt.jl`
- Modify: `test/test_dmt_kernel.jl`

**Interfaces:**
- Consumes: Task 6 kernel.
- Produces: `_dmt_bond_factorize(psi, bond, left_inds) -> (left_isometry, bond_matrix_mul, bond_matrix_adj, bond_values_or_R, left_link, right_link)` — internal only.

The kernel needs orthonormal left and right bases and the bond matrix expressed in them. It does
**not** need the Schmidt form. Replacing `svd(psi[bond], left_inds)` with a QR gives an identical
truncated operator (the two differ by unitaries on both sides, and every step of the rule is
covariant under them) and was measured at 1.06 s versus 4.42 s at `chi = 1600, d = 3` — 86% of the
bond step is this one factorization.

- [ ] **Step 1: Write the equivalence test**

Append to `test/test_dmt_kernel.jl`:

```julia
@testset "QR and SVD bond factorization give the same truncation" begin
  Random.seed!(4242)
  for d in (2, 3)
    nsites, chi = 6, 50
    maxdim = 2 * d^2 + 10
    sites = operator_siteinds(nsites; d = d)
    base = add(operator_basis_state(sites, fill(1, nsites)),
               0.3 * random_mps(ComplexF64, sites; linkdims = chi); maxdim = chi + 1, cutoff = 0.0)
    via_svd = copy(base)
    MPSToolkit._dmt_bond_truncate!(via_svd, 3; maxdim = maxdim, cutoff = 1e-14, factorize = :svd)
    via_qr = copy(base)
    MPSToolkit._dmt_bond_truncate!(via_qr, 3; maxdim = maxdim, cutoff = 1e-14, factorize = :qr)
    @test abs(inner(via_svd, via_qr)) / (norm(via_svd) * norm(via_qr)) ≈ 1.0 atol = 1e-9
    @test dim(linkind(via_qr, 3)) == dim(linkind(via_svd, 3))
  end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_dmt_kernel.jl`
Expected: FAIL — `_dmt_bond_truncate!` has no `factorize` keyword.

- [ ] **Step 3: Implement the QR path**

Add a `factorize::Symbol=:qr` keyword to `_dmt_bond_truncate!`. For `:qr`, replace the
`svd(psi[bond], left_inds)` block with a QR:

```julia
  if factorize === :qr
    q, r = qr(psi[bond], left_inds)
    left_isometry = q
    left_link = commonind(q, r)
    right_link = commonind(r, psi[bond + 1])
    right_block = psi[bond + 1]
    bond_matrix = matrix(r, left_link, right_link)
  else
    u, s, v = svd(psi[bond], left_inds)
    left_isometry = u
    left_link = commonind(u, s)
    right_link = commonind(v, s)
    right_block = v * psi[bond + 1]
    bond_matrix = Matrix(Diagonal(ComplexF64.(diag(matrix(s, left_link, right_link)))))
  end
```

No change is needed in `src/operator_space/dmt/lowrank.jl`: `_dmt_connector` and
`_dmt_complement_ops` already take a bond matrix and dispatch through `_bond_mul` / `_bond_adj`,
whose `Diagonal` method keeps the SVD path at `O(chi)` per product. Under `:qr` the bond matrix
is `UpperTriangular`, so complement matvecs fall back to `O(chi^2)`; that is still strongly
favorable against the `4.42 s -> 1.06 s` saved on the factorization itself.

Pass the QR factor as `UpperTriangular(bond_matrix)` rather than a bare `Matrix`, so BLAS uses
the triangular kernels.

- [ ] **Step 4: Write the benchmark script**

Create `dev/bench_dmt.jl`:

```julia
# Per-bond DMT cost across local dimension, bond dimension, and kernel options.
# Run: julia --project=. dev/bench_dmt.jl
using ITensors, ITensorMPS, LinearAlgebra, MPSToolkit, Printf, Random

Random.seed!(1)
BLAS.set_num_threads(max(Sys.CPU_THREADS ÷ 2, 1))
println("BLAS threads = ", BLAS.get_num_threads())

function bench_bond(d, chi, maxdim; factorize, truncation, reps = 3)
  nsites = 8
  sites = operator_siteinds(nsites; d = d)
  best = Inf
  for _ in 1:reps
    psi = random_mps(ComplexF64, sites; linkdims = chi)
    normalize!(psi)
    elapsed = @elapsed MPSToolkit._dmt_bond_truncate!(
      psi, 4; maxdim = maxdim, cutoff = 1e-14, factorize = factorize, truncation = truncation)
    best = min(best, elapsed)
  end
  return best
end

for d in (2, 3, 4), chi in (400, 800, 1600)
  maxdim = 2 * d^2 + 180
  for factorize in (:svd, :qr), truncation in (:dense, :random)
    t = bench_bond(d, chi, maxdim; factorize = factorize, truncation = truncation)
    @printf("d=%d chi=%-5d maxdim=%-4d factorize=%-4s truncation=%-7s %7.3f s\n",
            d, chi, maxdim, factorize, truncation, t)
  end
end
```

- [ ] **Step 5: Add the allocation guard**

The point of the basis-free kernel is that no `chi x chi` temporary is ever formed. Pin it, so a
future refactor cannot quietly reintroduce one:

```julia
@testset "the kernel allocates well below chi^2" begin
  d, nsites, chi = 2, 6, 400
  maxdim = 2 * d^2 + 60
  sites = operator_siteinds(nsites; d = d)
  psi = random_mps(ComplexF64, sites; linkdims = chi)
  normalize!(psi)
  MPSToolkit._dmt_bond_truncate!(copy(psi), 3; maxdim = maxdim, cutoff = 1e-14)  # warm up
  bytes = @allocated MPSToolkit._dmt_bond_truncate!(psi, 3; maxdim = maxdim, cutoff = 1e-14)
  # a single dense chi x chi ComplexF64 temporary is 16 * chi^2 bytes; allow a generous
  # multiple of the chi x maxdim working set but stay far below repeated chi^2 matrices.
  @test bytes < 8 * 16 * chi * maxdim
end
```

Run it with `truncation = :dense` too and record the number; if `:dense` cannot meet the bound
(it materializes the complement operator by design), scope the test to `truncation = :random`
and say so in a comment.

- [ ] **Step 6: Run tests and the benchmark**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS. The preservation tests must still be below `1e-11` — the QR path changes the
representation, not the guarantee.

Run: `julia --project=. dev/bench_dmt.jl`
Expected: `factorize=:qr` faster than `:svd` at every `(d, chi)` with `chi >= 400`. Record the
table in the commit message.

- [ ] **Step 7: Commit**

```bash
git add src/operator_space/dmt/bond.jl src/operator_space/dmt/lowrank.jl dev/bench_dmt.jl test/test_dmt_kernel.jl
git commit -m "perf(dmt): factorize the bond tensor by QR instead of SVD"
```

---

### Task 10: Stage B performance — randomized complement and fused gate step

**Files:**
- Modify: `src/operator_space/dmt.jl` (`dmt_step!`, `_dmt_window_truncate!`)
- Modify: `src/operator_space/dmt/bond.jl`
- Modify: `dev/bench_dmt.jl`
- Modify: `test/test_dmt_kernel.jl`

**Interfaces:**
- Consumes: Tasks 6 and 9.
- Produces: `truncation=:random` reachable from `DMTOptions` / `DMTGateEvolution` (already plumbed in Task 6); `gate_maxdim=0` meaning "apply the gate exactly".

- [ ] **Step 1: Write the failing test**

Append to `test/test_dmt_kernel.jl`:

```julia
@testset "randomized truncation preserves the guarantee and tracks the dense result" begin
  Random.seed!(31337)
  for d in (2, 3)
    nsites, chi = 6, 60
    maxdim = 2 * d^2 + 16
    sites = operator_siteinds(nsites; d = d)
    probes = diameter_probes(nsites, d, 3)
    base = add(operator_basis_state(sites, fill(1, nsites)),
               0.3 * random_mps(ComplexF64, sites; linkdims = chi); maxdim = chi + 1, cutoff = 0.0)
    before = operator_expectation_profile(base, probes; normalize = false)
    dense_state = copy(base)
    MPSToolkit._dmt_bond_truncate!(dense_state, 3; maxdim = maxdim, cutoff = 1e-14, truncation = :dense)
    random_state = copy(base)
    MPSToolkit._dmt_bond_truncate!(random_state, 3; maxdim = maxdim, cutoff = 1e-14, truncation = :random)
    # the guarantee is independent of how well the complement is approximated
    @test preservation_error(before, operator_expectation_profile(random_state, probes; normalize = false)) < 1e-11
    # and the approximation is good
    overlap = abs(inner(dense_state, random_state)) / (norm(dense_state) * norm(random_state))
    @test overlap > 0.99
  end
end

@testset "gate_maxdim = 0 applies the gate exactly" begin
  d, nsites = 2, 6
  sites = operator_siteinds(nsites; d = d)
  psi = random_mps(ComplexF64, sites; linkdims = 8)
  normalize!(psi)
  gate = operator_gate_from_hamiltonian(
    spinhalf_xyz_bond_hamiltonian(; Jx = 1.0, Jy = 0.5, Jz = 0.3), 0.1; d = d)
  exact = copy(psi)
  dmt_step!(exact, gate, 3; maxdim = 200, gate_maxdim = 0, cutoff = 1e-14)
  capped = copy(psi)
  dmt_step!(capped, gate, 3; maxdim = 200, gate_maxdim = 4096, cutoff = 1e-14)
  @test abs(inner(exact, capped)) / (norm(exact) * norm(capped)) ≈ 1.0 atol = 1e-10
end
```

`diameter_probes` and `preservation_error` come from `test/dmt_test_helpers.jl` (created in
Task 6); add `include("dmt_test_helpers.jl")` at the top of `test_dmt_kernel.jl`.

- [ ] **Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_dmt_kernel.jl`
Expected: FAIL — `gate_maxdim = 0` is currently rejected by the `gate_maxdim >= 1` validation.

- [ ] **Step 3: Implement exact gate application**

In `src/operator_space/dmt.jl`, allow `gate_maxdim == 0` in `dmt_step!`, `DMTOptions` and
`DMTGateEvolution` validation, meaning "no cap". Where `dmt_step!` currently calls

```julia
  tebd_evolve!(psi, gate, start; maxdim=Int(gate_maxdim), cutoff=0.0)
```

use

```julia
  gate_budget = Int(gate_maxdim) == 0 ? typemax(Int) : Int(gate_maxdim)
  tebd_evolve!(psi, gate, start; maxdim=gate_budget, cutoff=0.0)
```

Change the default `gate_maxdim` in `DMTOptions` and `DMTGateEvolution` from
`max(maxdim * 16, 64)` to `0`, and update their docstrings: the pre-truncation is a plain SVD
that discards small singular values *before* DMT sees them, which is exactly the error DMT
exists to avoid, so the default is now to not do it.

- [ ] **Step 4: Benchmark the combinations and pick the defaults**

Extend `dev/bench_dmt.jl` with a full-sweep benchmark that measures gate application and DMT
truncation separately:

```julia
function bench_sweep(d, nsites, maxdim, gate_maxdim; factorize, truncation, nstep = 1)
  sites = operator_siteinds(nsites; d = d)
  sx = d == 2 ? ComplexF64[0 1; 1 0] / 2 : ComplexF64[0 1 0; 1 0 1; 0 1 0] / sqrt(2)
  sz = d == 2 ? ComplexF64[1 0; 0 -1] / 2 : ComplexF64[1 0 0; 0 0 0; 0 0 -1]
  gate = operator_gate_from_hamiltonian(kron(sx, sx) + kron(sz, sz), 0.1; d = d)
  identity_matrix = Matrix{ComplexF64}(I, d, d)
  rho = add(operator_product_state(sites, fill(identity_matrix, nsites)),
            operator_local_sum_state(sites, sz, [j <= nsites ÷ 2 ? -0.25 : 0.25 for j in 1:nsites]);
            maxdim = 8, cutoff = 0.0)
  schedule = collect(1:(nsites - 1))
  evo = DMTGateEvolution(gate, 0.1; schedule = schedule, reverse_schedule = reverse(schedule),
                         nstep = nstep, maxdim = maxdim, gate_maxdim = gate_maxdim,
                         cutoff = 1e-14, truncation = truncation, normalize = false)
  return @elapsed dmt_evolve!(rho, evo)
end
```

Run it for `d in (2, 3)`, `maxdim in (2 d^2 + 100, 2 d^2 + 300)`, `gate_maxdim in (0, 2 * maxdim)`,
`truncation in (:dense, :random)`. Choose the shipped default for `truncation` from the measured
table: `:dense` if `:random` is not clearly faster at these sizes, `:random` otherwise. Record
the table in the commit message and set the default accordingly in `DMTOptions` and
`DMTGateEvolution`.

- [ ] **Step 5: Run the full suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS, with the preservation tests unchanged at `< 1e-11`.

- [ ] **Step 6: Commit**

```bash
git add src/operator_space/dmt.jl src/operator_space/dmt/bond.jl dev/bench_dmt.jl test/test_dmt_kernel.jl test/dmt_test_helpers.jl
git commit -m "perf(dmt): randomized complement truncation and exact gate application"
```

---

### Task 11: Documentation

**Files:**
- Modify: `docs/src/manual/dmt.md`
- Modify: `docs/src/manual/operator-space.md`
- Modify: `docs/src/manual/transport-methods.md`
- Modify: `docs/src/api.md`
- Modify: `README.md` if it mentions Pauli-only operator space

- [ ] **Step 1: Rewrite the DMT manual's contract sections**

In `docs/src/manual/dmt.md`:

- Replace the "identity-preserving truncation rule" section with the actual rule: the first
  `d^2` rows and columns of the bond matrix are preserved exactly, so every observable of
  diameter `<= preserve_diameter` survives truncation exactly. State the `2 d^(2n) + chi'` rank
  bound and cite arXiv:1707.01506 Sec. III.
- Replace the `connector_buffer` bullet in the knobs list with `preserve_diameter`, and add the
  budget table:

```markdown
| `preserve_diameter` | protected block `2 d^(2n)` | minimum `maxdim`, `d = 2` | `d = 3` | `d = 4` |
|:--|--:|--:|--:|--:|
| 3 | `2 d^2`  | 9  | 19  | 33  |
| 5 | `2 d^4`  | 33 | 163 | 513 |
```

- Add a **Migration** admonition: `connector_buffer` was removed; `maxdim` is now the total
  bond dimension inclusive of the protected block, so a run at fixed `maxdim` keeps
  `maxdim - 2 d^(2n)` complement directions.
- Delete the "Pauli sites only" bullet from Tips and pitfalls; replace it with a higher-spin
  note: the protected block grows as `d^(2n)`, so spin-1 needs a substantially larger budget
  than spin-1/2 for the same complement resolution (arXiv:2205.02853 ran `d = 3` at
  `chi = 128-256`).
- Add a short **Higher spin** section with a runnable snippet:

```julia
using MPSToolkit, ITensors, ITensorMPS
using EDKit: spin                       # any source of dense local Hamiltonians works

d, nsites, dt, maxdim = 3, 40, 0.1, 64
sites = operator_siteinds(nsites; d = d)

h = spin((1, "xx"), (1, "yy"), (1, "zz"); D = d)      # sparse input is accepted
gate = operator_gate_from_hamiltonian(h, dt; d = d)

sz = Matrix(spin("z"; D = d))
identity_matrix = Matrix{ComplexF64}(I, d, d)
rho = add(operator_product_state(sites, fill(identity_matrix, nsites)),
          operator_local_sum_state(sites, sz, [j <= nsites ÷ 2 ? -0.25 : 0.25 for j in 1:nsites]);
          maxdim = 8, cutoff = 0.0)

schedule = collect(1:(nsites - 1))
evo = DMTGateEvolution(gate, dt; schedule, reverse_schedule = reverse(schedule),
                       maxdim = maxdim, cutoff = 1e-12, normalize = false)
profile(state) = real.(operator_expectation_profile(state, [(x, sz) for x in 1:nsites];
                                                    normalize = false))
```

  State the EDKit convention explicitly: `spin(...; D = d)` composes multi-site strings with
  `kron` left-to-right, first character = leftmost site, matching MPSToolkit.

- [ ] **Step 2: Update the operator-space manual**

In `docs/src/manual/operator-space.md`, document the `operator_*` API as primary and `pauli_*` as
the `d = 2` case; document the generalized Gell-Mann basis (ordering, identity first, exact
reduction to `(I, X, Y, Z)/sqrt(2)`); and add the `d^(N/2)` trace-convention note replacing
`2^(N/2)`. Add the warning that at `d >= 3` a physical `S^z` is not a single basis element, so
`operator_product_state` / `operator_local_sum_state` are the builders to use.

- [ ] **Step 3: Update the remaining docs**

Add the new exported symbols to `docs/src/api.md`. In
`docs/src/manual/transport-methods.md`, update any `connector_buffer` mention. Grep to confirm:

```bash
grep -rn "connector_buffer\|dimension 4\|dimension-4\|Pauli sites only" docs README.md
```

Expected: no hits outside the migration notes.

- [ ] **Step 4: Build the docs and run the doc tests**

Run: `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'`
Expected: builds without errors; doctests pass.

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS, including `test_docs.jl` (which checks documentation portability) and
`test_docstrings.jl`.

- [ ] **Step 5: Commit**

```bash
git add docs README.md
git commit -m "docs: higher-spin DMT contract, budget semantics and EDKit interop"
```

---

## Verification checklist

Before opening a PR, confirm each spec success criterion:

1. `julia --project=. test/test_dmt_preservation.jl` and `test/test_dmt_higher_spin.jl` pass —
   diameter-`<=3` observables preserved below `1e-11` at `d = 2, 3, 4`, Hermitian and
   non-Hermitian.
2. The spin-1 Heisenberg melt testset conserves total `S^z` below `1e-9` with the bond dimension
   pinned at the budget.
3. `grep -rn "connector_buffer" --include=*.jl src` returns only the migration guards; every
   `pauli_*` name still exists and still works.
4. `julia --project=. dev/bench_dmt.jl` output is recorded, with Stage A (`:qr`) and Stage B
   (`:random`, `gate_maxdim = 0`) contributions separable.
5. `julia --project=. -e 'using Pkg; Pkg.test()'` is fully green.
