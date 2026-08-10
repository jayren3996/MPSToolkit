"""
    _mpo_unprimed_site_index(op, n)

Return the unprimed physical site index of MPO tensor `n`, inferred structurally (prime level
zero and not shared with a neighboring tensor) so the helper works for MPOs built from any
site-index tags.
"""
function _mpo_unprimed_site_index(op::MPO, n::Int)
  links = Index[]
  n > 1 && push!(links, commonind(op[n], op[n - 1]))
  n < length(op) && push!(links, commonind(op[n], op[n + 1]))
  remaining = filter(i -> plev(i) == 0 && all(l -> i != l, links), collect(inds(op[n])))
  length(remaining) == 1 || throw(ArgumentError("could not infer the unprimed MPO site index at site $(n)"))
  return only(remaining)
end

"""
    _validate_vectorization_sites(op, sites)

Validate that `op` and `sites` have matching lengths and that `sites` share a uniform local
dimension, and return that dimension.
"""
function _validate_vectorization_sites(op::MPO, sites)
  length(op) == length(sites) || throw(ArgumentError("operator MPO and operator-space sites must have matching lengths"))
  d = local_dimension(first(sites))
  for site in sites
    local_dimension(site) == d || throw(ArgumentError("operator-space sites must share a common local dimension"))
  end
  return d
end

"""
    operator_state_from_mpo(op, sites)

Vectorize a physical-space `MPO` on local-dimension-`d` sites into an operator-basis
operator-space `MPS`.

The output amplitudes follow the normalized-basis convention used throughout the operator-space
helpers: the coefficient of the basis product `α = (α_1, …, α_N)` is `tr(P_α† O)` with
`P_α = ⊗_j P_{α_j}`, the multi-site product of [`operator_basis_matrices`](@ref)`(d)`.

# Arguments
- `op`: Physical `MPO` acting on local dimension-`d` sites.
- `sites`: Operator-space site indices (dimension `d^2`), typically from
  [`operator_siteinds`](@ref).

# Returns
- An `MPS` over `sites` with the same bond dimensions as `op`.

# Notes
- The converter reuses the MPO link indices, so a bond-dimension-`k` MPO vectorizes at bond
  dimension `k`; for example the PXP constraint projector vectorizes at bond dimension 2
  ([`pauli_pxp_constraint_state`](@ref)).
"""
function operator_state_from_mpo(op::MPO, sites)
  d = _validate_vectorization_sites(op, sites)
  basis = operator_basis_matrices(d)
  tensors = ITensor[]
  for n in 1:length(op)
    phys = _mpo_unprimed_site_index(op, n)
    dim(phys) == d || throw(ArgumentError("operator_state_from_mpo requires physical dimension-$(d) sites"))
    conversion = ITensor(ComplexF64, sites[n], prime(phys), phys)
    for (alpha, basis_op) in enumerate(basis)
      for row in 1:d, column in 1:d
        value = conj(basis_op[row, column])
        iszero(value) && continue
        conversion[sites[n] => alpha, prime(phys) => row, phys => column] = value
      end
    end
    push!(tensors, op[n] * conversion)
  end
  return MPS(tensors)
end

"""
    pauli_state_from_mpo(op, sites)

Spin-1/2 case of [`operator_state_from_mpo`](@ref): vectorize a physical-space `MPO` on
spin-1/2 sites into a Pauli-basis operator-space `MPS`.

# Arguments
- `op`: Physical `MPO` acting on local dimension-2 sites.
- `sites`: Pauli-space site indices (dimension 4), typically from [`pauli_siteinds`](@ref).

# Returns
- An `MPS` over `sites` with the same bond dimensions as `op`.
"""
pauli_state_from_mpo(op::MPO, sites) = operator_state_from_mpo(op, sites)

"""
    operator_superoperator_mpo(op, sites)

Lift a physical-space `MPO` `M` on local-dimension-`d` sites to the operator-space `MPO`
implementing the two-sided sandwich `ρ ↦ M ρ M†` in the normalized operator basis.

# Arguments
- `op`: Physical `MPO` acting on local dimension-`d` sites.
- `sites`: Operator-space site indices (dimension `d^2`).

# Returns
- An `MPO` over `sites` with prime convention `(prime(site), site)`, directly usable with
  `ITensorMPS.apply`. A bond-dimension-`k` physical MPO lifts to bond dimension `k^2`.

# Notes
- The local tensor is `S[α', α] = tr(P_α'† M_loc P_α M_loc'†)` contracted over a ket copy and
  a conjugated bra copy of each MPO tensor, with the (ket, bra) link pairs fused by combiners.
- For a Hermitian projector `P` this gives the superoperator `ρ ↦ P ρ P`; see
  [`pauli_pxp_constraint_projector`](@ref).
"""
function operator_superoperator_mpo(op::MPO, sites)
  d = _validate_vectorization_sites(op, sites)
  nsites = length(op)
  basis = operator_basis_matrices(d)
  ket_links = [commonind(op[n], op[n + 1]) for n in 1:(nsites - 1)]
  bra_links = [sim(link) for link in ket_links]

  tensors = ITensor[]
  for n in 1:nsites
    phys = _mpo_unprimed_site_index(op, n)
    dim(phys) == d || throw(ArgumentError("operator_superoperator_mpo requires physical dimension-$(d) sites"))
    bra_out = sim(prime(phys))
    bra_in = sim(phys)

    old_indices = Index[prime(phys), phys]
    new_indices = Index[bra_out, bra_in]
    if n > 1
      push!(old_indices, ket_links[n - 1])
      push!(new_indices, bra_links[n - 1])
    end
    if n < nsites
      push!(old_indices, ket_links[n])
      push!(new_indices, bra_links[n])
    end
    bra = replaceinds(dag(op[n]), old_indices, new_indices)

    input_conversion = ITensor(ComplexF64, sites[n], phys, bra_in)
    output_conversion = ITensor(ComplexF64, prime(sites[n]), prime(phys), bra_out)
    for (alpha, basis_op) in enumerate(basis)
      for row in 1:d, column in 1:d
        value = basis_op[row, column]
        iszero(value) && continue
        input_conversion[sites[n] => alpha, phys => row, bra_in => column] = value
        output_conversion[prime(sites[n]) => alpha, prime(phys) => row, bra_out => column] = conj(value)
      end
    end

    push!(tensors, ((op[n] * input_conversion) * bra) * output_conversion)
  end

  for n in 1:(nsites - 1)
    fuse = combiner(ket_links[n], bra_links[n]; tags="OperatorSuperLink,n=$(n)")
    tensors[n] = tensors[n] * fuse
    tensors[n + 1] = tensors[n + 1] * fuse
  end
  return MPO(tensors)
end

"""
    pauli_superoperator_mpo(op, sites)

Spin-1/2 case of [`operator_superoperator_mpo`](@ref): lift a physical-space `MPO` on spin-1/2
sites to the Pauli-basis operator-space `MPO` implementing `ρ ↦ M ρ M†`.

# Arguments
- `op`: Physical `MPO` acting on local dimension-2 sites.
- `sites`: Pauli-space site indices (dimension 4).

# Returns
- An `MPO` over `sites` with prime convention `(prime(site), site)`.
"""
pauli_superoperator_mpo(op::MPO, sites) = operator_superoperator_mpo(op, sites)
