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

function _validate_vectorization_sites(op::MPO, sites)
  length(op) == length(sites) || throw(ArgumentError("MPO and Pauli sites must have matching lengths"))
  for site in sites
    dim(site) == 4 || throw(ArgumentError("Pauli operator-space sites must have local dimension 4"))
  end
  return nothing
end

"""
    pauli_state_from_mpo(op, sites)

Vectorize a physical-space `MPO` on spin-1/2 sites into a Pauli-basis operator-space `MPS`.

The output amplitudes follow the normalized-Pauli convention used throughout the
operator-space helpers: the coefficient of the basis product `α = (α_1, …, α_N)` is
`tr(P_α† O)` with `P_α = ⊗_j σ_{α_j} / √2`.

# Arguments
- `op`: Physical `MPO` acting on local dimension-2 sites.
- `sites`: Pauli-space site indices (dimension 4), typically from [`pauli_siteinds`](@ref).

# Returns
- An `MPS` over `sites` with the same bond dimensions as `op`.

# Notes
- The converter reuses the MPO link indices, so a bond-dimension-`k` MPO vectorizes at bond
  dimension `k`; for example the PXP constraint projector vectorizes at bond dimension 2
  ([`pauli_pxp_constraint_state`](@ref)).
"""
function pauli_state_from_mpo(op::MPO, sites)
  _validate_vectorization_sites(op, sites)
  basis = _pauli_basis_operators(1)
  tensors = ITensor[]
  for n in 1:length(op)
    phys = _mpo_unprimed_site_index(op, n)
    dim(phys) == 2 || throw(ArgumentError("pauli_state_from_mpo requires physical dimension-2 sites"))
    conversion = ITensor(ComplexF64, sites[n], prime(phys), phys)
    for (alpha, pauli) in enumerate(basis)
      for row in 1:2, column in 1:2
        value = conj(pauli[row, column])
        iszero(value) && continue
        conversion[sites[n] => alpha, prime(phys) => row, phys => column] = value
      end
    end
    push!(tensors, op[n] * conversion)
  end
  return MPS(tensors)
end

"""
    pauli_superoperator_mpo(op, sites)

Lift a physical-space `MPO` `M` on spin-1/2 sites to the operator-space `MPO` implementing
the two-sided sandwich `ρ ↦ M ρ M†` in the normalized Pauli basis.

# Arguments
- `op`: Physical `MPO` acting on local dimension-2 sites.
- `sites`: Pauli-space site indices (dimension 4).

# Returns
- An `MPO` over `sites` with prime convention `(prime(site), site)`, directly usable with
  `ITensorMPS.apply`. A bond-dimension-`k` physical MPO lifts to bond dimension `k^2`.

# Notes
- The local tensor is `S[α', α] = tr(P_α'† M_loc P_α M_loc'†)` contracted over a ket copy and
  a conjugated bra copy of each MPO tensor, with the (ket, bra) link pairs fused by combiners.
- For a Hermitian projector `P` this gives the superoperator `ρ ↦ P ρ P`; see
  [`pauli_pxp_constraint_projector`](@ref).
"""
function pauli_superoperator_mpo(op::MPO, sites)
  _validate_vectorization_sites(op, sites)
  nsites = length(op)
  basis = _pauli_basis_operators(1)
  ket_links = [commonind(op[n], op[n + 1]) for n in 1:(nsites - 1)]
  bra_links = [sim(link) for link in ket_links]

  tensors = ITensor[]
  for n in 1:nsites
    phys = _mpo_unprimed_site_index(op, n)
    dim(phys) == 2 || throw(ArgumentError("pauli_superoperator_mpo requires physical dimension-2 sites"))
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
    for (alpha, pauli) in enumerate(basis)
      for row in 1:2, column in 1:2
        value = pauli[row, column]
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
