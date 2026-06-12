"""
    pauli_gate_from_imaginary_time(h, dbeta)

Build the dense Pauli-basis superoperator implementing one two-sided imaginary-time slice

`ρ ↦ e^{-(dbeta/2) h} ρ e^{-(dbeta/2) h}`

for a dense local Hermitian `h`. This reuses [`pauli_gate`](@ref), whose construction
`G[α', α] = tr(P_α'† A P_α A†)` implements conjugation by any matrix `A`, not just a
unitary; for Hermitian `h` the slice `A = e^{-(dbeta/2) h}` is Hermitian and positive.

# Arguments
- `h`: Dense local Hermitian operator on one or more spin-1/2 sites.
- `dbeta`: Imaginary-time increment for this slice (each side picks up `dbeta / 2`).

# Returns
- A dense operator-space gate suitable for `tebd_evolve!` on Pauli sites.
"""
function pauli_gate_from_imaginary_time(h::AbstractMatrix, dbeta::Real)
  size(h, 1) == size(h, 2) || throw(ArgumentError("imaginary-time Hamiltonian must be square"))
  return pauli_gate(exp(-(dbeta / 2) * Matrix{ComplexF64}(h)))
end

"""
    pauli_gibbs_state(sites, terms, weights; nsteps=4, maxdim=64, cutoff=1e-12, initial_state=nothing)

Prepare the vectorized weighted Gibbs operator

`ρ ∝ e^{-K/2} ρ_0 e^{-K/2}`,  `K = Σ_j w_j h_j`,

by Trotterized imaginary-time TEBD in Pauli operator space.

A uniform weight vector `w_j = β` gives the Gibbs state `e^{-β H/2} ρ_0 e^{-β H/2}`; a
sign-split profile (e.g. `+β` on the left half, `-β` on the right, `0` on wall-straddling
terms) gives the energy domain wall `e^{-β H_L} ⊗ e^{+β H_R}` used as the initial state of
DMT energy-transport runs.

# Arguments
- `sites`: Pauli-space site indices (dimension 4).
- `terms`: Collection of `(start, h)` pairs of dense local Hermitian blocks, e.g. built from
  [`pxp_term_hamiltonian`](@ref) and [`pxp_term_support`](@ref).
- `weights`: One real weight `w_j` per term. Zero-weight terms are skipped.

# Keyword Arguments
- `nsteps`: Number of Trotter steps. Each step applies every term gate once in a forward
  sweep and once in a reverse sweep (symmetrized splitting), so each application sandwiches
  `e^{-(w_j / (4 nsteps)) h_j}` on both sides.
- `maxdim`, `cutoff`: Ordinary TEBD truncation budget used during preparation. Plain SVD
  truncation is appropriate here: at small weights the prepared state stays close to a
  product in operator space.
- `initial_state`: Optional starting operator MPS `ρ_0` over `sites`. Defaults to the
  vectorized identity (the unconstrained infinite-temperature state). PXP transport runs
  should pass [`pauli_pxp_constraint_state`](@ref), the constrained infinite-temperature
  state; every PXP term commutes with the constraint, so preparation stays in the sector up
  to truncation error.

# Returns
- A `normalize!`d operator-space `MPS`. The overall scale is arbitrary; physical
  expectations divide by the trace (see [`pauli_expectation_profile`](@ref)).
"""
function pauli_gibbs_state(
  sites,
  terms,
  weights;
  nsteps::Integer=4,
  maxdim::Integer=64,
  cutoff::Real=1e-12,
  initial_state=nothing,
)
  length(terms) == length(weights) || throw(ArgumentError("pauli_gibbs_state requires one weight per term"))
  nsteps >= 1 || throw(ArgumentError("pauli_gibbs_state requires nsteps >= 1"))
  maxdim >= 1 || throw(ArgumentError("pauli_gibbs_state requires maxdim >= 1"))
  cutoff >= 0 || throw(ArgumentError("pauli_gibbs_state requires cutoff >= 0"))

  rho = if isnothing(initial_state)
    pauli_basis_state(sites, fill(1, length(sites)))
  else
    length(initial_state) == length(sites) || throw(ArgumentError("initial_state length must match sites"))
    for n in 1:length(sites)
      siteind(initial_state, n) == sites[n] || throw(ArgumentError("initial_state site indices must match sites"))
    end
    copy(initial_state)
  end

  starts = Int[]
  gates = Matrix{ComplexF64}[]
  for ((start, h), weight) in zip(terms, weights)
    iszero(weight) && continue
    push!(starts, Int(start))
    # Each gate is applied twice per Trotter step (forward + reverse sweep), so the slice per
    # application is w / (2 nsteps), accumulating to e^{-(w/2) h} on each side overall.
    push!(gates, pauli_gate_from_imaginary_time(h, weight / (2 * nsteps)))
  end

  for _ in 1:nsteps
    for i in eachindex(starts)
      tebd_evolve!(rho, gates[i], starts[i]; maxdim=Int(maxdim), cutoff=cutoff)
    end
    for i in reverse(eachindex(starts))
      tebd_evolve!(rho, gates[i], starts[i]; maxdim=Int(maxdim), cutoff=cutoff)
    end
  end
  normalize!(rho)
  return rho
end
