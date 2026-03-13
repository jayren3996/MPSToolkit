"""
    chebyshev_moments(H, psi; order, maxdim=32, cutoff=1e-12, normalize_initial=true)

Compute Chebyshev moments `μ_n = ⟨ψ|T_n(H)|ψ⟩` for a finite `MPO` and `MPS`.

`H` is assumed to already be rescaled so its spectrum lies in `[-1, 1]`.
The returned vector is stored in the natural Julia order:

- `moments[1] = μ_0`
- `moments[2] = μ_1`
- ...
"""
function chebyshev_moments(
  H::MPO,
  psi::MPS;
  order::Integer,
  maxdim::Integer=32,
  cutoff::Real=1e-12,
  normalize_initial::Bool=true,
)
  order > 0 || throw(ArgumentError("Chebyshev order must be positive"))
  hascommoninds(siteinds, H, psi) || throw(ArgumentError("Hamiltonian and state must share site indices"))
  hascommoninds(siteinds, H, psi') || throw(ArgumentError("Hamiltonian and bra state must share site indices"))

  base_state = normalize_initial ? normalize(psi) : copy(psi)
  moments = Vector{Float64}(undef, order)
  moments[1] = real(inner(base_state, base_state))
  order == 1 && return moments

  tm_prev = base_state
  tm_curr = apply(H, base_state; maxdim=maxdim, cutoff=cutoff)
  moments[2] = real(inner(base_state, tm_curr))
  order == 2 && return moments

  doubled_h = 2 * H
  for n in 3:order
    previous_state = -tm_prev
    next_state = add(apply(doubled_h, tm_curr; maxdim=maxdim, cutoff=cutoff), previous_state; maxdim=maxdim, cutoff=cutoff)
    _optimize_chebyshev_vector!(next_state, doubled_h, tm_curr, previous_state)
    tm_prev = tm_curr
    tm_curr = next_state
    moments[n] = real(inner(base_state, tm_curr))
  end
  return moments
end

"""
    _optimize_chebyshev_vector!(psi, h2, current, previous; nsweep=5)

Refine the fitted Chebyshev vector `psi ≈ 2H * current - previous` by sweeping.
"""
function _optimize_chebyshev_vector!(psi::MPS, h2::MPO, current::MPS, previous::MPS; nsweep::Integer=5)
  nsites = length(psi)
  orthogonalize!(psi, 1)

  left_right_h = Vector{ITensor}(undef, nsites)
  left_right_prev = Vector{ITensor}(undef, nsites)
  env_h = ITensor(1.0)
  env_prev = ITensor(1.0)
  for site in nsites:-1:2
    left_right_h[site] = env_h = conj(psi[site]) * noprime(env_h * h2[site] * current[site])
    left_right_prev[site] = env_prev = conj(psi[site]) * env_prev * previous[site]
  end

  for _ in 1:nsweep
    psi[1] = noprime(current[1] * h2[1] * left_right_h[2]) + contract(previous[1], left_right_prev[2])
    orthogonalize!(psi, 2)
    left_right_h[1] = conj(psi[1]) * noprime(h2[1] * current[1])
    left_right_prev[1] = conj(psi[1]) * previous[1]

    for site in 2:(nsites - 1)
      psi[site] = noprime(left_right_h[site - 1] * current[site] * h2[site] * left_right_h[site + 1]) +
                  contract(left_right_prev[site - 1], previous[site], left_right_prev[site + 1])
      orthogonalize!(psi, site + 1)
      left_right_h[site] = conj(psi[site]) * noprime(left_right_h[site - 1] * h2[site] * current[site])
      left_right_prev[site] = contract(left_right_prev[site - 1], conj(psi[site]), previous[site])
    end

    psi[nsites] = noprime(left_right_h[nsites - 1] * current[nsites] * h2[nsites]) + contract(left_right_prev[nsites - 1], previous[nsites])
    orthogonalize!(psi, nsites - 1)
    left_right_h[nsites] = conj(psi[nsites]) * noprime(current[nsites] * h2[nsites])
    left_right_prev[nsites] = contract(conj(psi[nsites]), previous[nsites])

    for site in (nsites - 1):-1:2
      psi[site] = noprime(left_right_h[site - 1] * current[site] * h2[site] * left_right_h[site + 1]) +
                  contract(left_right_prev[site - 1], previous[site], left_right_prev[site + 1])
      orthogonalize!(psi, site - 1)
      left_right_h[site] = conj(psi[site]) * noprime(left_right_h[site + 1] * h2[site] * current[site])
      left_right_prev[site] = contract(left_right_prev[site + 1], conj(psi[site]), previous[site])
    end
  end
  return psi
end
