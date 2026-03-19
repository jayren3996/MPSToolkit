"""
    chebyshev_moments(H, psi;
      order,
      maxdim=32,
      cutoff=1e-12,
      normalize_initial=true,
      energy_cutoff=false,
      energy_cutoff_sweeps=5,
      krylovdim=30,
      window=1.0,
      energy_cutoff_tol=1e-12,
      energy_cutoff_verbose=false,
    )

Compute Chebyshev moments `μ_n = ⟨ψ|T_n(H)|ψ⟩` for a finite `MPO` and `MPS`.

`H` is assumed to already be rescaled so its spectrum lies in `[-1, 1]`.
The returned vector is stored in the natural Julia order:

- `moments[1] = μ_0`
- `moments[2] = μ_1`
- ...

# Arguments
- `H`: Rescaled `MPO` whose spectrum should already lie in `[-1, 1]`.
- `psi`: Initial state used to seed the Chebyshev recursion.

# Keyword Arguments
- `order`: Number of moments to compute.
- `maxdim`: Bond dimension cap used during MPO application and vector addition.
- `cutoff`: Truncation cutoff used during MPO application and vector addition.
- `normalize_initial`: If `true`, normalize the starting state before building moments.
- `energy_cutoff`: If `true`, apply [`energy_cutoff!`](@ref) after each recursion step.
- `energy_cutoff_sweeps`: Number of sweeps used by `energy_cutoff!`.
- `krylovdim`: Local Krylov subspace dimension used by `energy_cutoff!`.
- `window`: Allowed rescaled energy window for the cutoff projector.
- `energy_cutoff_tol`: Early-stop tolerance for the cutoff sweeps.
- `energy_cutoff_verbose`: If `true`, print per-sweep cutoff diagnostics.

# Returns
- A dense `Vector{Float64}` storing `μ_0, μ_1, ..., μ_{order-1}`.
"""
function chebyshev_moments(
  H::MPO,
  psi::MPS;
  order::Integer,
  maxdim::Integer=32,
  cutoff::Real=1e-12,
  normalize_initial::Bool=true,
  energy_cutoff::Bool=false,
  energy_cutoff_sweeps::Integer=5,
  krylovdim::Integer=30,
  window::Real=1.0,
  energy_cutoff_tol::Real=1e-12,
  energy_cutoff_verbose::Bool=false,
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
  if energy_cutoff
    energy_cutoff!(
      tm_curr,
      H;
      sweeps=energy_cutoff_sweeps,
      krylovdim=krylovdim,
      window=window,
      tol=energy_cutoff_tol,
      verbose=energy_cutoff_verbose,
    )
  end
  moments[2] = real(inner(base_state, tm_curr))
  order == 2 && return moments

  doubled_h = 2 * H
  for n in 3:order
    previous_state = -tm_prev
    next_state = add(apply(doubled_h, tm_curr; maxdim=maxdim, cutoff=cutoff), previous_state; maxdim=maxdim, cutoff=cutoff)
    _optimize_chebyshev_vector!(next_state, doubled_h, tm_curr, previous_state)
    if energy_cutoff
      energy_cutoff!(
        next_state,
        H;
        sweeps=energy_cutoff_sweeps,
        krylovdim=krylovdim,
        window=window,
        tol=energy_cutoff_tol,
        verbose=energy_cutoff_verbose,
      )
    end
    tm_prev = tm_curr
    tm_curr = next_state
    moments[n] = real(inner(base_state, tm_curr))
  end
  return moments
end

"""
    _optimize_chebyshev_vector!(psi, h2, current, previous; nsweep=5)

Refine the fitted Chebyshev vector `psi ≈ 2H * current - previous` by sweeping.

# Arguments
- `psi`: Current approximation to the next Chebyshev vector. Mutated in place.
- `h2`: The doubled Hamiltonian `2H`.
- `current`: Current Chebyshev vector.
- `previous`: Previous Chebyshev vector, already negated in the caller when appropriate.

# Keyword Arguments
- `nsweep`: Number of left-right/right-left optimization sweeps.

# Returns
- The mutated `psi`.
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

"""
    energy_cutoff!(psi, h; sweeps=5, krylovdim=30, window=1.0, tol=1e-12, verbose=false)

Apply the standalone CheMPS-style energy-window projection to an `MPS` with respect to the
effective local problem induced by the MPO `h`. This is intended for Chebyshev vectors after
the Hamiltonian has already been rescaled into the target Chebyshev window.

# Arguments
- `psi`: State to mutate in place.
- `h`: Rescaled Hamiltonian `MPO`.

# Keyword Arguments
- `sweeps`: Maximum number of left-right/right-left cutoff sweeps.
- `krylovdim`: Local Krylov subspace dimension used at each site update.
- `window`: Allowed rescaled energy window.
- `tol`: Early-stop tolerance on the accumulated sweep error estimate.
- `verbose`: If `true`, print per-sweep diagnostics.

# Returns
- The mutated `psi`.
"""
function energy_cutoff!(
  psi::MPS,
  h::MPO;
  sweeps::Integer=5,
  krylovdim::Integer=30,
  window::Real=1.0,
  tol::Real=1e-12,
  verbose::Bool=false,
)
  projector = ProjMPO(h)
  ITensorMPS.set_nsite!(projector, 1)

  error_estimate = 0.0
  for sweep in 1:Int(sweeps)
    error_estimate = _energy_cutoff_sweep!(projector, psi; krylovdim=Int(krylovdim), window=window)
    if verbose
      println("energy cutoff sweep $sweep: err = $error_estimate")
      flush(stdout)
    end
    error_estimate < tol && break
  end
  return psi
end

"""
    _energy_cutoff_sweep!(projector, psi; krylovdim, window)

Run one bidirectional energy-cutoff sweep and return its RMS error estimate.
"""
function _energy_cutoff_sweep!(projector::ProjMPO, psi::MPS; krylovdim::Integer, window::Real)
  nsites = length(psi)
  accumulated_error = 0.0

  for site in 2:nsites
    orthogonalize!(psi, site)
    position!(projector, psi, site)
    local_error, psi[site] = _krylov_energy_cutoff(projector, psi[site], krylovdim; window=window)
    accumulated_error += local_error
  end

  for site in (nsites - 1):-1:1
    orthogonalize!(psi, site)
    position!(projector, psi, site)
    local_error, psi[site] = _krylov_energy_cutoff(projector, psi[site], krylovdim; window=window)
    accumulated_error += local_error
  end

  return sqrt(accumulated_error / nsites)
end

"""
    _krylov_energy_cutoff(projector, tensor, krylovdim; window=1.0)

Project one local tensor onto the target energy window using a Krylov approximation.

# Returns
- `(error_estimate, projected_tensor)`.
"""
function _krylov_energy_cutoff(projector::ProjMPO, tensor::ITensor, krylovdim::Integer; window::Real=1.0)
  tensor_norm = norm(tensor)
  iszero(tensor_norm) && return 0.0, tensor

  local_krylovdim = min(krylovdim, length(tensor.tensor))
  tridiagonal, basis = _lanczos_tridiagonal(projector, tensor, local_krylovdim)
  eigenvalues, eigenvectors = eigen(tridiagonal)
  projected_coefficients, error_estimate = _project_energy_window(eigenvalues, eigenvectors; window=window)
  projected_tensor = sum(projected_coefficients[index] * basis[index] for index in eachindex(projected_coefficients))
  return error_estimate, tensor_norm * projected_tensor
end

"""
    _project_energy_window(eigenvalues, eigenvectors; window=1.0)

Project a local Krylov basis state onto the eigenmodes whose energies lie within `[-window, window]`.

# Returns
- `(projected_coefficients, removed_weight)`.
"""
function _project_energy_window(eigenvalues::AbstractVector, eigenvectors::AbstractMatrix; window::Real=1.0)
  coefficients = zeros(eltype(eigenvectors), size(eigenvectors, 1))
  coefficients[1] = one(eltype(coefficients))
  removed_weight = 0.0

  for (column, value) in enumerate(eigenvalues)
    if value < -window || value > window
      vector = eigenvectors[:, column]
      overlap = conj(vector[1])
      coefficients .-= overlap * vector
      removed_weight += abs2(overlap)
    end
  end

  return coefficients, removed_weight
end

"""
    _lanczos_tridiagonal(projector, tensor, krylovdim)

Build the Lanczos tridiagonal matrix and basis vectors for one local effective projector problem.

# Returns
- `(tridiagonal_matrix, basis_vectors)`.
"""
function _lanczos_tridiagonal(projector::ProjMPO, tensor::ITensor, krylovdim::Integer)
  basis = Vector{ITensor}(undef, krylovdim)
  diagonal = Vector{Float64}(undef, krylovdim)
  offdiagonal = Vector{Float64}(undef, max(krylovdim - 1, 0))

  basis[1] = normalize(tensor)
  if krylovdim == 1
    diagonal[1] = real(inner(basis[1], projector(basis[1])))
    return SymTridiagonal(diagonal, offdiagonal), basis
  end

  residual = projector(basis[1])
  diagonal[1] = real(inner(basis[1], residual))
  residual -= diagonal[1] * basis[1]
  offdiagonal[1] = norm(residual)
  offdiagonal[1] > 0 || return SymTridiagonal(diagonal[1:1], Float64[]), basis[1:1]
  basis[2] = residual / offdiagonal[1]

  for index in 2:(krylovdim - 1)
    residual = projector(basis[index])
    diagonal[index] = real(inner(basis[index], residual))
    residual -= offdiagonal[index - 1] * basis[index - 1] + diagonal[index] * basis[index]
    offdiagonal[index] = norm(residual)
    offdiagonal[index] > 0 || return SymTridiagonal(diagonal[1:index], offdiagonal[1:(index - 1)]), basis[1:index]
    basis[index + 1] = residual / offdiagonal[index]
  end

  diagonal[krylovdim] = real(inner(basis[krylovdim], projector(basis[krylovdim])))
  return SymTridiagonal(diagonal, offdiagonal), basis
end
