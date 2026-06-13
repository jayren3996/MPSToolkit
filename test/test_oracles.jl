using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Test

# Exact-diagonalization oracles for the core numerics. These guard the time-stepping
# (Trotter sign/factors, 2nd-order Strang, the TDVP `t = -i*time` convention), the periodic
# wraparound gate-ordering convention, and the Chebyshev/KPM normalization against
# independent dense references rather than internal self-consistency.

# Contract an MPS into a dense statevector with a fixed (big-endian) index order.
function _statevector(psi)
  sites = [siteind(psi, n) for n in 1:length(psi)]
  contracted = psi[1]
  for n in 2:length(psi)
    contracted = contracted * psi[n]
  end
  arr = Array(contracted, sites...)
  return vec(permutedims(arr, reverse(1:ndims(arr))))
end

# Embed a 2-site bond operator `h2` at bond `b` of an `N`-site chain (matching `_statevector`).
function _embed_bond(h2, b, N)
  id = Matrix{ComplexF64}(I, 2, 2)
  left = b == 1 ? Matrix{ComplexF64}(I, 1, 1) : foldl(kron, fill(id, b - 1))
  right_count = N - (b + 1)
  right = right_count == 0 ? Matrix{ComplexF64}(I, 1, 1) : foldl(kron, fill(id, right_count))
  return kron(left, h2, right)
end

# Gauge-invariant L2 distance between two statevectors (removes the global phase).
function _gauge_l2(reference, candidate)
  ref = reference / norm(reference)
  cand = candidate / norm(candidate)
  phase = angle(dot(ref, cand))
  return norm(cand - ref * exp(im * phase))
end

function _tfim_dense(N; J=1.0, g=1.3)
  Hfull = zeros(ComplexF64, 2^N, 2^N)
  for b in 1:(N - 1)
    Hfull .+= _embed_bond(spinhalf_tfim_bond_hamiltonian(N, b; J=J, g=g), b, N)
  end
  return Hermitian(Hfull)
end

function _tfim_mpo(sites; J=1.0, g=1.3)
  opsum = OpSum()
  for j in 1:(length(sites) - 1)
    opsum += -J, "Sz", j, "Sz", j + 1
  end
  for j in 1:length(sites)
    opsum += -g, "Sx", j
  end
  return MPO(opsum, sites)
end

@testset "exact-diagonalization oracles" begin
  @testset "TEBD matches exp(-iHt) and is 2nd-order (TFIM N=4)" begin
    N = 4
    J = 1.0
    g = 1.3
    Ttot = 0.2
    Hfull = _tfim_dense(N; J=J, g=g)
    sites = siteinds("S=1/2", N)
    init() = MPS(sites, n -> isodd(n) ? "Up" : "Dn")
    v_exact = exp(-im * Matrix(Hfull) * Ttot) * _statevector(init())
    local_ham = (bond, weight) -> weight * spinhalf_tfim_bond_hamiltonian(N, bond; J=J, g=g)

    errors = Float64[]
    for nsteps in (8, 16)
      psi = init()
      evolve!(psi, tebd_strang_evolution(N, Ttot / nsteps; local_hamiltonian=local_ham,
                                         nstep=nsteps, maxdim=64, cutoff=1e-14))
      push!(errors, _gauge_l2(v_exact, _statevector(psi)))
    end

    @test errors[2] < 1e-4            # faithful at the finer step (validates the -im sign)
    @test errors[1] / errors[2] > 3.0 # halving dt cuts the error ~4x => 2nd-order Strang
  end

  @testset "TDVP matches exp(-iHt) with t = -i*time (TFIM N=4)" begin
    N = 4
    sites = siteinds("S=1/2", N)
    H = _tfim_mpo(sites; J=1.0, g=1.3)
    Hfull = _tfim_dense(N; J=1.0, g=1.3)
    init() = MPS(sites, n -> isodd(n) ? "Up" : "Dn")
    time = 0.2
    nsteps = 20
    v_exact = exp(-im * Matrix(Hfull) * time) * _statevector(init())
    v_wrong = exp(+im * Matrix(Hfull) * time) * _statevector(init())

    psi = init()
    evolve!(psi, TDVPEvolution(H, -im * time; time_step=-im * time / nsteps, nsteps=nsteps,
                               solver_kwargs=(maxdim=32, cutoff=1e-14)))
    evolved = _statevector(psi)
    @test _gauge_l2(v_exact, evolved) < 1e-3
    @test _gauge_l2(v_wrong, evolved) > 1e-1   # the opposite sign must NOT match
  end

  @testset "periodic wraparound gate applies kron(op_N, op_1)" begin
    N = 4
    sites = siteinds("S=1/2", N)
    X = ComplexF64[0 1; 1 0]
    Id = Matrix{ComplexF64}(I, 2, 2)
    flipped(gate) = begin
      psi = MPS(sites, n -> "Up")
      tebd_evolve!(psi, gate, N; maxdim=16, cutoff=1e-14)
      findall(z -> z < 0, expect(psi, "Sz"))
    end
    # At the wraparound bond (bond == N) the first kron factor acts on physical site N and
    # the second on physical site 1.
    @test flipped(kron(X, Id)) == [N]
    @test flipped(kron(Id, X)) == [1]
  end

  @testset "Chebyshev KPM sum rule recovers mu_0" begin
    N = 4
    sites = siteinds("S=1/2", N)
    opsum = OpSum()
    for j in 1:(N - 1)
      opsum += "Sz", j, "Sz", j + 1
    end
    for j in 1:N
      opsum += 0.5, "Sx", j
    end
    H = MPO(opsum, sites) / 3.0   # rescale the spectrum strictly inside (-1, 1)
    psi = MPS(sites, n -> isodd(n) ? "Up" : "Dn")
    order = 60
    moments = chebyshev_moments(H, psi; order=order, maxdim=64, cutoff=0.0)

    # Chebyshev-Gauss quadrature is exact for the 1/(pi*sqrt(1-x^2)) KPM measure.
    nq = 4 * order
    nodes = [cos((j - 0.5) * pi / nq) for j in 1:nq]
    integral = sum(reconstruct_chebyshev(x, moments; kernel=nothing) * pi * sqrt(1 - x^2) for x in nodes) / nq
    @test integral ≈ moments[1] atol = 1e-8
    @test moments[1] ≈ real(inner(psi, psi)) atol = 1e-10

    # The m=0 sum rule is blind to the n>=1 terms: int T_n/(pi*sqrt(1-x^2)) = 0 for n>=1, so a
    # wrong factor on the `2 * sum` would still integrate to mu_0. Project the reconstructed
    # density onto T_m by the same Chebyshev-Gauss quadrature -- (1/nq) sum_q g(x_q) T_m(x_q)
    # recovers mu_m for m>=1 -- to pin the factor-of-2 and the kernel=nothing reconstruction.
    for m in (1, 2, 5)
      projected = sum(reconstruct_chebyshev(x, moments; kernel=nothing) * pi * sqrt(1 - x^2) * cos(m * acos(x)) for x in nodes) / nq
      @test projected ≈ moments[m + 1] atol = 1e-8
    end
  end
end
