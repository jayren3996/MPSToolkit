using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Test

function _dmt_test_state(nsites)
  sites = pauli_siteinds(nsites)
  labels = [isodd(n) ? "Z" : "I" for n in 1:nsites]
  state = pauli_basis_state(sites, labels)
  normalize!(state)
  return sites, state
end

function _identity_gate(window)
  dim = 4^window
  return Matrix{ComplexF64}(I, dim, dim)
end

function _link_dims(psi)
  return [dim(linkind(psi, bond)) for bond in 1:(length(psi) - 1)]
end

function _dense_pauli_coefficients(psi)
  sites = [siteind(psi, n) for n in 1:length(psi)]
  coeffs = ComplexF64[]
  for labels in Iterators.product(ntuple(_ -> 1:4, length(psi))...)
    push!(coeffs, inner(pauli_basis_state(sites, collect(labels)), psi))
  end
  return coeffs
end

function _manual_dmt_sweep!(psi, gates; maxdim, cutoff, gate_maxdim, connector_buffer)
  for (bond, gate) in enumerate(gates)
    dmt_step!(
      psi,
      gate,
      bond;
      maxdim=maxdim,
      cutoff=cutoff,
      direction=:R,
      gate_maxdim=gate_maxdim,
      connector_buffer=connector_buffer,
    )
  end
  for (offset, gate) in enumerate(reverse(gates))
    bond = length(gates) - offset + 1
    dmt_step!(
      psi,
      gate,
      bond;
      maxdim=maxdim,
      cutoff=cutoff,
      direction=:L,
      gate_maxdim=gate_maxdim,
      connector_buffer=connector_buffer,
    )
  end
  normalize!(psi)
  return psi
end

@testset "operator-space DMT" begin
  @test isdefined(MPSToolkit, :dmt_step!)
  @test isdefined(MPSToolkit, :dmt_evolve!)
  @test isdefined(MPSToolkit, :DMTGateEvolution)
  @test MPSToolkit.OperatorSpace.dmt_step! === MPSToolkit.dmt_step!
  @test MPSToolkit.OperatorSpace.dmt_evolve! === MPSToolkit.dmt_evolve!
  @test MPSToolkit.OperatorSpace.DMTGateEvolution === MPSToolkit.DMTGateEvolution

  @testset "DMTOptions validates truncation budgets" begin
    opts = DMTOptions()
    @test opts.maxdim == 30
    @test opts.cutoff == 1e-12
    @test opts.gate_maxdim == 480
    @test opts.connector_buffer == 8
    # gate_maxdim default follows the maxdim*16 formula shared with dmt_step!/DMTGateEvolution.
    @test DMTOptions(maxdim=8).gate_maxdim == max(8 * 16, 64)

    @test_throws ArgumentError DMTOptions(maxdim=0)
    @test_throws ArgumentError DMTOptions(cutoff=-1e-12)
    @test_throws ArgumentError DMTOptions(gate_maxdim=0)
    @test_throws ArgumentError DMTOptions(connector_buffer=-1)
    @test_throws ArgumentError DMTOptions(maxdim=2, connector_buffer=3)

    # DMTOptions is consumed by dmt_step! and forwards equivalently to the keyword form.
    sites = pauli_siteinds(4)
    via_opts = pauli_basis_state(sites, [2, 3, 4, 2])
    normalize!(via_opts)
    via_kwargs = pauli_basis_state(sites, [2, 3, 4, 2])
    normalize!(via_kwargs)
    gate = _identity_gate(2)
    step_opts = DMTOptions(maxdim=4, cutoff=1e-12, gate_maxdim=16, connector_buffer=2)
    dmt_step!(via_opts, gate, 1, step_opts)
    dmt_step!(via_kwargs, gate, 1; maxdim=4, cutoff=1e-12, gate_maxdim=16, connector_buffer=2)
    @test _dense_pauli_coefficients(via_opts) ≈ _dense_pauli_coefficients(via_kwargs) atol = 1e-12
  end

  @testset "DMT validates Pauli dimension for single-site gates" begin
    sites = siteinds("S=1", 3)   # dim-3 sites are not Pauli operator-space
    psi = MPS(sites, n -> "Up")
    gate = Matrix{ComplexF64}(I, 3, 3)
    @test_throws ArgumentError dmt_step!(psi, gate, 1; maxdim=8, connector_buffer=2)
  end

  @testset "identity DMT step preserves a product operator" begin
    _, psi = _dmt_test_state(4)
    reference = copy(psi)

    dmt_step!(psi, _identity_gate(2), 2; maxdim=8, cutoff=1e-12, direction=:R, gate_maxdim=64)

    @test inner(reference, psi) ≈ 1.0 atol = 1e-10
    @test dmt_step!(psi, _identity_gate(2), 2; maxdim=8, cutoff=1e-12, direction=:R, gate_maxdim=64) === psi
  end

  @testset "DMT truncates an enlarged first bond" begin
    sites = pauli_siteinds(3)
    terms = [
      pauli_basis_state(sites, [1, 1, 1]),
      pauli_basis_state(sites, [2, 2, 1]),
      pauli_basis_state(sites, [3, 3, 1]),
    ]
    psi = normalize(add(terms...; maxdim=8, cutoff=1e-14))

    @test dim(linkind(psi, 1)) > 1
    dmt_step!(psi, _identity_gate(2), 1; maxdim=1, cutoff=1e-12, direction=:R, gate_maxdim=8, connector_buffer=0)

    @test dim(linkind(psi, 1)) <= 1
    @test isfinite(real(inner(psi, psi)))
  end

  @testset "DMT validates connector buffer budget" begin
    @test_throws ArgumentError DMTGateEvolution(_identity_gate(2), 0.1; schedule=[1], maxdim=2, connector_buffer=3)

    _, psi = _dmt_test_state(3)
    @test_throws ArgumentError dmt_step!(psi, _identity_gate(2), 1; maxdim=2, connector_buffer=3)
    @test_throws ArgumentError dmt_step!(psi, _identity_gate(2), 1; maxdim=0)
  end

  @testset "DMT completes protected bases beyond local Pauli dimension" begin
    sites = pauli_siteinds(6)
    psi = random_mps(sites; linkdims=20)
    normalize!(psi)

    @test dim(linkind(psi, 3)) == 20
    dmt_step!(psi, _identity_gate(2), 3; maxdim=12, cutoff=1e-12, gate_maxdim=40, connector_buffer=8)

    @test dim(linkind(psi, 3)) == 12
  end

  @testset "DMT preserves identity and local Pauli data under truncation" begin
    sites = pauli_siteinds(5)
    terms = [
      pauli_basis_state(sites, fill(1, 5); coefficient=1.0),
      pauli_basis_state(sites, [1, 4, 1, 1, 1]; coefficient=0.3),
      pauli_basis_state(sites, fill(2, 5); coefficient=0.2),
    ]
    psi = add(terms...; maxdim=8, cutoff=1e-14)
    probes = [fill(1, 5), [1, 4, 1, 1, 1], [1, 1, 4, 1, 1], [1, 1, 1, 4, 1]]
    before = [inner(pauli_basis_state(sites, labels), psi) for labels in probes]

    dmt_step!(psi, _identity_gate(3), 2; maxdim=1, cutoff=1e-12, direction=:R, gate_maxdim=8, connector_buffer=0)

    after = [inner(pauli_basis_state(sites, labels), psi) for labels in probes]
    @test after ≈ before atol = 1e-12
  end

  @testset "complex DMT projection uses adjoint orthonormal bases" begin
    left_protected = ComplexF64[
      1 1+im
      im 2
      0.5-im -0.25
      0.1 0.3im
    ]
    right_protected = ComplexF64[
      1-im 0.2
      0.7 0.4im
      im 1
      0.3+0.2im -0.1
    ]
    left_basis = MPSToolkit._complete_orthonormal_basis(left_protected, 4)
    right_basis = MPSToolkit._complete_orthonormal_basis(right_protected, 4)
    singular_data = ComplexF64[
      1.0 0.2im 0.1 0.3
      -0.4im 0.8 0.2+0.1im 0.0
      0.1im -0.2 0.5 0.3im
      0.0 0.1 0.2im 0.2
    ]

    reduced = left_basis' * singular_data * right_basis
    repaired = left_basis * reduced * right_basis'

    @test left_basis' * left_basis ≈ Matrix{ComplexF64}(I, 4, 4) atol = 1e-12
    @test right_basis' * right_basis ≈ Matrix{ComplexF64}(I, 4, 4) atol = 1e-12
    @test repaired ≈ singular_data atol = 1e-12
  end

  @testset "orthonormal basis completion is orthonormal and span-preserving (large ambient)" begin
    # Guards the allocation-free BLAS Gram-Schmidt path used at the bond dimensions that dominate
    # DMT runtime. The completion must (1) be orthonormal, (2) keep col(protected) inside its
    # span with the connector direction first, and (3) leave a full square basis a unitary that
    # round-trips any operator exactly.
    for (ambient, ncol, target) in ((64, 3, 64), (128, 4, 40), (32, 0, 32), (48, 5, 20))
      protected = ncol == 0 ? zeros(ComplexF64, ambient, 0) : randn(ComplexF64, ambient, ncol)
      basis = MPSToolkit._complete_orthonormal_basis(protected, target)
      @test size(basis) == (ambient, target)
      @test basis' * basis ≈ Matrix{ComplexF64}(I, target, target) atol = 1e-10
      if ncol > 0 && target >= min(ncol, ambient)
        # Every protected column lies in the span of the returned basis.
        @test norm(protected - basis * (basis' * protected)) < 1e-10
      end
      if target == ambient
        x = randn(ComplexF64, ambient, ambient)
        @test basis * (basis' * x * basis) * basis' ≈ x atol = 1e-10
      end
    end
  end

  @testset "DMT truncates every internal bond of wider update windows" begin
    for span in (1, 2, 3, 4, 5)
      nsites = span + 3
      sites = pauli_siteinds(nsites)
      psi = random_mps(sites; linkdims=16)
      normalize!(psi)
      start = 2

      dmt_step!(psi, _identity_gate(span), start; maxdim=10, cutoff=1e-12, direction=:R, gate_maxdim=32, connector_buffer=4)

      for bond in start:(start + span - 2)
        @test dim(linkind(psi, bond)) <= 10
      end
    end
  end

  @testset "invalid DMT calls do not mutate the state" begin
    sites = pauli_siteinds(5)
    psi = random_mps(sites; linkdims=8)
    normalize!(psi)
    reference = copy(psi)
    reference_dims = _link_dims(reference)

    @test_throws ArgumentError dmt_step!(psi, _identity_gate(3), 2; direction=:bad, maxdim=4, gate_maxdim=16)
    @test _link_dims(psi) == reference_dims
    @test inner(reference, psi) ≈ inner(reference, reference) atol = 1e-12

    @test_throws ArgumentError dmt_step!(psi, _identity_gate(3), 4; direction=:R, maxdim=4, gate_maxdim=16)
    @test _link_dims(psi) == reference_dims
    @test inner(reference, psi) ≈ inner(reference, reference) atol = 1e-12

    @test_throws ArgumentError dmt_step!(psi, _identity_gate(2), length(psi); direction=:R, maxdim=4, gate_maxdim=16)
    @test _link_dims(psi) == reference_dims
    @test inner(reference, psi) ≈ inner(reference, reference) atol = 1e-12

    bad_sites = [Index(3, "NotPauli,n=$n") for n in 1:4]
    bad_psi = random_mps(bad_sites; linkdims=4)
    bad_reference_dims = _link_dims(bad_psi)
    @test_throws ArgumentError dmt_step!(bad_psi, Matrix{ComplexF64}(I, 9, 9), 2; maxdim=2, gate_maxdim=4)
    @test _link_dims(bad_psi) == bad_reference_dims
  end

  @testset "scheduled DMT evolution matches explicit two-site sweep" begin
    _, manual = _dmt_test_state(5)
    scheduled = copy(manual)
    x = ComplexF64[0 1; 1 0]
    zz = ComplexF64[1 0; 0 -1]
    gate = pauli_gate(exp(-0.1im * kron(x, zz)))
    gates = [gate, gate, gate, gate]

    _manual_dmt_sweep!(manual, gates; maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
    @test dmt_evolve!(scheduled, evo) === scheduled

    @test inner(manual, scheduled) ≈ 1.0 atol = 1e-8
  end

  @testset "scheduled DMT evolution matches explicit three-site sweep" begin
    _, manual = _dmt_test_state(6)
    scheduled = copy(manual)
    x = ComplexF64[0 1; 1 0]
    z = ComplexF64[1 0; 0 -1]
    gate = pauli_gate(exp(-0.05im * kron(kron(x, z), x)))
    gates = [gate, gate, gate, gate]

    _manual_dmt_sweep!(manual, gates; maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
    evo = DMTGateEvolution(gates, 0.05; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=4, cutoff=1e-12, gate_maxdim=64, connector_buffer=4)
    dmt_evolve!(scheduled, evo)

    @test inner(manual, scheduled) ≈ 1.0 atol = 1e-8
  end

  @testset "reverse DMT sweep uses gates associated with original bonds" begin
    sites = pauli_siteinds(4)
    state_a = pauli_basis_state(sites, [2, 1, 1, 1])
    state_b = pauli_basis_state(sites, [1, 1, 2, 1])
    manual = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))
    scheduled = copy(manual)

    gate1 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 2.0 for n in 1:16])
    gate2 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 0.5 for n in 1:16])
    gate3 = Matrix{ComplexF64}(I, 16, 16)
    gates = [gate1, gate2, gate3]

    _manual_dmt_sweep!(manual, gates; maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 3], reverse_schedule=[3, 2, 1], maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    dmt_evolve!(scheduled, evo)

    @test abs(inner(manual, scheduled)) ≈ 1.0 atol = 1e-10
  end

  @testset "reverse DMT sweep maps repeated forward schedule entries by position" begin
    sites = pauli_siteinds(3)
    state_a = pauli_basis_state(sites, [2, 1, 1])
    state_b = pauli_basis_state(sites, [1, 2, 1])
    manual = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))
    scheduled = copy(manual)

    gate1 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 2.0 for n in 1:16])
    gate2 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 0.5 for n in 1:16])
    gate3 = Diagonal(ComplexF64[isodd(n) ? 1.0 : 1.5 for n in 1:16])
    gates = [gate1, gate2, gate3]

    for (bond, gate) in zip([1, 2, 1], gates)
      dmt_step!(manual, gate, bond; maxdim=8, cutoff=1e-12, direction=:R, gate_maxdim=16, connector_buffer=0)
    end
    for (bond, gate) in zip([1, 2, 1], reverse(gates))
      dmt_step!(manual, gate, bond; maxdim=8, cutoff=1e-12, direction=:L, gate_maxdim=16, connector_buffer=0)
    end
    normalize!(manual)

    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 1], reverse_schedule=[1, 2, 1], maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    dmt_evolve!(scheduled, evo)

    @test abs(inner(manual, scheduled)) ≈ 1.0 atol = 1e-10
    ambiguous = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 1], reverse_schedule=[2, 1, 1], maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    @test_throws ArgumentError dmt_evolve!(copy(scheduled), ambiguous)
  end

  @testset "reverse DMT sweep maps callable gates by forward position" begin
    sites = pauli_siteinds(4)
    state_a = pauli_basis_state(sites, [2, 1, 1, 1])
    state_b = pauli_basis_state(sites, [1, 1, 2, 1])
    manual = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))
    scheduled = copy(manual)

    sx = ComplexF64[0 1; 1 0]
    sy = ComplexF64[0 -im; im 0]
    sz = ComplexF64[1 0; 0 -1]
    gate1 = pauli_gate(exp(-0.07im * (kron(sx, sz) + 0.3 * kron(sz, sx))))
    gate2 = pauli_gate(exp(-0.11im * (kron(sy, sx) + 0.2 * kron(sz, sz))))
    gate3 = pauli_gate(exp(-0.05im * (kron(sz, sy) + 0.4 * kron(sx, sx))))
    gates = [gate1, gate2, gate3]

    _manual_dmt_sweep!(manual, gates; maxdim=32, cutoff=1e-12, gate_maxdim=64, connector_buffer=0)
    evo = DMTGateEvolution((bond, index) -> gates[index], 0.1; schedule=[1, 2, 3], reverse_schedule=[3, 2, 1], maxdim=32, cutoff=1e-12, gate_maxdim=64, connector_buffer=0)
    dmt_evolve!(scheduled, evo)

    @test _dense_pauli_coefficients(scheduled) ≈ _dense_pauli_coefficients(manual) atol = 1e-10
  end

  @testset "custom repeated reverse DMT schedule rejects ambiguous callable indices" begin
    sites = pauli_siteinds(3)
    state_a = pauli_basis_state(sites, [2, 1, 1])
    state_b = pauli_basis_state(sites, [1, 2, 1])
    scheduled = normalize(add(state_a, state_b; maxdim=4, cutoff=1e-14))

    gate = Diagonal(ComplexF64[isodd(n) ? 1.0 : 1.2 for n in 1:16])
    forward_schedule = [1, 2, 1]
    reverse_schedule = [2, 1, 1]

    evo = DMTGateEvolution((bond, index) -> gate, 0.1; schedule=forward_schedule, reverse_schedule=reverse_schedule, maxdim=8, cutoff=1e-12, gate_maxdim=16, connector_buffer=0)
    @test_throws ArgumentError dmt_evolve!(scheduled, evo)
  end

  @testset "reduced-matrix truncation enforces maxdim for traceless operators" begin
    n = 16
    cb = 4
    chi = 4
    make_block() = ComplexF64[cis(2 * pi * (i - 1) * (j - 1) / n) for i in 1:n, j in 1:n]
    trailing = (cb + 1):n

    # Zero identity overlap (traceless operator): truncation must still happen so that
    # maxdim is enforced rather than silently skipped by an exact `== 0` guard.
    zero_conn = make_block()
    zero_conn[1, 1] = 0.0 + 0.0im
    original = copy(zero_conn)
    MPSToolkit._mat_trunc!(zero_conn, chi; connector_buffer=cb)
    @test norm(zero_conn - original) > 1e-8
    @test rank(zero_conn[trailing, trailing]; atol=1e-8) <= chi

    # Near-singular identity overlap must not blow up the rank-1 connector.
    near_singular = make_block()
    near_singular[1, 1] = 1e-200 + 0.0im
    MPSToolkit._mat_trunc!(near_singular, chi; connector_buffer=cb)
    @test all(isfinite, near_singular)
    @test norm(near_singular) < 1e3

    # A well-conditioned identity direction is still preserved exactly (unchanged behavior):
    # the protected connector rows must survive truncation.
    well = make_block()
    MPSToolkit._mat_trunc!(well, chi; connector_buffer=cb)
    @test well[1:cb, :] ≈ make_block()[1:cb, :] atol = 1e-10
  end

  @testset "multi-bond DMT window matches single-bond sequence (canonical gauge)" begin
    # A span-S gate window truncates its S-1 internal bonds. The faithful DMT result is each of
    # those bonds truncated as an independent single-bond (span-2) DMT update in canonical
    # gauge. This is a regression guard for a gauge bug in which the multi-bond window reused
    # cached environments and truncated later bonds with the orthogonality center left on the
    # wrong site -- mildly wrong for the :R forward sweep and severely wrong for the :L reverse
    # sweep, where the bond SVD degenerated to s ≈ I and the truncation discarded information
    # indiscriminately. random_mps gives generic (high operator-entanglement) states where the
    # bug is visible; simple low-rank states are not a sufficient discriminator.
    function _single_bond_sequence!(psi, start, span, direction; maxdim, cutoff, connector_buffer)
      bonds = collect(start:(start + span - 2))
      direction === :L && (bonds = reverse(bonds))
      for bond in bonds
        MPSToolkit._dmt_window_truncate!(psi, bond, 2; maxdim=maxdim, cutoff=cutoff,
          direction=direction, connector_buffer=connector_buffer)
      end
      return psi
    end

    for direction in (:R, :L), span in (3, 4, 5)
      nsites = span + 3
      sites = pauli_siteinds(nsites)
      psi = random_mps(ComplexF64, sites; linkdims=24)
      normalize!(psi)

      windowed = copy(psi)
      MPSToolkit._dmt_window_truncate!(windowed, 2, span; maxdim=6, cutoff=1e-12,
        direction=direction, connector_buffer=2)

      sequential = copy(psi)
      _single_bond_sequence!(sequential, 2, span, direction; maxdim=6, cutoff=1e-12, connector_buffer=2)

      @test abs(inner(windowed, sequential)) / (norm(windowed) * norm(sequential)) ≈ 1.0 atol = 1e-10
    end
  end
end
