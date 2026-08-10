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

function _manual_dmt_sweep!(psi, gates; maxdim, cutoff, gate_maxdim)
  for (bond, gate) in enumerate(gates)
    dmt_step!(
      psi,
      gate,
      bond;
      maxdim=maxdim,
      cutoff=cutoff,
      direction=:R,
      gate_maxdim=gate_maxdim,
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
    @test opts.gate_maxdim == 0
    @test opts.preserve_diameter == 3
    @test opts.truncation == :dense
    # gate_maxdim defaults to 0 = "apply the gate exactly, no pre-truncation", independently of
    # maxdim, and matches DMTGateEvolution / dmt_step!.
    @test DMTOptions(maxdim=8).gate_maxdim == 0
    @test DMTGateEvolution(_identity_gate(2), 0.1; schedule=[1]).gate_maxdim == 0

    @test_throws ArgumentError DMTOptions(maxdim=0)
    @test_throws ArgumentError DMTOptions(cutoff=-1e-12)
    # 0 is the "no cap" sentinel, so only a negative budget is an error.
    @test DMTOptions(gate_maxdim=0).gate_maxdim == 0
    @test_throws ArgumentError DMTOptions(gate_maxdim=-1)
    @test_throws ArgumentError DMTOptions(preserve_diameter=0)
    @test_throws ArgumentError DMTOptions(preserve_diameter=4)
    @test_throws ArgumentError DMTOptions(truncation=:bogus)

    # DMTOptions is consumed by dmt_step! and forwards equivalently to the keyword form.
    sites = pauli_siteinds(4)
    via_opts = pauli_basis_state(sites, [2, 3, 4, 2])
    normalize!(via_opts)
    via_kwargs = pauli_basis_state(sites, [2, 3, 4, 2])
    normalize!(via_kwargs)
    gate = _identity_gate(2)
    step_opts = DMTOptions(maxdim=12, cutoff=1e-12, gate_maxdim=16)
    dmt_step!(via_opts, gate, 1, step_opts)
    dmt_step!(via_kwargs, gate, 1; maxdim=12, cutoff=1e-12, gate_maxdim=16)
    @test _dense_pauli_coefficients(via_opts) ≈ _dense_pauli_coefficients(via_kwargs) atol = 1e-12
  end

  @testset "generic evolve! dispatch carries the normalize field" begin
    sites = pauli_siteinds(6)
    psi0 = random_mps(sites; linkdims=20)
    normalize!(psi0)
    gates = [_identity_gate(2) for _ in 1:5]
    schedule = collect(1:5)
    mk(nrm) = DMTGateEvolution(
      gates, 0.1; schedule=schedule, reverse_schedule=reverse(schedule),
      nstep=1, maxdim=12, cutoff=1e-12, gate_maxdim=64, normalize=nrm,
    )

    # Generic evolve! returns the mutated MPS and matches dmt_evolve! at the same setting.
    a = copy(psi0)
    @test evolve!(a, mk(true)) === a
    b = copy(psi0)
    dmt_evolve!(b, mk(true))
    @test abs(inner(a, b)) ≈ 1.0 atol = 1e-10

    # The normalize field is honored: truncation sheds norm, restored only when normalize=true.
    c = copy(psi0)
    evolve!(c, mk(false))
    @test norm(a) ≈ 1.0 atol = 1e-10
    @test norm(c) < 0.999
    # An explicit keyword still overrides the field.
    d = copy(psi0)
    evolve!(d, mk(true); normalize=false)
    @test norm(d) ≈ norm(c) atol = 1e-8
  end

  @testset "DMT rejects sites whose dimension is not a perfect square" begin
    sites = siteinds("S=1", 3)   # dim 3 is not d^2 for any integer d
    psi = MPS(sites, n -> "Up")
    gate = Matrix{ComplexF64}(I, 3, 3)
    @test_throws ArgumentError dmt_step!(psi, gate, 1; maxdim=12)
  end

  @testset "identity DMT step preserves a product operator" begin
    _, psi = _dmt_test_state(4)
    reference = copy(psi)

    dmt_step!(psi, _identity_gate(2), 2; maxdim=12, cutoff=1e-12, direction=:R, gate_maxdim=64)

    @test inner(reference, psi) ≈ 1.0 atol = 1e-10
    @test dmt_step!(psi, _identity_gate(2), 2; maxdim=12, cutoff=1e-12, direction=:R, gate_maxdim=64) === psi
  end

  @testset "DMT truncates an enlarged bond" begin
    sites = pauli_siteinds(5)
    psi = random_mps(ComplexF64, sites; linkdims=16)
    normalize!(psi)

    @test dim(linkind(psi, 2)) > 9
    dmt_step!(psi, _identity_gate(2), 2; maxdim=9, cutoff=1e-12, direction=:R, gate_maxdim=32)

    @test dim(linkind(psi, 2)) <= 9
    @test isfinite(real(inner(psi, psi)))
  end

  @testset "connector_buffer raises a migration error" begin
    @test_throws ArgumentError DMTOptions(maxdim=30, connector_buffer=8)
    @test_throws ArgumentError DMTGateEvolution(_identity_gate(2), 0.1; schedule=[1],
      maxdim=30, connector_buffer=8)

    _, psi = _dmt_test_state(3)
    @test_throws ArgumentError dmt_step!(psi, _identity_gate(2), 1; maxdim=30, connector_buffer=8)
    @test_throws ArgumentError dmt_step!(psi, _identity_gate(2), 1; maxdim=0)
  end

  @testset "maxdim is the total budget" begin
    sites = operator_siteinds(6; d=2)
    psi = random_mps(ComplexF64, sites; linkdims=40)
    normalize!(psi)
    MPSToolkit._dmt_bond_truncate!(psi, 3; maxdim=20, cutoff=1e-14)
    @test dim(linkind(psi, 3)) <= 20
    # 2 d^2 = 8 of those 20 are the protected block
    @test dim(linkind(psi, 3)) >= 8
  end

  @testset "maxdim is enforced on a genuinely traceless operator" begin
    # The tight case against the budget. With zero trace the connector is disabled (a = 0), so
    # B = S: `B QR` keeps its first column and `QL' B` its first row, both protected blocks are
    # full rank, and the reinstated matrix saturates 2 d^2 + (maxdim - 2 d^2) = maxdim exactly.
    # This is why the budget cannot be widened to `maxdim - 2 d^2 + 1` -- that would overflow
    # here and `_dmt_refactor` would clip a genuinely protected direction.
    sites = pauli_siteinds(7)
    noise = random_mps(ComplexF64, sites; linkdims=24)
    identity_state = pauli_basis_state(sites, fill(1, 7))
    # Subtracting the identity component zeroes the trace exactly while leaving the left and
    # right identity environments individually nonzero, which is what makes q0' S r0 vanish.
    traceless = add(noise, -inner(identity_state, noise) * identity_state; maxdim=32, cutoff=0.0)
    @test abs(pauli_trace(traceless)) < 1e-10 * norm(traceless)
    @test dim(linkind(traceless, 4)) > 20

    MPSToolkit._dmt_bond_truncate!(traceless, 4; maxdim=20, cutoff=0.0)
    @test dim(linkind(traceless, 4)) <= 20
    @test abs(pauli_trace(traceless)) < 1e-10 * norm(traceless)
    @test all(isfinite, [real(inner(traceless, traceless)), imag(inner(traceless, traceless))])
  end

  @testset "DMT truncates a generic bond into the total budget" begin
    sites = pauli_siteinds(6)
    psi = random_mps(sites; linkdims=20)
    normalize!(psi)

    @test dim(linkind(psi, 3)) == 20
    dmt_step!(psi, _identity_gate(2), 3; maxdim=12, cutoff=1e-12, gate_maxdim=40)

    # On a TRACE-CARRYING bond the reinstated matrix has RANK maxdim - 1, not maxdim.
    # `B = S - C` annihilates the identity directions, so `QL' B` has a zero first ROW and
    # `B QR` a zero first COLUMN: the rank-one trace connector shares one direction with the
    # two protected blocks instead of adding one, giving
    # 1 + (d^2 - 1) + (d^2 - 1) + (maxdim - 2 d^2) = maxdim - 1.
    # The reported DIMENSION only reaches that rank when `cutoff > 0`; at `cutoff = 0.0` the
    # surplus direction survives at ~1e-16 and the bond stays at maxdim. Either way the budget
    # chi' = maxdim - 2 d^2 never overflows maxdim -- and on a traceless bond, where the
    # connector is disabled, the rank does saturate maxdim (see the traceless testset above).
    @test dim(linkind(psi, 3)) == 11
  end

  @testset "DMT preserves identity and local Pauli data under truncation" begin
    sites = pauli_siteinds(5)
    terms = [
      pauli_basis_state(sites, fill(1, 5); coefficient=1.0),
      pauli_basis_state(sites, [1, 4, 1, 1, 1]; coefficient=0.3),
      pauli_basis_state(sites, fill(2, 5); coefficient=0.2),
    ]
    # A random component pushes the internal bonds past the maxdim below, so the truncation
    # actually fires; without it maxdim >= the new floor would leave the state untouched and
    # the preservation assertion would be vacuous.
    psi = add(terms..., 0.15 * random_mps(ComplexF64, sites; linkdims=16); maxdim=20, cutoff=1e-14)
    probes = [fill(1, 5), [1, 4, 1, 1, 1], [1, 1, 4, 1, 1], [1, 1, 1, 4, 1]]
    before = [inner(pauli_basis_state(sites, labels), psi) for labels in probes]

    @test maximum(dim(linkind(psi, b)) for b in 2:3) > 9
    dmt_step!(psi, _identity_gate(3), 2; maxdim=9, cutoff=1e-12, direction=:R, gate_maxdim=32)

    after = [inner(pauli_basis_state(sites, labels), psi) for labels in probes]
    @test after ≈ before atol = 1e-12
  end

  @testset "DMT preserves total S^z for a structured traceless operator" begin
    # Total S^z commutes with every XXZ bond gate, so for the traceless domain-wall operator
    # D = sum_j sign(j - kink) sigma^z_j the charge sum_x tr(D(t) sigma^z_x) is conserved at 0
    # exactly. Each summand is a diameter-1 observable, so the faithful kernel conserves it to
    # machine precision rather than merely to the truncation floor -- the 1e-2 tolerance this
    # testset used to carry is exactly what hid the clipped-protected-block defect. A structured
    # (evolved-Pauli) operator is the discriminating input; random unit-test inputs are not.
    nsites = 12
    sites = pauli_siteinds(nsites)
    state = pauli_domain_wall_state(sites; kink=nsites ÷ 2)
    gate = pauli_gate_from_hamiltonian(spinhalf_xyz_bond_hamiltonian(; Jx=1.0, Jy=1.0, Jz=1.0), 0.1)
    schedule = collect(1:(nsites - 1))
    evo = DMTGateEvolution(gate, 0.1; schedule=schedule, reverse_schedule=reverse(schedule),
      maxdim=24, cutoff=1e-10, gate_maxdim=96, normalize=false)
    Z = pauli_matrices().Z
    total_charge(psi) = sum(real.(pauli_expectation_profile(psi, [(x, Z) for x in 1:nsites]; normalize=false)))
    @test abs(total_charge(state)) < 1e-10           # exact at t = 0
    for _ in 1:8
      evolve!(state, evo)
      @test abs(total_charge(state)) < 1e-10         # conserved EXACTLY by the faithful kernel
    end
  end

  @testset "DMT truncates every internal bond of wider update windows" begin
    for span in (1, 2, 3, 4, 5)
      nsites = span + 3
      sites = pauli_siteinds(nsites)
      psi = random_mps(sites; linkdims=16)
      normalize!(psi)
      start = 2

      dmt_step!(psi, _identity_gate(span), start; maxdim=10, cutoff=1e-12, direction=:R, gate_maxdim=32)

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

    # An under-floor maxdim must be rejected BEFORE the gate is applied. The kernel-side check
    # alone is too late: `dmt_step!` reaches it only after `tebd_evolve!` has mutated and
    # re-gauged the state, which used to leave link dims [16,16,16,16,4] -> [4,16,16,16,4]
    # behind on a call that threw.
    wide = random_mps(ComplexF64, sites; linkdims=16)
    normalize!(wide)
    wide_reference = copy(wide)
    wide_dims = _link_dims(wide)
    @test_throws ArgumentError dmt_step!(wide, _identity_gate(2), 2; maxdim=5, gate_maxdim=64)
    @test _link_dims(wide) == wide_dims
    @test inner(wide_reference, wide) ≈ inner(wide_reference, wide_reference) atol = 1e-12
    # ... including for a span-1 gate, which never reaches the kernel at all.
    @test_throws ArgumentError dmt_step!(wide, _identity_gate(1), 2; maxdim=5, gate_maxdim=64)
    @test _link_dims(wide) == wide_dims

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

    _manual_dmt_sweep!(manual, gates; maxdim=12, cutoff=1e-12, gate_maxdim=64)
    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=12, cutoff=1e-12, gate_maxdim=64)
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

    _manual_dmt_sweep!(manual, gates; maxdim=12, cutoff=1e-12, gate_maxdim=64)
    evo = DMTGateEvolution(gates, 0.05; schedule=[1, 2, 3, 4], reverse_schedule=[4, 3, 2, 1], maxdim=12, cutoff=1e-12, gate_maxdim=64)
    dmt_evolve!(scheduled, evo)

    @test inner(manual, scheduled) ≈ 1.0 atol = 1e-8
  end

  @testset "threaded env cache reproduces the rebuild bit-for-bit (PXP hard case, verify on)" begin
    # Guards the perf-1 environment cache. With _DMT_VERIFY_ENVS on, _dmt_bond_truncate! throws
    # if any threaded environment differs from the from-scratch rebuild (index identity + norm of
    # difference), so reaching the end of dmt_evolve! proves the cached path is byte-faithful. The
    # PXP energy-transport schedule is the adversarial input: non-monotonic (bond 1 revisited),
    # mixed-span (2 and 3), overlapping multi-bond windows, run forward and reverse over 2 sweeps
    # on a generic high-operator-entanglement state (random_mps) where truncation actually fires.
    old = MPSToolkit._DMT_VERIFY_ENVS[]
    MPSToolkit._DMT_VERIFY_ENVS[] = true
    try
      nsites = 8
      psites = pauli_siteinds(nsites)
      terms = [(first(pxp_term_support(nsites, j)), Matrix(pxp_term_hamiltonian(nsites, j))) for j in 1:nsites]
      gates = [pauli_gate_from_hamiltonian(h, 0.05) for (_, h) in terms]
      schedule = [start for (start, _) in terms]
      evo = DMTGateEvolution(gates, 0.05; schedule=schedule, reverse_schedule=reverse(schedule),
        nstep=2, maxdim=12, cutoff=1e-12, gate_maxdim=64)

      psi = random_mps(ComplexF64, psites; linkdims=24)
      normalize!(psi)
      @test dmt_evolve!(psi, evo) === psi   # no env-verify assertion fired

      # And explicitly: the threaded result equals the rebuild path (cache=nothing) — same state
      # AND same trimmed bond dimensions (the canonical-form output, not merely the physical state).
      MPSToolkit._DMT_VERIFY_ENVS[] = false
      base = random_mps(ComplexF64, psites; linkdims=24)
      normalize!(base)
      threaded = copy(base)
      dmt_evolve!(threaded, evo)
      rebuilt = copy(base)
      for _ in 1:evo.nstep
        for (i, b) in pairs(evo.schedule)
          dmt_step!(rebuilt, MPSToolkit._gate_for_step(evo.gate, b, i), b; maxdim=evo.maxdim,
            cutoff=evo.cutoff, direction=:R, gate_maxdim=evo.gate_maxdim,
            preserve_diameter=evo.preserve_diameter, truncation=evo.truncation)
        end
        for (i, b) in pairs(evo.reverse_schedule)
          g = MPSToolkit._reverse_gate_for_step(evo.gate, true, evo.schedule, b, i)
          dmt_step!(rebuilt, g, b; maxdim=evo.maxdim, cutoff=evo.cutoff, direction=:L,
            gate_maxdim=evo.gate_maxdim, preserve_diameter=evo.preserve_diameter,
            truncation=evo.truncation)
        end
      end
      normalize!(rebuilt)
      @test abs(inner(threaded, rebuilt)) ≈ 1.0 atol = 1e-10
      @test _link_dims(threaded) == _link_dims(rebuilt)
    finally
      MPSToolkit._DMT_VERIFY_ENVS[] = old
    end
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

    _manual_dmt_sweep!(manual, gates; maxdim=12, cutoff=1e-12, gate_maxdim=16)
    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 3], reverse_schedule=[3, 2, 1], maxdim=12, cutoff=1e-12, gate_maxdim=16)
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
      dmt_step!(manual, gate, bond; maxdim=12, cutoff=1e-12, direction=:R, gate_maxdim=16)
    end
    for (bond, gate) in zip([1, 2, 1], reverse(gates))
      dmt_step!(manual, gate, bond; maxdim=12, cutoff=1e-12, direction=:L, gate_maxdim=16)
    end
    normalize!(manual)

    evo = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 1], reverse_schedule=[1, 2, 1], maxdim=12, cutoff=1e-12, gate_maxdim=16)
    dmt_evolve!(scheduled, evo)

    @test abs(inner(manual, scheduled)) ≈ 1.0 atol = 1e-10
    ambiguous = DMTGateEvolution(gates, 0.1; schedule=[1, 2, 1], reverse_schedule=[2, 1, 1], maxdim=12, cutoff=1e-12, gate_maxdim=16)
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

    _manual_dmt_sweep!(manual, gates; maxdim=32, cutoff=1e-12, gate_maxdim=64)
    evo = DMTGateEvolution((bond, index) -> gates[index], 0.1; schedule=[1, 2, 3], reverse_schedule=[3, 2, 1], maxdim=32, cutoff=1e-12, gate_maxdim=64)
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

    evo = DMTGateEvolution((bond, index) -> gate, 0.1; schedule=forward_schedule, reverse_schedule=reverse_schedule, maxdim=12, cutoff=1e-12, gate_maxdim=16)
    @test_throws ArgumentError dmt_evolve!(scheduled, evo)
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
    function _single_bond_sequence!(psi, start, span, direction; maxdim, cutoff)
      bonds = collect(start:(start + span - 2))
      direction === :L && (bonds = reverse(bonds))
      for bond in bonds
        MPSToolkit._dmt_window_truncate!(psi, bond, 2; maxdim=maxdim, cutoff=cutoff,
          direction=direction)
      end
      return psi
    end

    for direction in (:R, :L), span in (3, 4, 5)
      nsites = span + 3
      sites = pauli_siteinds(nsites)
      psi = random_mps(ComplexF64, sites; linkdims=24)
      normalize!(psi)

      windowed = copy(psi)
      MPSToolkit._dmt_window_truncate!(windowed, 2, span; maxdim=12, cutoff=1e-12,
        direction=direction)

      sequential = copy(psi)
      _single_bond_sequence!(sequential, 2, span, direction; maxdim=12, cutoff=1e-12)

      @test abs(inner(windowed, sequential)) / (norm(windowed) * norm(sequential)) ≈ 1.0 atol = 1e-10
    end
  end
end
