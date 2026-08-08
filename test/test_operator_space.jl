using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using SparseArrays
using Test

function _basis_product_mps(sites, labels::Vector{Int})
  tensors = ITensor[]
  for (site, label) in zip(sites, labels)
    tensor = ITensor(site)
    tensor[site => label] = 1.0
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

function _dense_operator_from_basis(labels::Vector{Int})
  local_ops = (
    ComplexF64[1 0; 0 1],
    ComplexF64[0 1; 1 0],
    ComplexF64[0 -im; im 0],
    ComplexF64[1 0; 0 -1],
  )
  dense = local_ops[labels[1]]
  for label in labels[2:end]
    dense = kron(dense, local_ops[label])
  end
  return dense
end

@testset "operator-space projector MPOs" begin
  sites = pauli_siteinds(4)

  @testset "Pauli-basis state helpers" begin
    basis_state = pauli_basis_state(sites, [4, 1, 1, 1])
    @test inner(basis_state, _basis_product_mps(sites, [4, 1, 1, 1])) ≈ 1.0

    mixed_labels = pauli_basis_state(sites, [:I, " z ", 2, "x"])
    @test inner(mixed_labels, _basis_product_mps(sites, [1, 4, 2, 2])) ≈ 1.0
    @test_throws ArgumentError pauli_siteinds(0)
    @test_throws ArgumentError pauli_basis_state(sites, [:I])
    @test_throws ArgumentError pauli_basis_state(sites, [:Q, :I, :I, :I])

    total_sz = pauli_total_sz_state(sites)
    expected_coeff = 2.0^(length(sites) / 2 - 1)
    expected = sum(expected_coeff * _basis_product_mps(sites, [j == n ? 4 : 1 for j in 1:length(sites)]) for n in 1:length(sites))
    @test inner(expected, total_sz) ≈ length(sites) * expected_coeff^2
    @test inner(total_sz, total_sz) ≈ length(sites) * expected_coeff^2

    overridden = pauli_total_sz_state(sites; coefficient=0.5)
    @test inner(pauli_basis_state(sites, [4, 1, 1, 1]), overridden) ≈ 0.5
  end

  @testset "pauli_total_sz_state uses normalized Pauli-basis coefficients" begin
    for nsites in (1, 2, 4)
      local_sites = pauli_siteinds(nsites)
      total_sz = pauli_total_sz_state(local_sites)
      expected_coeff = 2.0^(nsites / 2 - 1)

      for target in 1:nsites
        labels = fill(1, nsites)
        labels[target] = 4
        probe = pauli_basis_state(local_sites, labels)
        @test inner(probe, total_sz) ≈ expected_coeff atol = 1e-12
      end
    end
  end

  @testset "pauli_domain_wall_state encodes a signed single-Z sum" begin
    for nsites in (2, 4, 5)
      local_sites = pauli_siteinds(nsites)
      kink = nsites ÷ 2
      dw = pauli_domain_wall_state(local_sites)
      # The per-site overlap with the single-Z probe is the signed coefficient
      # c_j = (j <= kink ? -1 : +1) -- i.e. tr(σ^z_j · D) up to the basis normalization.
      for target in 1:nsites
        labels = fill(1, nsites)
        labels[target] = 4
        probe = pauli_basis_state(local_sites, labels)
        expected = target <= kink ? -1.0 : 1.0
        @test inner(probe, dw) ≈ expected atol = 1e-12
      end
      # Compact bond-dimension-2 encoding of the operator sum.
      @test maxlinkdim(dw) == 2
      # Traceless: there is no identity-string component.
      @test inner(pauli_basis_state(local_sites, fill(1, nsites)), dw) ≈ 0.0 atol = 1e-12
    end

    # Custom kink position and coefficient.
    s = pauli_siteinds(6)
    dw = pauli_domain_wall_state(s; kink=2, coefficient=0.5)
    for target in 1:6
      labels = fill(1, 6)
      labels[target] = 4
      @test inner(pauli_basis_state(s, labels), dw) ≈ (target <= 2 ? -0.5 : 0.5) atol = 1e-12
    end
  end

  @testset "Pauli-basis gate helpers" begin
    hbond = kron(ComplexF64[1 0; 0 -1], ComplexF64[1 0; 0 -1])
    gate = pauli_gate_from_hamiltonian(hbond, 0.0)
    state = pauli_basis_state(pauli_siteinds(2), [2, 3])

    evolved = copy(state)
    tebd_evolve!(evolved, gate, 1; maxdim=16, cutoff=0.0)

    @test inner(state, evolved) ≈ 1.0
    @test size(gate) == (16, 16)

    swap = ComplexF64[1 0 0 0;
                      0 0 1 0;
                      0 1 0 0;
                      0 0 0 1]
    swap_gate = pauli_gate(swap)
    swap_sites = pauli_siteinds(2)
    swapped = pauli_basis_state(swap_sites, [2, 4])
    tebd_evolve!(swapped, swap_gate, 1; maxdim=16, cutoff=0.0)
    @test inner(pauli_basis_state(swap_sites, [4, 2]), swapped) ≈ 1.0 atol = 1e-10

    rotation = exp(-0.37im * ComplexF64[0 1; 1 0])
    rotation_gate = pauli_gate(rotation)
    @test rotation_gate' * rotation_gate ≈ Matrix{ComplexF64}(I, 4, 4) atol = 1e-10
    @test_throws ArgumentError pauli_gate(zeros(ComplexF64, 2, 3))
    @test_throws ArgumentError pauli_gate(ones(ComplexF64, 1, 1))
    @test_throws ArgumentError pauli_gate(zeros(ComplexF64, 3, 3))
  end

  @testset "pauli_gate_from_hamiltonian multi-site correctness" begin
    σx = ComplexF64[0 1; 1 0]
    σy = ComplexF64[0 -im; im 0]
    σz = ComplexF64[1 0; 0 -1]
    id2 = Matrix{ComplexF64}(I, 2, 2)
    zero_jumps = Matrix{ComplexF64}[]

    @testset "2-site: super-operator is real for Hermitian H" begin
      # Transverse-field Ising-like Hamiltonian on two sites
      h = kron(σx, id2) + kron(id2, σx) + 0.5 * kron(σz, σz)
      dt = 0.3
      gate = pauli_gate_from_hamiltonian(h, dt)
      @test maximum(abs.(imag.(gate))) < 1e-10
    end

    @testset "2-site: matches pauli_gate_from_lindbladian with no jumps" begin
      # For any Hermitian H, pauli_gate_from_hamiltonian(H, dt) must equal
      # pauli_gate_from_lindbladian(H, [], dt) because the latter computes
      # exp(dt * -i[H, ·]), which is exactly the coherent super-operator.
      h = kron(σx, id2) + kron(id2, σx) + 0.5 * kron(σz, σz)
      dt = 0.3
      gate_ham = pauli_gate_from_hamiltonian(h, dt)
      gate_lin = pauli_gate_from_lindbladian(h, zero_jumps, dt)
      @test gate_ham ≈ gate_lin atol = 1e-10
    end

    @testset "2-site: XY bond super-operator matches lindbladian path" begin
      # Non-symmetric coupling that was broken in the original implementation
      # due to the kron(vec(P₁), vec(P₂)) vs vec(P₁ ⊗ P₂) basis mismatch.
      h = kron(σx, σx) + kron(σy, σy)
      dt = 0.4
      gate_ham = pauli_gate_from_hamiltonian(h, dt)
      gate_lin = pauli_gate_from_lindbladian(h, zero_jumps, dt)
      @test gate_ham ≈ gate_lin atol = 1e-10
    end

    @testset "3-site: super-operator matches lindbladian path" begin
      h = kron(σx, id2, id2) + kron(id2, σz, σz) + 0.2 * kron(σy, σy, id2)
      dt = 0.15
      gate_ham = pauli_gate_from_hamiltonian(h, dt)
      gate_lin = pauli_gate_from_lindbladian(h, Matrix{ComplexF64}[], dt)
      @test gate_ham ≈ gate_lin atol = 1e-10
    end

    @testset "tebd_evolve! with coherent 2-site gate preserves Hermiticity" begin
      # Apply a non-trivial coherent gate to a Pauli-basis MPS built from a
      # Hermitian density matrix, and confirm that reading back the Pauli
      # coefficients gives real values (i.e. the gate preserves Hermiticity).
      #
      # Use h = σx ⊗ σz (NOT swap-symmetric) and start from ρ = |↑↑⟩⟨↑↑|. The
      # gate `exp(-im*dt*h)` rotates |↑↑⟩ → cos(dt)|↑↑⟩ - im sin(dt)|↓↑⟩, so
      # the evolved density matrix has non-trivial mixing and the gate is not
      # a no-op on this state.
      h = kron(σx, σz)
      dt = 0.5
      gate = pauli_gate_from_hamiltonian(h, dt)

      # ρ = |↑↑⟩⟨↑↑| = (I+Z)/2 ⊗ (I+Z)/2. In the normalized Pauli basis
      # P_α = σ_α/√2, the nonzero coefficients are (I,I),(I,Z),(Z,I),(Z,Z).
      pauli_sites2 = pauli_siteinds(2)
      uu_labels = [(1,1), (1,4), (4,1), (4,4)]
      ρ_mps = sum(0.5 * pauli_basis_state(pauli_sites2, [a, b]) for (a, b) in uu_labels)
      evolved = copy(ρ_mps)
      tebd_evolve!(evolved, gate, 1; maxdim=16, cutoff=0.0)

      # Reading each Pauli coefficient via inner product with the corresponding
      # basis state should yield a real number (Hermiticity of ρ).
      max_imag = 0.0
      for a in 1:4, b in 1:4
        coeff = inner(pauli_basis_state(pauli_sites2, [a, b]), evolved)
        max_imag = max(max_imag, abs(imag(coeff)))
      end
      @test max_imag < 1e-10
    end
  end

  @testset "Pauli-basis Lindblad helpers" begin
    gamma = 0.3
    dt = 0.5
    jump = sqrt(gamma) .* ComplexF64[1 0; 0 -1]
    generator = pauli_lindblad_generator(zeros(ComplexF64, 2, 2), [jump])
    gate = pauli_gate_from_lindbladian(zeros(ComplexF64, 2, 2), [jump], dt)

    expected_generator = Diagonal(ComplexF64[0.0, -2gamma, -2gamma, 0.0])
    expected_gate = Diagonal(exp.(dt .* diag(expected_generator)))

    @test generator ≈ expected_generator atol = 1e-10
    @test generator[1, :] ≈ zeros(ComplexF64, 4) atol = 1e-10
    @test gate ≈ expected_gate atol = 1e-10
    @test pauli_gate_from_lindbladian(zeros(ComplexF64, 2, 2), ComplexF64[0 0; 0 0], dt) ≈ Matrix{ComplexF64}(I, 4, 4) atol = 1e-10
    @test_throws ArgumentError pauli_lindblad_generator(ones(ComplexF64, 1, 1), Matrix{ComplexF64}[])
    @test_throws ArgumentError pauli_lindblad_generator(zeros(ComplexF64, 2, 2), [zeros(ComplexF64, 3, 3)])

    sites1 = pauli_siteinds(1)
    identity_density = pauli_basis_state(sites1, ["I"]; coefficient=2.0^(-0.5))
    identity_operator = pauli_basis_state(sites1, ["I"])
    damping_gate = pauli_gate_from_lindbladian(zeros(ComplexF64, 2, 2), [ComplexF64[0 0; 1 0]], 0.2)
    evolved = copy(identity_density)
    tebd_evolve!(evolved, damping_gate, 1; maxdim=4, cutoff=0.0)
    @test abs(sqrt(2.0) * inner(identity_operator, evolved) - 1.0) ≤ 1e-10
  end

  @testset "DAOE projector keeps short strings and damps long strings" begin
    projector = pauli_daoe_projector(sites; lstar=1, gamma=0.7)
    short = _basis_product_mps(sites, [2, 1, 1, 1])
    long = _basis_product_mps(sites, [2, 3, 4, 1])

    projected_short = apply(projector, short; maxdim=16, cutoff=0.0)
    projected_long = apply(projector, long; maxdim=16, cutoff=0.0)

    @test abs(inner(short, projected_short) - 1.0) ≤ 1e-10
    @test abs(inner(long, projected_long) - exp(-0.7 * 2)) ≤ 1e-10
  end

  @testset "FDAOE projector respects fermionic weights in the Pauli basis" begin
    projector = fdaoe_projector(sites; wstar=2, gamma=0.4)
    single_z = _basis_product_mps(sites, [4, 1, 1, 1])
    adjacent_xx = _basis_product_mps(sites, [2, 2, 1, 1])
    double_z = _basis_product_mps(sites, [4, 4, 1, 1])
    jw_tail = _basis_product_mps(sites, [2, 4, 4, 2])

    projected_single_z = apply(projector, single_z; maxdim=16, cutoff=0.0)
    projected_adjacent_xx = apply(projector, adjacent_xx; maxdim=16, cutoff=0.0)
    projected_double_z = apply(projector, double_z; maxdim=16, cutoff=0.0)
    projected_jw_tail = apply(projector, jw_tail; maxdim=16, cutoff=0.0)

    @test abs(inner(single_z, projected_single_z) - 1.0) ≤ 1e-10
    @test abs(inner(adjacent_xx, projected_adjacent_xx) - 1.0) ≤ 1e-10
    @test abs(inner(double_z, projected_double_z) - exp(-0.4 * 2)) ≤ 1e-10
    @test abs(inner(jw_tail, projected_jw_tail) - 1.0) ≤ 1e-10
  end

  @testset "DAOE damping weight matches the exp(-gamma*max(size-lstar,0)) oracle" begin
    # Operator "size" in DAOE is the Hamming weight -- the number of non-identity Pauli sites,
    # independent of their span (a gapped weight-k string damps identically to a contiguous one).
    # Sweep size and lstar against the closed form on a 5-site chain so every size is reachable,
    # including gapped/edge-separated strings the 4-string smoke tests above never exercise.
    oracle_sites = pauli_siteinds(5)
    daoe_weight(labels, lstar, gamma) = exp(-gamma * max(count(!=(1), labels) - lstar, 0))
    test_strings = [
      [1, 1, 1, 1, 1],   # identity, size 0
      [4, 1, 1, 1, 1],   # size 1
      [2, 1, 3, 1, 1],   # size 2, gapped
      [2, 1, 3, 1, 4],   # size 3, edge-separated
      [2, 3, 4, 2, 3],   # size 5, full
    ]
    for lstar in (0, 1, 2), gamma in (0.3, 0.9)
      projector = pauli_daoe_projector(oracle_sites; lstar=lstar, gamma=gamma)
      for labels in test_strings
        state = _basis_product_mps(oracle_sites, labels)
        projected = apply(projector, state; maxdim=16, cutoff=0.0)
        @test inner(state, projected) ≈ daoe_weight(labels, lstar, gamma) atol = 1e-10
      end
    end
  end

  @testset "DAOE and FDAOE damping preserve the operator trace" begin
    # The weight-0 identity component carries tr(O); DAOE/FDAOE damp only the higher-weight
    # tail, leaving the identity coefficient -- and hence the trace -- untouched. Regression for
    # the corrected daoe.md guidance: these are size-damping channels, not trace-changing maps.
    trace_sites = pauli_siteinds(4)
    operator = 3.0 * _basis_product_mps(trace_sites, [1, 1, 1, 1])
    for (labels, coeff) in (([2, 2, 1, 1], 1.5), ([4, 3, 2, 1], 0.8), ([2, 3, 4, 2], 0.4))
      operator = add(operator, coeff * _basis_product_mps(trace_sites, labels); cutoff=0.0)
    end
    trace_before = pauli_trace(operator)
    for projector in (pauli_daoe_projector(trace_sites; lstar=1, gamma=0.7),
                      pauli_fdaoe_projector(trace_sites; wstar=2, gamma=0.4))
      projected = apply(projector, operator; maxdim=32, cutoff=0.0)
      @test pauli_trace(projected) ≈ trace_before atol = 1e-10
    end
  end

  @testset "DAOE/FDAOE projectors validate cutoffs and gamma" begin
    @test_throws ArgumentError pauli_daoe_projector(sites; lstar=-1, gamma=0.1)
    @test_throws ArgumentError pauli_daoe_projector(sites; lstar=2, gamma=-0.1)
    @test_throws ArgumentError pauli_fdaoe_projector(sites; wstar=-1, gamma=0.1)
    @test_throws ArgumentError pauli_fdaoe_projector(sites; wstar=2, gamma=-0.1)
    # The fdaoe_projector alias shares the implementation and its validation.
    @test_throws ArgumentError fdaoe_projector(sites; wstar=-1, gamma=0.1)
    @test pauli_fdaoe_projector(sites; wstar=2, gamma=0.4) isa MPO
  end
end

@testset "generic operator-space state builders" begin
  @testset "pauli_* wrappers are exactly the d = 2 case" begin
    @test dim(first(operator_siteinds(3; d = 2))) == 4
    @test dim(first(operator_siteinds(3; d = 3))) == 9

    sites = pauli_siteinds(4)
    a = pauli_basis_state(sites, ["I", "Z", "X", "I"])
    b = operator_basis_state(sites, [1, 4, 2, 1])
    @test inner(a, b) ≈ 1.0 atol = 1e-14

    # pauli_total_sz_state/pauli_domain_wall_state keep their pre-existing bond-dimension-2
    # amplitudes (identity_amplitude = 1), while operator_local_sum_state now builds the
    # *literal* operator sum (identity_amplitude = sqrt(d)) -- see the dedicated pin testset
    # below. Built from the same weights/coeffs, the two therefore differ by exactly
    # sqrt(d)^(nsites - 1), not by a factor of 1.
    scale = sqrt(2.0)^(4 - 1)
    total_via_wrapper = pauli_total_sz_state(sites)
    total_via_generic = operator_local_sum_state(
      sites, pauli_matrices().Z / sqrt(2), fill(2.0^(4 / 2 - 1), 4))
    @test inner(total_via_wrapper, total_via_generic) ≈ scale * inner(total_via_wrapper, total_via_wrapper) atol = 1e-12

    dw_wrapper = pauli_domain_wall_state(sites; kink = 2)
    dw_generic = operator_local_sum_state(
      sites, pauli_matrices().Z / sqrt(2), [j <= 2 ? -1.0 : 1.0 for j in 1:4])
    @test inner(dw_wrapper, dw_generic) ≈ scale * inner(dw_wrapper, dw_wrapper) atol = 1e-12
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
    # Coefficient of "S^z on site j, identity elsewhere" is coeffs[j] * tr(P_mu' * S^z) times
    # sqrt(d)^(nsites - 1): operator_local_sum_state now inserts a literal identity (amplitude
    # sqrt(d) on basis element 1, not 1) on every other site, so the literal operator it
    # represents is exactly sum_j coeffs[j] * S^z_j -- see the dedicated literal-sum testset
    # below for the direct check against operator_product_state.
    literal_scale = sqrt(d)^(length(sites) - 1)
    for j in 1:4, mu in 2:d^2
      labels = [k == j ? mu : 1 for k in 1:4]
      probe = operator_basis_state(sites, labels)
      @test inner(probe, state) ≈ literal_scale * coeffs[j] * tr(basis[mu]' * sz) atol = 1e-12
    end
  end

  @testset "operator_local_sum_state equals the literal operator sum" begin
    # Pin operator_local_sum_state against literal per-term product states built with
    # operator_product_state, for both a traceless and a non-traceless local operator, at two
    # different local dimensions. This is the exact composition later transport tasks rely on
    # (e.g. rho = I + sum_j c_j S^z_j), so the two builders must agree on what "identity
    # elsewhere" means.
    cases = (
      (d = 2, op = ComplexF64[2 0; 0 1], nsites = 4),      # non-traceless, d = 2
      (d = 2, op = ComplexF64[1 0; 0 -1], nsites = 4),     # traceless (sigma^z), d = 2
      (d = 3, op = ComplexF64.(diagm([2, 1, 0])), nsites = 3),   # non-traceless, d = 3
      (d = 3, op = ComplexF64.(diagm([1, 0, -1])), nsites = 3),  # traceless (spin-1 S^z), d = 3
    )
    for case in cases
      d, op, nsites = case.d, case.op, case.nsites
      sites = operator_siteinds(nsites; d = d)
      identity_matrix = Matrix{ComplexF64}(I, d, d)
      coeffs = [1.0, -2.0, 3.0, -4.0][1:nsites]

      single(j) = operator_product_state(sites, [k == j ? op : identity_matrix for k in 1:nsites];
                                          coefficient = coeffs[j])
      reference = single(1)
      for j in 2:nsites
        reference = add(reference, single(j); maxdim = 16, cutoff = 0.0)
      end
      built = operator_local_sum_state(sites, op, coeffs)
      @test inner(reference, built) ≈ inner(reference, reference) atol = 1e-10
      @test norm(built) ≈ norm(reference) atol = 1e-10
    end
  end

  @testset "pauli_total_sz_state/pauli_domain_wall_state pin identity_amplitude = 1" begin
    # These wrappers deliberately keep their pre-refactor bond-dimension-2 numerics (unlike
    # operator_local_sum_state, which now builds the literal operator sum). Pin that convention
    # directly: built from the *same* normalized op (pauli_matrices().Z / sqrt(2)), the wrapper
    # and operator_local_sum_state must differ by exactly identity_amplitude^(nsites - 1), i.e.
    # 1^(nsites-1) vs sqrt(2)^(nsites-1) -- nothing else.
    for nsites in (1, 4)
      sites = pauli_siteinds(nsites)
      kink = min(2, nsites)
      coefficient = 0.75
      scale = 2.0^((nsites - 1) / 2)

      dw = pauli_domain_wall_state(sites; kink = kink, coefficient = coefficient)
      dw_weights = [j <= kink ? -coefficient : coefficient for j in 1:nsites]
      dw_reference = operator_local_sum_state(sites, pauli_matrices().Z / sqrt(2), dw_weights)
      @test inner(dw, dw_reference) ≈ scale * inner(dw, dw) atol = 1e-12
      @test norm(dw_reference) ≈ scale * norm(dw) atol = 1e-12

      total = pauli_total_sz_state(sites)
      total_weights = fill(2.0^(nsites / 2 - 1), nsites)
      total_reference = operator_local_sum_state(sites, pauli_matrices().Z / sqrt(2), total_weights)
      @test inner(total, total_reference) ≈ scale * inner(total, total) atol = 1e-12
      @test norm(total_reference) ≈ scale * norm(total) atol = 1e-12

      # nsites == 1 has no "unoccupied" site, so identity_amplitude never enters: scale == 1
      # and the wrapper and operator_local_sum_state must coincide exactly.
      if nsites == 1
        @test scale ≈ 1.0 atol = 1e-14
      end
    end
  end

  @testset "pauli_total_sz_state/pauli_domain_wall_state reject non-spin-1/2 sites" begin
    # Regression: these wrappers used to route through operator_local_sum_state, which derives
    # d from `sites` and rejects a size mismatch against the fixed 2x2 Z/sqrt(2) operator. After
    # the identity_amplitude fix they call _local_sum_state directly and must therefore validate
    # d == 2 themselves -- otherwise a d = 3 site silently builds a bogus bond-dimension-2 MPS
    # (every valid site dimension d^2 >= 4 keeps the tensor[site => mu] writes for mu in 1:4 in
    # bounds, so no natural BoundsError would catch this).
    @test_throws ArgumentError pauli_total_sz_state(operator_siteinds(4; d = 3))
    @test_throws ArgumentError pauli_domain_wall_state(operator_siteinds(4; d = 3))
  end

  @testset "label validation is dimension aware" begin
    @test_throws ArgumentError operator_basis_state(operator_siteinds(2; d = 3), [10, 1])
    # Pauli letter labels are only meaningful at d = 2.
    @test_throws ArgumentError operator_basis_state(operator_siteinds(2; d = 3), ["X", "I"])
  end
end

@testset "generic operator-space gates" begin
  @testset "pauli_gate matches the explicit d = 2 definition" begin
    # Ground truth independent of operator_gate's own implementation: recompute
    # G[a, b] = tr(P_a' * u * P_b * u') directly from operator_basis_matrices(2), the same
    # definition operator_gate's docstring states. pauli_gate(u) is literally
    # operator_gate(u; d=2) (see gates.jl), so comparing pauli_gate against operator_gate
    # would be tautological -- this checks the convention itself, independently.
    x = ComplexF64[0 1; 1 0]
    z = ComplexF64[1 0; 0 -1]
    dt = 0.3
    u = exp(-dt * im * kron(x, z))
    # Genuinely non-Hermitian: a real-time evolution gate is unitary, not Hermitian in
    # general. A Hermitian u would not distinguish the correct A*P*A' from a wrong
    # conjugation such as A*P*A, so pin this so a future edit can't silently make u Hermitian.
    @test norm(u - u') > 0
    basis = operator_basis_matrices(2)
    two_site = [kron(a, b) for a in basis for b in basis]
    expected = ComplexF64[
      tr(two_site[row]' * u * two_site[col] * u') for row in eachindex(two_site), col in eachindex(two_site)
    ]
    @test pauli_gate(u) ≈ expected atol = 1e-12
  end

  @testset "pauli_gate_from_hamiltonian delegates to operator_gate_from_hamiltonian" begin
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
    d = 3
    sz = sparse(ComplexF64[1 0 0; 0 0 0; 0 0 -1])
    h = kron(sz, sz)
    @test h isa AbstractSparseMatrix
    dense_gate = operator_gate_from_hamiltonian(Matrix(h), 0.1; d = d)
    @test operator_gate_from_hamiltonian(h, 0.1; d = d) ≈ dense_gate atol = 1e-12
  end

  @testset "operator_gate densifies sparse input directly" begin
    # The test above only exercises sparse input through operator_gate_from_hamiltonian,
    # whose exp(...) densifies before operator_gate ever sees it. Call operator_gate itself
    # on a sparse matrix to exercise its own Matrix{ComplexF64}(op) densification.
    d = 3
    sz = sparse(ComplexF64[1 0 0; 0 0 0; 0 0 -1])
    @test sz isa AbstractSparseMatrix
    @test operator_gate(sz; d = d) ≈ operator_gate(Matrix(sz); d = d) atol = 1e-12
  end

  @testset "span inference rejects incompatible sizes" begin
    # A dense physical matrix of size d^span x d^span acts on `span` sites of local dimension
    # d: 3^2 = 9 is two sites, 3^4 = 81 is four sites (matches the "d = 3 gate conjugates
    # correctly" testset above, which independently pins a 9x9 two-site Hamiltonian to a
    # d^4 x d^4 gate).
    @test MPSToolkit._operator_span(9, 3) == 2
    @test MPSToolkit._operator_span(81, 3) == 4
    @test_throws ArgumentError MPSToolkit._operator_span(8, 3)
  end
end

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

    # Independent ground truth via dense linear algebra, using the same _dense_operator_from_basis
    # helper the label-based tests at the top of this file use (not pauli_trace/operator_trace
    # themselves -- comparing those against each other is tautological, since pauli_trace(rho) IS
    # literally operator_trace(rho), a one-line delegation that cannot fail regardless of whether
    # the underlying math is right; this is the same defect class fixed in Task 3's F2 finding).
    # rho = pauli_basis_state(sites, fill(1, 4)) + pauli_domain_wall_state(sites; kink = 2) is,
    # in the normalized-Pauli convention, literal_scale * I_16 + literal_scale * sum_j coeffs[j] *
    # Z_j (embedded, identity elsewhere), with literal_scale = (1/sqrt(2))^4 from the four
    # P_1 = I/sqrt(2) (resp. P_4 = Z/sqrt(2)) factors, and coeffs = [-1, -1, 1, 1] for kink = 2,
    # coefficient = 1.0 (default): sites 1:kink get -1, the rest get +1.
    literal_scale = (1 / sqrt(2))^4
    coeffs = [-1.0, -1.0, 1.0, 1.0]
    dense_rho = literal_scale * _dense_operator_from_basis(fill(1, 4))
    for (j, c) in enumerate(coeffs)
      labels = fill(1, 4)
      labels[j] = 4
      dense_rho += literal_scale * c * _dense_operator_from_basis(labels)
    end
    expected_trace = tr(dense_rho)
    expected_profile = ComplexF64[
      tr(dense_rho * _dense_operator_from_basis([j == x ? 4 : 1 for j in 1:4])) for x in 1:4
    ]

    @test operator_trace(rho) ≈ expected_trace atol = 1e-10
    @test operator_expectation_profile(rho, terms; normalize = false) ≈ expected_profile atol = 1e-10

    # Honest delegation checks: pauli_trace/pauli_expectation_profile ARE operator_trace/
    # operator_expectation_profile here, so these only confirm the wrapper forwards its
    # arguments unchanged, not that the underlying math is right -- that is the pinned
    # comparison above. Exact equality, not approximate: it is the same call.
    @test pauli_trace(rho) == operator_trace(rho)
    @test pauli_expectation_profile(rho, terms; normalize = false) ==
          operator_expectation_profile(rho, terms; normalize = false)
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

@testset "spin-1/2-only helpers reject higher local dimension" begin
  # DAOE and the PXP constraint helpers key off the Pauli (I, X, Y, Z) label ordering and the
  # two-level PXP blockade construction; they are out of scope for the higher-spin
  # generalization and must say so rather than silently producing nonsense at d = 3.
  sites = operator_siteinds(4; d = 3)
  @test_throws ArgumentError pauli_daoe_projector(sites; lstar = 2, gamma = 0.5)
  @test_throws ArgumentError pauli_fdaoe_projector(sites; wstar = 2, gamma = 0.5)
  @test_throws ArgumentError pauli_pxp_constraint_state(sites)
  @test_throws ArgumentError pauli_pxp_constraint_projector(sites)
end
