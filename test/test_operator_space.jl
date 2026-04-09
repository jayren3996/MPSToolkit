using ITensors
using ITensorMPS
using LinearAlgebra

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

    total_sz = pauli_total_sz_state(sites)
    expected = sum(0.5 * _basis_product_mps(sites, [j == n ? 4 : 1 for j in 1:length(sites)]) for n in 1:length(sites))
    @test inner(expected, total_sz) ≈ length(sites) / 4
    @test inner(total_sz, total_sz) ≈ length(sites) / 4
  end

  @testset "Pauli-basis gate helpers" begin
    hbond = kron(ComplexF64[1 0; 0 -1], ComplexF64[1 0; 0 -1])
    gate = pauli_gate_from_hamiltonian(hbond, 0.0)
    state = pauli_basis_state(pauli_siteinds(2), [2, 3])

    evolved = copy(state)
    tebd_evolve!(evolved, gate, 1; maxdim=16, cutoff=0.0)

    @test inner(state, evolved) ≈ 1.0
    @test size(gate) == (16, 16)
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
    @test gate ≈ expected_gate atol = 1e-10
    @test pauli_gate_from_lindbladian(zeros(ComplexF64, 2, 2), ComplexF64[0 0; 0 0], dt) ≈ Matrix{ComplexF64}(I, 4, 4) atol = 1e-10

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
end
