"""
    _pxp_physical_sites(nsites)

Create the internal physical dimension-2 site indices used to build PXP constraint objects in
operator space.
"""
function _pxp_physical_sites(nsites::Integer)
  return [Index(2, "PXPSite,n=$(n)") for n in 1:Int(nsites)]
end

function _pauli_pxp_constraint_state_real(sites)
  nsites = length(sites)
  basis = operator_basis_matrices(2)
  links = [Index(2, "PXPConstraintStateLink,n=$(n)") for n in 1:(nsites - 1)]
  tensors = ITensor[]
  for n in 1:nsites
    tensor = nsites == 1 ? ITensor(Float64, sites[n]) :
      n == 1 ? ITensor(Float64, sites[n], links[n]) :
      n == nsites ? ITensor(Float64, links[n - 1], sites[n]) :
      ITensor(Float64, links[n - 1], sites[n], links[n])
    previous_states = n == 1 ? (1:1) : (1:2)
    for previous in previous_states, local_state in 1:2, alpha in eachindex(basis)
      next_state, allowed = _pxp_constraint_transition(previous, local_state)
      coefficient = allowed * real(basis[alpha][local_state, local_state])
      iszero(coefficient) && continue
      if nsites == 1
        tensor[sites[n] => alpha] += coefficient
      elseif n == 1
        tensor[sites[n] => alpha, links[n] => next_state] += coefficient
      elseif n == nsites
        tensor[links[n - 1] => previous, sites[n] => alpha] += coefficient
      else
        tensor[links[n - 1] => previous, sites[n] => alpha,
          links[n] => next_state] += coefficient
      end
    end
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

_pxp_pair_state(row_state::Integer, column_state::Integer) =
  2 * (Int(row_state) - 1) + Int(column_state)

function _pxp_constraint_pair_transfer(alpha, beta, basis)
  transfer = zeros(ComplexF64, 4, 4)
  for previous_row in 1:2, previous_column in 1:2, row in 1:2, column in 1:2
    next_row, allowed_row = _pxp_constraint_transition(previous_row, row)
    next_column, allowed_column = _pxp_constraint_transition(previous_column, column)
    transfer[_pxp_pair_state(previous_row, previous_column),
      _pxp_pair_state(next_row, next_column)] +=
      conj(basis[alpha][row, column]) * basis[beta][row, column] *
      allowed_row * allowed_column
  end
  return transfer
end

function _pxp_pair_basis_transform(basis)
  transform = Matrix{ComplexF64}(undef, 4, 4)
  for row in 1:2, column in 1:2, alpha in eachindex(basis)
    transform[_pxp_pair_state(row, column), alpha] = basis[alpha][row, column]
  end
  @assert transform' * transform ≈ Matrix{ComplexF64}(I, 4, 4)
  return transform
end

function _pauli_pxp_constraint_projector_real(sites)
  nsites = length(sites)
  basis = operator_basis_matrices(2)
  transform = _pxp_pair_basis_transform(basis)
  initial = ComplexF64[1, 0, 0, 0]
  terminal = ones(ComplexF64, 4)
  links = [Index(4, "PXPConstraintProjectorLink,n=$(n)") for n in 1:(nsites - 1)]
  tensors = ITensor[]
  for n in 1:nsites
    tensor = nsites == 1 ? ITensor(Float64, prime(sites[n]), sites[n]) :
      n == 1 ? ITensor(Float64, prime(sites[n]), sites[n], links[n]) :
      n == nsites ? ITensor(Float64, links[n - 1], prime(sites[n]), sites[n]) :
      ITensor(Float64, links[n - 1], prime(sites[n]), sites[n], links[n])
    for alpha in eachindex(basis), beta in eachindex(basis)
      old_transfer = _pxp_constraint_pair_transfer(alpha, beta, basis)
      transformed = if nsites == 1
        transpose(initial) * old_transfer * terminal
      elseif n == 1
        transpose(initial) * old_transfer * transform
      elseif n == nsites
        transform' * old_transfer * terminal
      else
        transform' * old_transfer * transform
      end
      @assert norm(imag(transformed)) <= 32eps(Float64)
      values = real(transformed)
      if nsites == 1
        tensor[prime(sites[n]) => alpha, sites[n] => beta] = only(values)
      elseif n == 1
        for right in 1:4
          tensor[prime(sites[n]) => alpha, sites[n] => beta,
            links[n] => right] = values[right]
        end
      elseif n == nsites
        for left in 1:4
          tensor[links[n - 1] => left, prime(sites[n]) => alpha,
            sites[n] => beta] = values[left]
        end
      else
        for left in 1:4, right in 1:4
          tensor[links[n - 1] => left, prime(sites[n]) => alpha,
            sites[n] => beta, links[n] => right] = values[left, right]
        end
      end
    end
    push!(tensors, tensor)
  end
  return MPO(tensors)
end

function _realify_pxp_superoperator_mpo(op::MPO, sites)
  nsites = length(op)
  length(sites) == nsites || throw(ArgumentError("MPO and site counts must match"))
  transform = _pxp_pair_basis_transform(operator_basis_matrices(2))
  old_links = [commonind(op[n], op[n + 1]) for n in 1:(nsites - 1)]
  all(dim(link) == 4 for link in old_links) || throw(ArgumentError(
    "PXP superoperator realification requires bond-dimension-4 pair links"))
  new_links = [Index(4, "PXPRealSuperLink,n=$(n)") for n in 1:(nsites - 1)]
  tensors = ITensor[]
  for n in 1:nsites
    tensor = nsites == 1 ? ITensor(Float64, prime(sites[n]), sites[n]) :
      n == 1 ? ITensor(Float64, prime(sites[n]), sites[n], new_links[n]) :
      n == nsites ? ITensor(Float64, new_links[n - 1], prime(sites[n]), sites[n]) :
      ITensor(Float64, new_links[n - 1], prime(sites[n]), sites[n], new_links[n])
    for alpha in 1:4, beta in 1:4
      transformed = if nsites == 1
        op[n][prime(sites[n]) => alpha, sites[n] => beta]
      elseif n == 1
        values = ComplexF64[op[n][prime(sites[n]) => alpha, sites[n] => beta,
          old_links[n] => right] for right in 1:4]
        transpose(values) * transform
      elseif n == nsites
        values = ComplexF64[op[n][old_links[n - 1] => left,
          prime(sites[n]) => alpha, sites[n] => beta] for left in 1:4]
        transform' * values
      else
        values = ComplexF64[op[n][old_links[n - 1] => left,
          prime(sites[n]) => alpha, sites[n] => beta, old_links[n] => right]
          for left in 1:4, right in 1:4]
        transform' * values * transform
      end
      @assert norm(imag(transformed)) <= 64eps(Float64)
      values = real(transformed)
      if nsites == 1
        tensor[prime(sites[n]) => alpha, sites[n] => beta] = values
      elseif n == 1
        for right in 1:4
          tensor[prime(sites[n]) => alpha, sites[n] => beta,
            new_links[n] => right] = values[right]
        end
      elseif n == nsites
        for left in 1:4
          tensor[new_links[n - 1] => left, prime(sites[n]) => alpha,
            sites[n] => beta] = values[left]
        end
      else
        for left in 1:4, right in 1:4
          tensor[new_links[n - 1] => left, prime(sites[n]) => alpha,
            sites[n] => beta, new_links[n] => right] = values[left, right]
        end
      end
    end
    push!(tensors, tensor)
  end
  return MPO(tensors)
end

"""
    pauli_pxp_constraint_state(sites)

Return the vectorized PXP constraint projector `P_G = Π_j (1 - n_j n_{j+1})` as a Pauli-basis
operator-space `MPS` with bond dimension 2.

Up to normalization this is the infinite-temperature density matrix of the **constrained**
sector, `ρ_∞ ∝ P_G`, and is the natural `initial_state` for thermal preparation with
[`pauli_gibbs_state`](@ref) in PXP transport calculations: the unconstrained identity state
carries weight on blockade-violating configurations that the PXP Hamiltonian never couples
to physical dynamics.

# Arguments
- `sites`: Pauli-space site indices (dimension 4), typically from [`pauli_siteinds`](@ref).

# Returns
- An `MPS` over `sites` whose Pauli amplitudes are `tr(P_α† P_G)`.

# Notes
- Defined for spin-1/2 (local dimension 2) operator space only; throws `ArgumentError` for
  any other local dimension, since the PXP blockade constraint is a two-level construction.
- Constructed directly in the real I/Z coefficient automaton, so every core is `Float64` and no
  intermediate `ComplexF64` physical MPO is vectorized.
"""
function pauli_pxp_constraint_state(sites)
  all(local_dimension(site) == 2 for site in sites) ||
    throw(ArgumentError("this helper is defined for spin-1/2 (local dimension 2) operator space only"))
  return _pauli_pxp_constraint_state_real(sites)
end

"""
    pauli_pxp_constraint_projector(sites)

Return the operator-space MPO implementing the two-sided PXP constraint projection
`ρ ↦ P_G ρ P_G` in the normalized Pauli basis, with bond dimension 4.

Every PXP term commutes with `P_G`, so exact evolution never leaves the constrained sector;
truncation (DMT or plain SVD) does leak weight out of it. Periodically applying this MPO
("checkpoints") removes the leaked weight — see [`constrained_dmt_evolve!`](@ref).

# Arguments
- `sites`: Pauli-space site indices (dimension 4).

# Returns
- An `MPO` over `sites`, idempotent up to truncation, that fixes any vectorized operator
  supported in the constrained sector and annihilates operators supported on
  blockade-violating configurations.

# Notes
- Defined for spin-1/2 (local dimension 2) operator space only; throws `ArgumentError` for
  any other local dimension, since the PXP blockade constraint is a two-level construction.
- The paired computational virtual states are transformed to the normalized Pauli/Hermitian
  basis before real storage. This keeps bond dimension 4 and makes every core `Float64`; it is
  not a per-core downcast of an arbitrary complex gauge.
"""
function pauli_pxp_constraint_projector(sites)
  all(local_dimension(site) == 2 for site in sites) ||
    throw(ArgumentError("this helper is defined for spin-1/2 (local dimension 2) operator space only"))
  return _pauli_pxp_constraint_projector_real(sites)
end

"""
    PXPControlledGate

An exact local PXP superoperator MPO together with its bound center, support and duration.
Construct one with [`pauli_pxp_controlled_gate`](@ref).
"""
struct PXPControlledGate{TM}
  mpo::TM
  center::Int
  start::Int
  span::Int
  duration::Float64
  omega::Float64
  mu::Float64
end

function _two_product_mpo(sites, first_factors, second_factors)
  length(sites) == length(first_factors) == length(second_factors) ||
    throw(ArgumentError("site and factor counts must match"))
  nsites = length(sites)
  nsites >= 1 || throw(ArgumentError("product-sum MPO requires at least one site"))
  if nsites == 1
    return MPO([_dense_local_operator([sites[1]], first_factors[1] + second_factors[1])])
  end
  links = [Index(2, "PXPControlledLink,n=$(n)") for n in 1:(nsites - 1)]
  tensors = ITensor[]
  for n in 1:nsites
    local_tensor = if n == 1
      ITensor(ComplexF64, prime(sites[n]), dag(sites[n]), links[n])
    elseif n == nsites
      ITensor(ComplexF64, links[n - 1], prime(sites[n]), dag(sites[n]))
    else
      ITensor(ComplexF64, links[n - 1], prime(sites[n]), dag(sites[n]), links[n])
    end
    for row in 1:2, column in 1:2
      if n == 1
        local_tensor[prime(sites[n]) => row, sites[n] => column, links[n] => 1] =
          first_factors[n][row, column]
        local_tensor[prime(sites[n]) => row, sites[n] => column, links[n] => 2] =
          second_factors[n][row, column]
      elseif n == nsites
        local_tensor[links[n - 1] => 1, prime(sites[n]) => row, sites[n] => column] =
          first_factors[n][row, column]
        local_tensor[links[n - 1] => 2, prime(sites[n]) => row, sites[n] => column] =
          second_factors[n][row, column]
      else
        local_tensor[links[n - 1] => 1, prime(sites[n]) => row,
          sites[n] => column, links[n] => 1] = first_factors[n][row, column]
        local_tensor[links[n - 1] => 2, prime(sites[n]) => row,
          sites[n] => column, links[n] => 2] = second_factors[n][row, column]
      end
    end
    push!(tensors, local_tensor)
  end
  return MPO(tensors)
end

"""
    pauli_pxp_controlled_gate(sites, center, duration; omega=1.0, mu=0.0)

Build the exact local PXP plus chemical-potential evolution gate as an operator-space MPO bound
to `sites`. For `V=exp(-im*duration*mu*n)` and
`W=exp(-im*duration*(omega*X+mu*n))`, its physical representation is
`I_left ⊗ V ⊗ I_right + P_left ⊗ (W-V) ⊗ P_right`, omitting a missing edge projector.
The physical MPO bond dimension is 2 and the returned superoperator bond dimension is at most 4.

The paired virtual links are transformed from the computational pair basis to the normalized
Pauli/Hermitian basis before storage. Thus every returned core is `Float64` without increasing
the bond dimension or discarding an imaginary part from an arbitrary gauge.
"""
function pauli_pxp_controlled_gate(sites, center::Integer, duration::Real;
                                   omega::Real=1.0, mu::Real=0.0)
  nsites = length(sites)
  isfinite(duration) || throw(ArgumentError("PXP controlled gate requires finite duration"))
  isfinite(omega) || throw(ArgumentError("PXP controlled gate requires finite omega"))
  isfinite(mu) || throw(ArgumentError("PXP controlled gate requires finite mu"))
  all(local_dimension(site) == 2 for site in sites) || throw(ArgumentError(
    "PXP controlled gates require spin-1/2 operator-space sites"))
  support = pxp_term_support(nsites, center)
  local_sites = collect(sites[support])
  physical_sites = _pxp_physical_sites(length(support))
  paulis = pauli_matrices()
  projector = Matrix{ComplexF64}((paulis.I + paulis.Z) / 2)
  number = Matrix{ComplexF64}((paulis.I - paulis.Z) / 2)
  identity_matrix = Matrix{ComplexF64}(I, 2, 2)
  onsite = exp(-im * Float64(duration) * Float64(mu) * number)
  controlled = exp(-im * Float64(duration) *
    (Float64(omega) * Matrix{ComplexF64}(paulis.X) + Float64(mu) * number))
  center_offset = Int(center) - first(support) + 1
  first_factors = [offset == center_offset ? onsite : identity_matrix
                   for offset in 1:length(support)]
  second_factors = [offset == center_offset ? controlled - onsite : projector
                    for offset in 1:length(support)]
  physical_mpo = _two_product_mpo(physical_sites, first_factors, second_factors)
  complex_superoperator = pauli_superoperator_mpo(physical_mpo, local_sites)
  superoperator = _realify_pxp_superoperator_mpo(complex_superoperator, local_sites)
  return PXPControlledGate(superoperator, Int(center), first(support), length(support),
    Float64(duration), Float64(omega), Float64(mu))
end

"""One immutable entry of a compiled PXP DMT product-formula plan."""
struct PXPDMTPlanEntry
  id::Int
  center::Int
  start::Int
  span::Int
  duration::Float64
  layer::Int
  direction::Symbol
end

"""
    PXPDMTPlan

A compiled PXP product-formula plan, including inspectable entries, its physical step `tau`, a
stable algorithm signature, and the [`DMTGateEvolution`](@ref) used to execute it.
"""
struct PXPDMTPlan{TE}
  scheme::Symbol
  tau::Float64
  entries::Vector{PXPDMTPlanEntry}
  signature::String
  evolution::TE
end

function _pxp_dmt_config_signature(maxdim, cutoff, preserve_diameter, preserve_operators,
                                   truncation, normalize)
  preserve_digest = bytes2hex(sha256(codeunits(repr(preserve_operators))))
  return "maxdim=$(Int(maxdim)):cutoff=$(Float64(cutoff)):" *
    "diameter=$(Int(preserve_diameter)):preserve=$(preserve_digest):" *
    "truncation=$(truncation):normalize=$(Bool(normalize))"
end

function _pxp_parity_dmt_plan(
  scheme::Symbol,
  sites,
  tau::Real,
  layers;
  omega::Real,
  mu::Real,
  nstep::Integer,
  maxdim::Integer,
  cutoff::Real,
  preserve_diameter::Integer,
  preserve_operators,
  truncation::Symbol,
  normalize::Bool,
)
  length(sites) >= 2 || throw(ArgumentError("PXP $(scheme) plan requires at least two sites"))
  isfinite(tau) && tau > 0 || throw(ArgumentError("PXP $(scheme) plan requires finite tau > 0"))
  isfinite(omega) || throw(ArgumentError("PXP $(scheme) plan requires finite omega"))
  isfinite(mu) || throw(ArgumentError("PXP $(scheme) plan requires finite mu"))
  all(local_dimension(site) == 2 for site in sites) || throw(ArgumentError(
    "PXP $(scheme) plan requires spin-1/2 operator-space sites"))
  _validate_dmt_budget_dimension(2, maxdim, preserve_diameter, preserve_operators)
  entries = PXPDMTPlanEntry[]
  for (layer, (parity, duration)) in enumerate(layers)
    centers = parity === :odd ? (1:2:length(sites)) : (2:2:length(sites))
    for center in centers
      support = pxp_term_support(length(sites), center)
      push!(entries, PXPDMTPlanEntry(length(entries) + 1, center, first(support),
        length(support), Float64(duration), layer, :R))
    end
  end
  gates = [pauli_pxp_controlled_gate(sites, entry.center, entry.duration; omega=omega, mu=mu)
    for entry in entries]
  schedule = [entry.start for entry in entries]
  evolution = DMTGateEvolution(gates, tau; schedule=schedule, reverse_schedule=Int[],
    nstep=nstep, maxdim=maxdim, cutoff=cutoff, gate_maxdim=0,
    preserve_diameter=preserve_diameter, preserve_operators=preserve_operators,
    truncation=truncation, gate_backend=:controlled, normalize=normalize)
  layer_signature = join(("$(parity):$(Float64(duration))" for (parity, duration) in layers), ",")
  config_signature = _pxp_dmt_config_signature(maxdim, cutoff, preserve_diameter,
    preserve_operators, truncation, normalize)
  signature = "pxp-$(lowercase(String(scheme)))-v2:N=$(length(sites)):tau=$(Float64(tau)):" *
              "omega=$(Float64(omega)):mu=$(Float64(mu)):$(layer_signature):direction=R:" *
              config_signature
  return PXPDMTPlan(scheme, Float64(tau), entries, signature, evolution)
end

"""
    pxp_s2_dmt_plan(sites, tau; omega=1.0, nstep=1, kwargs...)

Compile the second-order parity schedule `O(tau/2) E(tau) O(tau/2)` using exact
[`PXPControlledGate`](@ref) entries. The generated DMT evolution has no reverse sweep; every
entry records `direction=:R`, and one physical step attempts exactly `3(length(sites)-1)` bond
updates.
"""
function pxp_s2_dmt_plan(
  sites,
  tau::Real;
  omega::Real=1.0,
  mu::Real=0.0,
  nstep::Integer=1,
  maxdim::Integer=30,
  cutoff::Real=1e-12,
  preserve_diameter::Integer=3,
  preserve_operators=nothing,
  truncation::Symbol=:dense,
  normalize::Bool=true,
)
  layers = ((:odd, tau / 2), (:even, tau), (:odd, tau / 2))
  return _pxp_parity_dmt_plan(:S2, sites, tau, layers; omega=omega, mu=mu, nstep=nstep,
    maxdim=maxdim, cutoff=cutoff, preserve_diameter=preserve_diameter,
    preserve_operators=preserve_operators, truncation=truncation, normalize=normalize)
end

"""
    pxp_s4_dmt_plan(sites, tau; omega=1.0, nstep=1, kwargs...)

Compile the fourth-order seven-layer Yoshida composition of S2. The coefficients are evaluated
from `cbrt(2)` rather than decimal approximations. Negative middle-layer durations are preserved
as gate-construction data, and one physical step attempts exactly `7(length(sites)-1)` DMT bond
updates.
"""
function pxp_s4_dmt_plan(
  sites,
  tau::Real;
  omega::Real=1.0,
  mu::Real=0.0,
  nstep::Integer=1,
  maxdim::Integer=30,
  cutoff::Real=1e-12,
  preserve_diameter::Integer=3,
  preserve_operators=nothing,
  truncation::Symbol=:dense,
  normalize::Bool=true,
)
  root2 = cbrt(2.0)
  w1 = inv(2.0 - root2)
  w0 = -root2 * w1
  layers = ((:odd, tau * w1 / 2),
            (:even, tau * w1),
            (:odd, tau * (w1 + w0) / 2),
            (:even, tau * w0),
            (:odd, tau * (w0 + w1) / 2),
            (:even, tau * w1),
            (:odd, tau * w1 / 2))
  return _pxp_parity_dmt_plan(:S4, sites, tau, layers; omega=omega, mu=mu, nstep=nstep,
    maxdim=maxdim, cutoff=cutoff, preserve_diameter=preserve_diameter,
    preserve_operators=preserve_operators, truncation=truncation, normalize=normalize)
end

"""Return the number of internal-bond DMT updates attempted per physical plan step."""
pxp_dmt_attempted_updates(plan::PXPDMTPlan) = sum(entry.span - 1 for entry in plan.entries)

"""An exact odd or even PXP layer represented by a full-chain operator-space MPO."""
struct PXPParityLayer{TM}
  mpo::TM
  parity::Symbol
  duration::Float64
  omega::Float64
  mu::Float64
  signature::String
end

function _identity_operator_mpo(sites)
  return MPO([delta(prime(site), dag(site)) for site in sites])
end

function _embed_local_mpo(op::MPO, sites, start::Integer)
  tensors = [delta(prime(site), dag(site)) for site in sites]
  for offset in 1:length(op)
    tensors[Int(start) + offset - 1] = op[offset]
  end
  return MPO(tensors)
end

"""
    pauli_pxp_parity_layer(sites, parity, duration; omega=1.0)

Build the exact full-chain product of all mutually commuting PXP gates on one center parity.
The returned real operator-space MPO has bond dimension at most 4 and is bound to the supplied
site indices. Construction cost is preparation work; reuse the returned layer across steps.
"""
function pauli_pxp_parity_layer(sites, parity::Symbol, duration::Real;
                                omega::Real=1.0, mu::Real=0.0)
  parity in (:odd, :even) || throw(ArgumentError("PXP layer parity must be :odd or :even"))
  length(sites) >= 2 || throw(ArgumentError("PXP parity layer requires at least two sites"))
  isfinite(duration) || throw(ArgumentError("PXP parity layer requires finite duration"))
  isfinite(omega) || throw(ArgumentError("PXP parity layer requires finite omega"))
  isfinite(mu) || throw(ArgumentError("PXP parity layer requires finite mu"))
  layer = _identity_operator_mpo(sites)
  centers = parity === :odd ? (1:2:length(sites)) : (2:2:length(sites))
  for center in centers
    gate = pauli_pxp_controlled_gate(sites, center, duration; omega=omega, mu=mu)
    embedded = _embed_local_mpo(gate.mpo, sites, gate.start)
    layer = apply(embedded, layer; maxdim=4, cutoff=0.0)
  end
  maxlinkdim(layer) <= 4 || error("PXP parity-layer MPO exceeded exact bond dimension 4")
  signature = "pxp-layer-v1:N=$(length(sites)):parity=$(parity):" *
              "duration=$(Float64(duration)):omega=$(Float64(omega)):mu=$(Float64(mu))"
  return PXPParityLayer(layer, parity, Float64(duration), Float64(omega), Float64(mu), signature)
end

"""A layerwise PXP DMT plan that applies each parity MPO before one directed DMT sweep."""
struct PXPLayerDMTPlan{TL,TO}
  scheme::Symbol
  tau::Float64
  layers::TL
  options::TO
  nstep::Int
  normalize::Bool
  direction::Symbol
  application_backend::Symbol
  signature::String
end

"""
    pxp_layer_dmt_plan(sites, tau; scheme=:S2, kwargs...)

Compile an experimental S2 or S4 layerwise plan. Each exact parity-layer MPO is applied in one
operation and followed by a directed DMT sweep over all `L-1` bonds.
"""
function pxp_layer_dmt_plan(
  sites,
  tau::Real;
  scheme::Symbol=:S2,
  omega::Real=1.0,
  mu::Real=0.0,
  nstep::Integer=1,
  maxdim::Integer=30,
  cutoff::Real=1e-12,
  preserve_diameter::Integer=3,
  preserve_operators=nothing,
  truncation::Symbol=:dense,
  normalize::Bool=true,
  direction::Symbol=:R,
)
  nstep >= 1 || throw(ArgumentError("PXP layerwise plan requires nstep >= 1"))
  direction in (:R, :L) || throw(ArgumentError("PXP layerwise direction must be :R or :L"))
  gate_plan = scheme === :S2 ? pxp_s2_dmt_plan(sites, tau; omega=omega, mu=mu,
    nstep=nstep, maxdim=maxdim, cutoff=cutoff, preserve_diameter=preserve_diameter,
    preserve_operators=preserve_operators, truncation=truncation, normalize=normalize) :
    scheme === :S4 ? pxp_s4_dmt_plan(sites, tau; omega=omega, mu=mu,
    nstep=nstep, maxdim=maxdim, cutoff=cutoff, preserve_diameter=preserve_diameter,
    preserve_operators=preserve_operators, truncation=truncation, normalize=normalize) :
    throw(ArgumentError("PXP layerwise scheme must be :S2 or :S4"))
  layer_entries = [first(filter(entry -> entry.layer == layer, gate_plan.entries))
    for layer in 1:maximum(entry.layer for entry in gate_plan.entries)]
  layers = [pauli_pxp_parity_layer(sites, isodd(entry.center) ? :odd : :even,
    entry.duration; omega=omega, mu=mu) for entry in layer_entries]
  options = DMTOptions(maxdim=maxdim, cutoff=cutoff, gate_maxdim=0,
    preserve_diameter=preserve_diameter, preserve_operators=preserve_operators,
    truncation=truncation)
  config_signature = _pxp_dmt_config_signature(maxdim, cutoff, preserve_diameter,
    preserve_operators, truncation, normalize)
  signature = "pxp-layerwise-$(lowercase(String(scheme)))-v2:N=$(length(sites)):" *
              "tau=$(Float64(tau)):omega=$(Float64(omega)):mu=$(Float64(mu)):" *
              "direction=$(direction):application=naive:" * config_signature
  return PXPLayerDMTPlan(scheme, Float64(tau), layers, options, Int(nstep),
    Bool(normalize), direction, :naive, signature)
end

"""Return the number of directed bond truncations attempted per layerwise physical step."""
pxp_dmt_attempted_updates(plan::PXPLayerDMTPlan) =
  length(plan.layers) * (length(first(plan.layers).mpo) - 1)
