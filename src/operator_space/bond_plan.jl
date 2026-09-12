"""One immutable gate application in a compiled nearest-neighbor DMT plan."""
struct BondDMTPlanEntry
  id::Int
  bond::Int
  term::Int
  duration::Float64
  layer::Int
  direction::Symbol
end

"""
    BondDMTPlan

A compiled second- or fourth-order product formula for open-boundary nearest-neighbor
Hamiltonians. Use [`bond_dmt_plan`](@ref) to construct one.
"""
struct BondDMTPlan{TS,TG,TO}
  order::Int
  composition::Symbol
  tau::Float64
  sites::TS
  entries::Vector{BondDMTPlanEntry}
  gates::TG
  options::TO
  nstep::Int
  normalize::Bool
  direction::Symbol
  signature::String
end


function _bond_formula_layers(tau::Float64, order::Int, composition::Symbol)
  order == 2 && return ((:odd, tau / 2), (:even, tau), (:odd, tau / 2))
  order == 4 || throw(ArgumentError("bond DMT order must be 2 or 4"))
  w = inv(2.0 - cbrt(2.0))
  v = 1.0 - 2.0w
  if composition === :merged7
    return ((:odd, tau * w / 2), (:even, tau * w),
      (:odd, tau * (w + v) / 2), (:even, tau * v),
      (:odd, tau * (v + w) / 2), (:even, tau * w), (:odd, tau * w / 2))
  elseif composition === :yoshida9
    return ((:odd, tau * w / 2), (:even, tau * w), (:odd, tau * w / 2),
      (:odd, tau * v / 2), (:even, tau * v), (:odd, tau * v / 2),
      (:odd, tau * w / 2), (:even, tau * w), (:odd, tau * w / 2))
  end
  throw(ArgumentError("fourth-order bond DMT composition must be :merged7 or :yoshida9"))
end

function _bond_term_eigensystems(bond_terms, d::Int)
  unique_terms = Matrix{ComplexF64}[]
  term_ids = Int[]
  eigensystems = Any[]
  for (bond, term) in pairs(bond_terms)
    size(term) == (d^2, d^2) || throw(ArgumentError(
      "bond term $(bond) must have size $((d^2, d^2)), got $(size(term))"))
    dense = Matrix{ComplexF64}(term)
    norm(dense - dense') <= sqrt(eps(Float64)) * max(norm(dense), 1.0) ||
      throw(ArgumentError("bond term $(bond) must be Hermitian"))
    term_id = findfirst(existing -> existing == dense, unique_terms)
    if isnothing(term_id)
      push!(unique_terms, dense)
      push!(eigensystems, eigen(Hermitian(dense)))
      term_id = length(unique_terms)
    end
    push!(term_ids, term_id)
  end
  return term_ids, eigensystems
end

function _bond_plan_signature(order, composition, nsites, tau, direction, options, normalize,
                              bond_terms)
  term_digest = bytes2hex(sha256(codeunits(repr(bond_terms))))
  preserve_digest = bytes2hex(sha256(codeunits(repr(options.preserve_operators))))
  return "bond-s$(order)-v1:N=$(nsites):tau=$(tau):composition=$(composition):" *
    "direction=$(direction):backend=$(options.gate_backend):maxdim=$(options.maxdim):" *
    "cutoff=$(options.cutoff):diameter=$(options.preserve_diameter):" *
    "preserve=$(preserve_digest):truncation=$(options.truncation):normalize=$(normalize):" *
    "terms=$(term_digest)"
end

"""
    bond_dmt_plan(sites, bond_terms, tau; order=4, composition=:merged7, kwargs...)

Compile a DMT product-formula plan for an open chain. `bond_terms[b]` is the Hermitian two-site
Hamiltonian on bond `b`, including any desired allocation of on-site fields. Here `tau` is the
duration of one complete product-formula step, not the duration of an individual sublayer.

Second order uses odd-even-odd Strang splitting. Fourth order uses either the seven-layer merged
Yoshida composition or the unmerged nine-layer composition; these remain distinct
finite-bond-dimension algorithms because every gate application is followed by DMT.

# Backend recommendation
- Keep the default `gate_backend=:qr` for production runs and unknown input rank profiles.
- Use `gate_backend=:direct` only for dense two-site gates with `truncation=:dense`, after an A/B
  check on a representative checkpoint. It skips the full gated-center SVD, but still performs
  the complement SVD, final low-rank refactorization, and every intermediate DMT truncation.
- `:fused` retains the full gated-center SVD and is primarily useful as a comparison path.
- `:product` retains the ITensorMPS gate-application path for compatibility.

The default remains `order=4`, `composition=:merged7`, `gate_backend=:qr`. Use `order=2` as a
diagnostic/control; choose `tau` from an observable-error-versus-wall convergence study rather
than assuming that fourth order is always cheaper at finite `maxdim`.

# Example
```julia
plan = bond_dmt_plan(sites, bond_terms, 0.05;
  maxdim=64, normalize=false) # recommended starting point: S4 merged7 + QR

direct_plan = bond_dmt_plan(sites, bond_terms, 0.05;
  maxdim=64, normalize=false, gate_backend=:direct, truncation=:dense)
```
"""
function bond_dmt_plan(
  sites,
  bond_terms,
  tau::Real;
  order::Integer=4,
  composition::Symbol=:merged7,
  nstep::Integer=1,
  maxdim::Integer=30,
  cutoff::Real=1e-12,
  preserve_diameter::Integer=3,
  preserve_operators=nothing,
  truncation::Symbol=:dense,
  # QR is intentionally the stable plan default. `:direct` remains explicit until production
  # checkpoint, peak-RSS, and broader rank-profile evidence justify reconsidering this choice.
  gate_backend::Symbol=:qr,
  normalize::Bool=true,
  direction::Symbol=:R,
)
  nsites = length(sites)
  nsites >= 2 || throw(ArgumentError("bond DMT plan requires at least two sites"))
  length(bond_terms) == nsites - 1 || throw(ArgumentError(
    "bond DMT plan requires one term for each of the $(nsites - 1) open-chain bonds"))
  isfinite(tau) && tau > 0 || throw(ArgumentError("bond DMT plan requires finite tau > 0"))
  nstep >= 1 || throw(ArgumentError("bond DMT plan requires nstep >= 1"))
  direction in (:R, :L) || throw(ArgumentError("bond DMT direction must be :R or :L"))
  order in (2, 4) || throw(ArgumentError("bond DMT order must be 2 or 4"))
  order == 2 && composition !== :merged7 && throw(ArgumentError(
    "composition only applies to fourth-order bond DMT plans"))
  d = local_dimension(first(sites))
  all(local_dimension(site) == d for site in sites) || throw(ArgumentError(
    "bond DMT sites must have a uniform local dimension"))
  options = DMTOptions(maxdim=maxdim, cutoff=cutoff, gate_maxdim=0,
    preserve_diameter=preserve_diameter, preserve_operators=preserve_operators,
    truncation=truncation, gate_backend=gate_backend)
  _validate_dmt_budget_dimension(d, maxdim, preserve_diameter, preserve_operators)
  gate_backend in (:product, :qr, :fused, :direct) || throw(ArgumentError(
    "bond DMT gate_backend must be :product, :qr, :fused, or :direct"))
  gate_backend === :direct && truncation !== :dense && throw(ArgumentError(
    "bond DMT gate_backend=:direct currently requires truncation=:dense"))
  layers = _bond_formula_layers(Float64(tau), Int(order), composition)
  term_ids, eigensystems = _bond_term_eigensystems(bond_terms, d)
  entries = BondDMTPlanEntry[]
  gates = Matrix{Float64}[]
  gate_cache = Dict{Tuple{Int,Float64},Matrix{Float64}}()
  for (layer, (parity, duration)) in enumerate(layers)
    bonds = parity === :odd ? collect(1:2:(nsites - 1)) : collect(2:2:(nsites - 1))
    direction === :L && reverse!(bonds)
    for bond in bonds
      term_id = term_ids[bond]
      key = (term_id, Float64(duration))
      gate = get!(gate_cache, key) do
        decomposition = eigensystems[term_id]
        unitary = decomposition.vectors *
          Diagonal(exp.(-im * Float64(duration) .* decomposition.values)) *
          decomposition.vectors'
        operator_gate(unitary; d=d)
      end
      push!(entries, BondDMTPlanEntry(length(entries) + 1, bond, term_id,
        Float64(duration), layer, direction))
      push!(gates, gate)
    end
  end
  signature = _bond_plan_signature(order, composition, nsites, Float64(tau), direction,
    options, normalize, bond_terms)
  return BondDMTPlan(Int(order), composition, Float64(tau), collect(sites), entries, gates,
    options, Int(nstep), normalize, direction, signature)
end

"""Return the number of two-site DMT updates attempted per physical plan step."""
bond_dmt_attempted_updates(plan::BondDMTPlan) = length(plan.entries)

function _validate_bond_dmt_plan(psi::MPS, plan::BondDMTPlan)
  length(psi) == length(plan.sites) ||
    throw(ArgumentError("bond DMT plan and state must have matching lengths"))
  for site in eachindex(plan.sites)
    siteind(psi, site) == plan.sites[site] ||
      throw(ArgumentError("bond DMT plan site indices do not match the target state"))
  end
  _validate_dmt_budget(psi, plan.options.maxdim, plan.options.preserve_diameter,
    plan.options.preserve_operators)
  for (entry, gate) in zip(plan.entries, plan.gates)
    _validate_dmt_step(psi, gate, entry.bond, 2, entry.direction, plan.options.maxdim,
      plan.options.preserve_diameter, plan.options.preserve_operators)
  end
  return nothing
end

function dmt_evolve!(psi::MPS, plan::BondDMTPlan; normalize::Bool=plan.normalize)
  _validate_bond_dmt_plan(psi, plan)
  cache = _DMTEnvCache(psi)
  for _ in 1:plan.nstep, (entry, gate) in zip(plan.entries, plan.gates)
    dmt_step!(psi, gate, entry.bond; maxdim=plan.options.maxdim,
      cutoff=plan.options.cutoff, direction=entry.direction, gate_maxdim=0,
      preserve_diameter=plan.options.preserve_diameter,
      preserve_operators=plan.options.preserve_operators,
      truncation=plan.options.truncation, gate_backend=plan.options.gate_backend, cache=cache)
  end
  normalize && normalize!(psi)
  return psi
end

evolve!(psi::MPS, plan::BondDMTPlan; normalize::Bool=plan.normalize) =
  dmt_evolve!(psi, plan; normalize=normalize)
