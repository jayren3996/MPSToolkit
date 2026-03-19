"""
    pauli_siteinds(nsites; tagprefix="PauliSpace")

Construct a vector of dimension-4 site indices suitable for vectorized local Pauli bases.
The basis ordering is `(I, X, Y, Z)`.

# Arguments
- `nsites`: Number of operator-space sites.

# Keyword Arguments
- `tagprefix`: Prefix used when naming the generated `Index` tags.

# Returns
- A vector of length `nsites` containing dimension-4 `Index` objects.
"""
function pauli_siteinds(nsites::Integer; tagprefix::AbstractString="PauliSpace")
  nsites < 1 && throw(ArgumentError("number of Pauli-basis sites must be positive"))
  return [Index(4, "$(tagprefix),n=$(n)") for n in 1:Int(nsites)]
end

"""
    _daoe_nontrivial(local_state)

Return `true` if the local Pauli-basis label contributes to the DAOE Pauli weight.
"""
function _daoe_nontrivial(local_state::Int)
  return local_state != 1
end

"""
    _daoe_transition(weight, local_state, lstar, gamma)

Return the next DAOE virtual weight and local decay factor for one Pauli-basis symbol.

# Returns
- `(next_weight, coeff)` for the diagonal MPO transition rule.
"""
function _daoe_transition(weight::Int, local_state::Int, lstar::Int, gamma::Real)
  if !_daoe_nontrivial(local_state)
    return weight, 1.0
  elseif weight < lstar
    return weight + 1, 1.0
  else
    return lstar, exp(-gamma)
  end
end

"""
    _fdaoe_parity(state, wstar)

Return the fermion-parity flag encoded by an FDAOE virtual state.
"""
function _fdaoe_parity(state::Int, wstar::Int)
  if state <= wstar
    return isodd(state - 1)
  else
    return state == wstar + 2
  end
end

"""
    _fdaoe_additional_weight(local_state, odd_parity)

Return the fermionic weight increment associated with one local Pauli-basis symbol.
"""
function _fdaoe_additional_weight(local_state::Int, odd_parity::Bool)
  if local_state == 1
    return odd_parity ? 2 : 0
  elseif local_state == 4
    return odd_parity ? 0 : 2
  else
    return 1
  end
end

"""
    _fdaoe_flips_parity(local_state)

Return `true` if the local Pauli-basis symbol flips the Jordan-Wigner parity tracker.
"""
function _fdaoe_flips_parity(local_state::Int)
  return local_state == 2 || local_state == 3
end

"""
    _fdaoe_target_state(weight, odd_parity, wstar)

Map a fermionic weight/parity pair onto an FDAOE virtual-state index.
"""
function _fdaoe_target_state(weight::Int, odd_parity::Bool, wstar::Int)
  if weight < wstar
    return weight + 1
  else
    return odd_parity ? (wstar + 2) : (wstar + 1)
  end
end

"""
    _fdaoe_transition(state, local_state, wstar, gamma)

Return the next FDAOE virtual state and local decay factor for one Pauli-basis symbol.

# Returns
- `(next_state, coeff)` for the diagonal MPO transition rule.
"""
function _fdaoe_transition(state::Int, local_state::Int, wstar::Int, gamma::Real)
  odd_parity = _fdaoe_parity(state, wstar)
  additional_weight = _fdaoe_additional_weight(local_state, odd_parity)
  next_parity = _fdaoe_flips_parity(local_state) ? !odd_parity : odd_parity

  if state <= wstar
    current_weight = state - 1
    next_weight = current_weight + additional_weight
    next_state = _fdaoe_target_state(next_weight, next_parity, wstar)
    decay = next_weight > wstar ? exp(-gamma * (next_weight - wstar)) : 1.0
    return next_state, decay
  else
    next_state = _fdaoe_target_state(wstar, next_parity, wstar)
    return next_state, exp(-gamma * additional_weight)
  end
end

"""
    _diagonal_projector_mpo(sites, bond_dim, transition_rule)

Build a diagonal MPO from a virtual-state transition rule on the local Pauli basis.

# Arguments
- `sites`: Pauli-space site indices.
- `bond_dim`: MPO bond dimension required by the automaton state space.
- `transition_rule`: Function `(state, local_state) -> (next_state, coeff)`.

# Returns
- A diagonal `MPO` whose action is fully determined by the supplied automaton.
"""
function _diagonal_projector_mpo(sites, bond_dim::Int, transition_rule::Function)
  nsites = length(sites)
  nsites < 1 && throw(ArgumentError("projector MPO requires at least one site"))

  linkinds = [Index(bond_dim, "OperatorProjectorLink,n=$(n)") for n in 1:(nsites - 1)]
  tensors = ITensor[]

  for n in 1:nsites
    site = sites[n]
    left = n > 1 ? linkinds[n - 1] : nothing
    right = n < nsites ? linkinds[n] : nothing
    tensor = if isnothing(left) && isnothing(right)
      ITensor(prime(site), site)
    elseif isnothing(left)
      ITensor(right, prime(site), site)
    elseif isnothing(right)
      ITensor(left, prime(site), site)
    else
      ITensor(left, right, prime(site), site)
    end

    if nsites == 1
      for local_state in 1:dim(site)
        _, coeff = transition_rule(1, local_state)
        tensor[prime(site) => local_state, site => local_state] = coeff
      end
    elseif n == 1
      for local_state in 1:dim(site)
        next_state, coeff = transition_rule(1, local_state)
        tensor[right => next_state, prime(site) => local_state, site => local_state] = coeff
      end
    elseif n == nsites
      for state in 1:bond_dim
        for local_state in 1:dim(site)
          _, coeff = transition_rule(state, local_state)
          tensor[left => state, prime(site) => local_state, site => local_state] += coeff
        end
      end
    else
      for state in 1:bond_dim
        for local_state in 1:dim(site)
          next_state, coeff = transition_rule(state, local_state)
          tensor[left => state, right => next_state, prime(site) => local_state, site => local_state] += coeff
        end
      end
    end

    push!(tensors, tensor)
  end

  return MPO(tensors)
end

"""
    pauli_daoe_projector(sites; lstar, gamma)

Construct the DAOE projector MPO that damps Pauli strings by their Pauli weight beyond cutoff `lstar`.
The local operator basis is assumed to be ordered as `(I, X, Y, Z)`.

# Arguments
- `sites`: Pauli-space site indices.

# Keyword Arguments
- `lstar`: Pauli-weight cutoff. Strings heavier than `lstar` acquire exponential damping.
- `gamma`: Damping strength.

# Returns
- A diagonal `MPO` implementing the DAOE projector.
"""
function pauli_daoe_projector(sites; lstar::Integer, gamma::Real)
  lstar < 0 && throw(ArgumentError("lstar must be nonnegative"))
  bond_dim = Int(lstar) + 1
  transition_rule = (state, local_state) -> begin
    next_weight, coeff = _daoe_transition(state - 1, local_state, Int(lstar), gamma)
    return next_weight + 1, coeff
  end
  return _diagonal_projector_mpo(sites, bond_dim, transition_rule)
end

"""
    fdaoe_projector(sites; wstar, gamma)

Construct the FDAOE projector MPO that damps Pauli-basis operator strings by fermionic Majorana weight beyond cutoff `wstar`.
The local operator basis is assumed to be ordered as `(I, X, Y, Z)`.

# Arguments
- `sites`: Pauli-space site indices.

# Keyword Arguments
- `wstar`: Fermionic Majorana-weight cutoff.
- `gamma`: Damping strength applied beyond the cutoff.

# Returns
- A diagonal `MPO` implementing the fermionic DAOE projector.
"""
function fdaoe_projector(sites; wstar::Integer, gamma::Real)
  wstar < 0 && throw(ArgumentError("wstar must be nonnegative"))
  bond_dim = Int(wstar) + 2
  transition_rule = (state, local_state) -> _fdaoe_transition(state, local_state, Int(wstar), gamma)
  return _diagonal_projector_mpo(sites, bond_dim, transition_rule)
end
