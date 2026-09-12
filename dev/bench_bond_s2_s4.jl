using ITensors
using ITensorMPS
using LinearAlgebra
using MPSToolkit
using Statistics

const N = parse(Int, get(ENV, "MPSTOOLKIT_S24_NSITES", "6"))
const TOTAL_TIME = parse(Float64, get(ENV, "MPSTOOLKIT_S24_TOTAL_TIME", "0.4"))
const STEPS = parse.(Int, split(get(ENV, "MPSTOOLKIT_S24_STEPS", "1,2,4,8"), ','))
const CHIS = parse.(Int, split(get(ENV, "MPSTOOLKIT_S24_CHIS", "9,16"), ','))
const REPEATS = parse(Int, get(ENV, "MPSTOOLKIT_S24_REPEATS", "3"))
const THRESHOLDS = parse.(Float64, split(get(ENV, "MPSTOOLKIT_S24_THRESHOLDS", "1e-3,1e-4,1e-5,1e-6"), ','))
const BLAS_THREADS = parse(Int, get(ENV, "MPSTOOLKIT_S24_BLAS_THREADS", "1"))
BLAS.set_num_threads(BLAS_THREADS)

const I2 = Matrix{ComplexF64}(I, 2, 2)
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]

embed(op, start) = kron(start == 1 ? ones(ComplexF64, 1, 1) : foldl(kron, fill(I2, start - 1)),
  Matrix{ComplexF64}(op),
  N - start - round(Int, log2(size(op, 1))) + 1 == 0 ? ones(ComplexF64, 1, 1) :
    foldl(kron, fill(I2, N - start - round(Int, log2(size(op, 1))) + 1)))

function initial_density(model)
  factors = Matrix{ComplexF64}[]
  for site in 1:N
    sign = site <= N ÷ 2 ? 1.0 : -1.0
    axis = model === :xxz ? Z : (0.7Z + 0.3X) / sqrt(0.58)
    push!(factors, I2 + 0.08sign * axis)
  end
  return factors, foldl(kron, factors)
end

function model_data(model)
  if model === :mfi
    terms = [spinhalf_mixed_field_ising_bond_hamiltonian(N, bond;
      J=0.7, gx=1.1, gz=-0.3) for bond in 1:(N - 1)]
    probes = [(bond, terms[bond]) for bond in 1:(N - 1)]
    transport = (N ÷ 2, terms[N ÷ 2])
  else
    term = spinhalf_xyz_bond_hamiltonian(Jx=1.0, Jy=1.0, Jz=0.6)
    terms = fill(term, N - 1)
    probes = [(site, Z) for site in 1:N]
    transport = (N ÷ 2, kron(X, Y) - kron(Y, X))
  end
  hamiltonian = sum(embed(terms[bond], bond) for bond in 1:(N - 1))
  return terms, probes, transport, hamiltonian
end

function dense_observables(rho, probes, transport, hamiltonian)
  trace_value = tr(rho)
  values = real.([tr(rho * embed(op, start)) / trace_value for (start, op) in probes])
  flow = real(tr(rho * embed(transport[2], transport[1])) / trace_value)
  energy = real(tr(rho * hamiltonian) / trace_value)
  return (profile=values, transport=flow, trace=trace_value, energy=energy)
end

function mps_observables(rho, probes, transport, terms)
  profile = real.(pauli_expectation_profile(rho, probes))
  flow = real(pauli_expectation(rho, transport[2], transport[1]))
  trace_value = pauli_trace(rho)
  energy = sum(real.(pauli_expectation_profile(rho,
    [(bond, terms[bond]) for bond in 1:(N - 1)])))
  return (profile=profile, transport=flow, trace=trace_value, energy=energy)
end

function observable_error(actual, exact)
  profile = maximum(abs.(actual.profile .- exact.profile))
  flow = abs(actual.transport - exact.transport)
  energy = abs(actual.energy - exact.energy)
  trace_error = abs(actual.trace - exact.trace) / max(abs(exact.trace), 1.0)
  return maximum((profile, flow, energy, trace_error)), (profile, flow, energy, trace_error)
end

function active_diagnostics(base, plan)
  state = copy(base)
  active = 0
  for _ in 1:plan.nstep, (entry, gate) in zip(plan.entries, plan.gates)
    MPSToolkit._exact_gate_qr!(state, gate, entry.bond, 2, entry.direction)
    active += dim(linkind(state, entry.bond)) > plan.options.maxdim
    MPSToolkit._dmt_window_truncate!(state, entry.bond, 2; maxdim=plan.options.maxdim,
      cutoff=plan.options.cutoff, direction=entry.direction,
      preserve_diameter=plan.options.preserve_diameter,
      preserve_operators=plan.options.preserve_operators,
      truncation=plan.options.truncation, cache=nothing)
  end
  return active
end

println("# Fixed-T S2/S4 observable-error benchmark")
println("Julia=$(VERSION), CPU=$(Sys.CPU_NAME), BLAS=$(BLAS.get_config()), BLAS_threads=$(BLAS.get_num_threads())")
println("N=$(N), T=$(TOTAL_TIME), steps=$(STEPS), chis=$(CHIS), repeats=$(REPEATS), backend=:qr")
println("Inputs are generated product density operators; exact references use dense exp(-i T H).")
println("model,order,chi,nstep,tau,error,profile_error,transport_error,energy_error,trace_error,wall_s,wall_mad_s,attempted,active")

rows = NamedTuple[]
for model in (:mfi, :xxz)
  terms, probes, transport, hamiltonian = model_data(model)
  factors, rho0 = initial_density(model)
  unitary = exp(-im * TOTAL_TIME * hamiltonian)
  exact = dense_observables(unitary * rho0 * unitary', probes, transport, hamiltonian)
  sites = pauli_siteinds(N)
  base = operator_product_state(sites, factors)
  for order in (2, 4), nstep in STEPS
    tau = TOTAL_TIME / nstep
    plan = bond_dmt_plan(sites, terms, tau; order=order, nstep=1, maxdim=9,
      cutoff=0.0, gate_backend=:qr, normalize=false)
    step_unitary = Matrix{ComplexF64}(I, 2^N, 2^N)
    for entry in plan.entries
      local_unitary = exp(-im * entry.duration * terms[entry.bond])
      step_unitary = embed(local_unitary, entry.bond) * step_unitary
    end
    trotter = step_unitary^nstep
    ed_actual = dense_observables(trotter * rho0 * trotter', probes, transport, hamiltonian)
    error, parts = observable_error(ed_actual, exact)
    println(join((model, order, "ED", nstep, tau, error, parts..., 0.0, 0.0,
      nstep * bond_dmt_attempted_updates(plan), 0), ','))
  end
  for order in (2, 4), chi in CHIS, nstep in STEPS
    tau = TOTAL_TIME / nstep
    plan = bond_dmt_plan(sites, terms, tau; order=order, nstep=nstep, maxdim=chi,
      cutoff=0.0, gate_backend=:qr, normalize=false)
    warm = copy(base); dmt_evolve!(warm, plan)
    times = Float64[]
    result = nothing
    for _ in 1:REPEATS
      state = copy(base)
      GC.gc()
      push!(times, @elapsed dmt_evolve!(state, plan))
      result = state
    end
    actual = mps_observables(result, probes, transport, terms)
    error, parts = observable_error(actual, exact)
    attempted = nstep * bond_dmt_attempted_updates(plan)
    active = active_diagnostics(base, plan)
    row = (model=model, order=order, chi=chi, nstep=nstep, tau=tau, error=error,
      wall=median(times), attempted=attempted, active=active)
    push!(rows, row)
    println(join((model, order, chi, nstep, tau, error, parts..., median(times),
      median(abs.(times .- median(times))), attempted, active), ','))
  end
end

println("# largest acceptable tau (and its measured wall)")
println("model,threshold,order,chi,tau,wall_s,error")
for model in (:mfi, :xxz), threshold in THRESHOLDS, order in (2, 4), chi in CHIS
  candidates = filter(row -> row.model === model && row.order == order && row.chi == chi &&
    row.error <= threshold, rows)
  isempty(candidates) && continue
  chosen = candidates[argmax(getfield.(candidates, :tau))]
  println(join((model, threshold, order, chi, chosen.tau, chosen.wall, chosen.error), ','))
end
