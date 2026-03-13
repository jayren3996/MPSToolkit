# XYZ ScarFinder spiral benchmark
#
# What this shows:
# - The XYZ ring Hamiltonian is built inline with `OpSum`, including the boundary bond
#   between the last and first sites.
# - `tebd_evolution_from_hamiltonians(hamiltonians, dt; schedule, maxdim, cutoff)` builds a
#   TEBD evolution from dense local bond Hamiltonians. Here `schedule` includes `length(sites)`
#   to apply a two-site gate across the periodic boundary.
# - `BondDimTruncation(1; cutoff)` makes the ScarFinder projection aggressively product-state
#   preserving, which is the regime where this spiral benchmark shows attraction.
# - `scarfinder!(psi, evolution, truncation, niter; refine=false)` repeatedly applies the
#   TEBD step and the explicit projection.
#
# The exact scar family for these couplings has unit cell `L = 5`. The script starts from a
# different 5-site product texture, runs ScarFinder, and reports the best overlap with the 5
# translated exact spiral states before and after the projected TEBD loop.

using ITensors
using ITensorMPS
using MPSToolkit

function bloch_product_state(sites, unitcell_vectors)
  tensors = ITensor[]
  for (j, site) in enumerate(sites)
    x, y, z = unitcell_vectors[mod1(j, length(unitcell_vectors))]
    theta = acos(clamp(z, -1.0, 1.0))
    phi = atan(y, x)
    tensor = ITensor(ComplexF64, site)
    tensor[site => 1] = cos(theta / 2)
    tensor[site => 2] = exp(1im * phi) * sin(theta / 2)
    push!(tensors, tensor)
  end
  return MPS(tensors)
end

function translated_unitcell(vectors, shift)
  return [vectors[mod1(j + shift, length(vectors))] for j in 1:length(vectors)]
end

function periodic_xyz_mpo(sites; Jx::Real, Jy::Real, Jz::Real)
  opsum = OpSum()
  for j in 1:(length(sites) - 1)
    opsum += Jx, "Sx", j, "Sx", j + 1
    opsum += Jy, "Sy", j, "Sy", j + 1
    opsum += Jz, "Sz", j, "Sz", j + 1
  end
  opsum += Jx, "Sx", length(sites), "Sx", 1
  opsum += Jy, "Sy", length(sites), "Sy", 1
  opsum += Jz, "Sz", length(sites), "Sz", 1
  return MPO(opsum, sites)
end

Jx = 0.8780104503990049
Jy = 1.0
Jz = 0.28915290771398805

# This exact 5-site unit cell is the commensurate `L = 5` spiral for the couplings above.
exact_unitcell = [
  (0.23132232617119042, 0.8179028885211242, 0.5268062702394030),
  (-0.6306767299128184, 0.5256546655019309, 0.5709063276777785),
  (-0.6306767299128180, -0.5256546655019314, 0.5709063276777785),
  (0.23132232617119097, -0.8179028885211240, 0.5268062702394030),
  (0.8000000000000000, 0.0, 0.6000000000000000),
]

# This seed has the same unit-cell length, but it is not the exact spiral.
seed_unitcell = [
  (-0.19252040748341362, 0.9201628467561563, 0.3409343458087303),
  (-0.20606685230234045, 0.00038614392239448976, 0.9785378394702356),
  (0.20524657952707037, -0.5516411477238218, 0.8084342185548449),
  (0.6812069132317764, -0.07000713019667017, 0.7287359899763847),
  (0.6474235949791325, -0.2840403261925360, 0.7072225828978695),
]

sites = siteinds("S=1/2", 10)
psi = bloch_product_state(sites, seed_unitcell)
spiral_family = [bloch_product_state(sites, translated_unitcell(exact_unitcell, shift)) for shift in 0:4]
hamiltonian = periodic_xyz_mpo(sites; Jx=Jx, Jy=Jy, Jz=Jz)

# `schedule` is the periodic odd-even-odd TEBD order on the 10-site ring.
# The entry `10` denotes the boundary bond `(10, 1)`.
schedule = [1, 3, 5, 7, 9, 2, 4, 6, 8, 10, 1, 3, 5, 7, 9]
# `weights` are the Strang prefactors: half steps on the outer layers and full steps in the middle.
weights = [fill(0.5, 5); fill(1.0, 5); fill(0.5, 5)]
# Each dense local Hamiltonian is the XYZ bond term scaled by the corresponding Strang weight.
bond_hamiltonians = [weight * spinhalf_xyz_bond_hamiltonian(; Jx=Jx, Jy=Jy, Jz=Jz) for weight in weights]
# `dt=0.1` is the TEBD time step for one full sweep.
# `schedule` gives the left bond index for each gate application.
# `maxdim=1` and `cutoff=1e-12` keep the evolution in a strongly projected product-state regime.
evolution = tebd_evolution_from_hamiltonians(
  bond_hamiltonians,
  0.1;
  schedule=schedule,
  maxdim=1,
  cutoff=1e-12,
)
# `maxdim=1` makes the explicit ScarFinder projection rank-1 on every bond.
# `cutoff` discards singular values below the threshold after each projection step.
truncation = BondDimTruncation(1; cutoff=1e-12)

initial_best_overlap = maximum(abs2(inner(reference, psi)) for reference in spiral_family)
initial_energy = energy_density(psi, hamiltonian)

# `niter=150` applies 150 projected TEBD sweeps.
# `refine=false` keeps the benchmark focused on the core ScarFinder loop itself.
scarfinder!(psi, evolution, truncation, 150; refine=false)

final_best_overlap = maximum(abs2(inner(reference, psi)) for reference in spiral_family)
final_energy = energy_density(psi, hamiltonian)

println("initial best spiral overlap = ", initial_best_overlap)
println("final best spiral overlap = ", final_best_overlap)
println("initial energy density = ", initial_energy)
println("final energy density = ", final_energy)
println("final bond entropy = ", bond_entropy(psi, 5))
