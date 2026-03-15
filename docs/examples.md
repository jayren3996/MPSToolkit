# Examples

The example tree is grouped by workflow so related scripts stay together:

- `examples/benchmarks/`
- `examples/chebyshev/`
- `examples/daoe/`
- `examples/open_systems/`
- `examples/operator_space/`
- `examples/scarfinder/`
- `examples/tebd/`

## Running Examples

From the package root:

```bash
julia --project=. examples/daoe/projectors.jl
julia --project=. examples/chebyshev/two_peak_spectrum.jl
julia --project=. examples/operator_space/tfim_local_z.jl
julia --project=. examples/operator_space/tfim_total_sz.jl
julia --project=. examples/operator_space/tfim_string.jl
julia --project=. examples/operator_space/tfim_autocorrelator.jl
julia --project=. examples/operator_space/tfim_entanglement.jl
julia --project=. examples/operator_space/custom_hamiltonians.jl
julia --project=. examples/operator_space/xyz_local_z.jl
julia --project=. examples/scarfinder/xyz_spiral.jl
```

Open the notebook benchmark in Jupyter or VS Code:

```bash
jupyter lab examples/benchmarks/pbc_tdvp_vs_tebd.ipynb
jupyter lab examples/chebyshev/two_peak_spectrum.ipynb
jupyter lab examples/chebyshev/tfim_local_spectrum.ipynb
jupyter lab examples/chebyshev/energy_cutoff_comparison.ipynb
jupyter lab examples/tebd/xxz_tebd_vs_ed.ipynb
jupyter lab examples/tebd/disordered_xxz_mbl_tebd.ipynb
jupyter lab examples/tebd/scheduler_patterns.ipynb
jupyter lab examples/tebd/tebd_helper_apis.ipynb
jupyter lab examples/open_systems/pauli_lindblad_tebd.ipynb
jupyter lab examples/open_systems/boundary_driven_xxz_steady_state.ipynb
jupyter lab examples/operator_space/dmt_scheduler.ipynb
jupyter lab examples/operator_space/operator_tebd_helper_apis.ipynb
```

## Workflow Index

### Chebyshev

- [two_peak_spectrum.jl](/Users/ren/Codex/MPSToolkit/examples/chebyshev/two_peak_spectrum.jl)
  Minimal Chebyshev spectral reconstruction on a diagonal spin chain, showing how moment generation, Jackson damping, and physical-frequency rescaling fit together.
- [two_peak_spectrum.ipynb](/Users/ren/Codex/MPSToolkit/examples/chebyshev/two_peak_spectrum.ipynb)
  Notebook version of the same workflow, with a method overview, model explanation, and plots of the reconstructed spectrum and moment sequence.
- [tfim_local_spectrum.ipynb](/Users/ren/Codex/MPSToolkit/examples/chebyshev/tfim_local_spectrum.ipynb)
  Local TFIM spectral function from the probe state `S^z_j |gs⟩`, comparing the Chebyshev reconstruction against exact finite-size spectral sticks.
- [energy_cutoff_comparison.ipynb](/Users/ren/Codex/MPSToolkit/examples/chebyshev/energy_cutoff_comparison.ipynb)
  Notebook comparison for the Chebyshev energy-window cutoff, starting from a small exact low-energy reference and then moving to a larger truncation-limited disordered-TFIM calculation.

### DAOE

- [projectors.jl](/Users/ren/Codex/MPSToolkit/examples/daoe/projectors.jl)
  DAOE/FDAOE projector application on explicit Pauli strings.

### Operator Space

- [tfim_local_z.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_local_z.jl)
  Local-operator TFIM evolution with `tebd_strang_evolution`.
- [tfim_total_sz.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_total_sz.jl)
  Extensive conserved-operator workflow.
- [tfim_string.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_string.jl)
  Pauli-string initialization and evolution.
- [tfim_autocorrelator.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_autocorrelator.jl)
  Autocorrelation tracing from an explicit TEBD-plus-measurement loop.
- [tfim_entanglement.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/tfim_entanglement.jl)
  Operator-entanglement growth.
- [custom_hamiltonians.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/custom_hamiltonians.jl)
  Manual local-Hamiltonian matrices through the generic TEBD API.
- [xyz_local_z.jl](/Users/ren/Codex/MPSToolkit/examples/operator_space/xyz_local_z.jl)
  Same operator workflow on a second model using `spinhalf_xyz_bond_hamiltonian`.
- [dmt_scheduler.ipynb](/Users/ren/Codex/MPSToolkit/examples/operator_space/dmt_scheduler.ipynb)
  Notebook walkthrough of scheduled operator-space DMT, showing `DMTGateEvolution`, `dmt_evolve!`, and the low-level `dmt_step!` relation.
- [operator_tebd_helper_apis.ipynb](/Users/ren/Codex/MPSToolkit/examples/operator_space/operator_tebd_helper_apis.ipynb)
  Notebook guide to the packed TEBD helper APIs in operator space, including `map_hamiltonian=pauli_gate_from_hamiltonian`.
- [pauli_lindblad_tebd.ipynb](/Users/ren/Codex/MPSToolkit/examples/open_systems/pauli_lindblad_tebd.ipynb)
  Lindblad TEBD notebook in the Pauli basis, using local dephasing jumps on a disordered XXZ chain.
- [boundary_driven_xxz_steady_state.ipynb](/Users/ren/Codex/MPSToolkit/examples/open_systems/boundary_driven_xxz_steady_state.ipynb)
  Boundary-driven XXZ steady-state notebook in the Pauli basis, comparing easy-plane and easy-axis magnetization profiles.

### TEBD

- [xxz_tebd_vs_ed.ipynb](/Users/ren/Codex/MPSToolkit/examples/tebd/xxz_tebd_vs_ed.ipynb)
  Small-chain TEBD benchmark against exact dense XXZ evolution using `EDKit` for `MPS`-to-vector conversion and exact entanglement.
- [disordered_xxz_mbl_tebd.ipynb](/Users/ren/Codex/MPSToolkit/examples/tebd/disordered_xxz_mbl_tebd.ipynb)
  Disordered-XXZ TEBD notebook comparing weak and strong disorder through imbalance, entanglement growth, and bond-dimension demand.
- [scheduler_patterns.ipynb](/Users/ren/Codex/MPSToolkit/examples/tebd/scheduler_patterns.ipynb)
  Tutorial notebook for scheduler-driven TEBD, covering standard sweeps, 2-site and 3-site brick schedules, and a shallow mixed-span random circuit.
- [tebd_helper_apis.ipynb](/Users/ren/Codex/MPSToolkit/examples/tebd/tebd_helper_apis.ipynb)
  Notebook walkthrough of the packed physical-spin TEBD helpers, from `local_gates_from_hamiltonians` up through `tebd_strang_evolution`.

### Benchmarks

- [pbc_tdvp_vs_tebd.ipynb](/Users/ren/Codex/MPSToolkit/examples/benchmarks/pbc_tdvp_vs_tebd.ipynb)
  Periodic-ring benchmark notebook comparing TDVP, TEBD with an explicit boundary gate, and exact evolution across both high-accuracy and truncation-limited regimes.

### ScarFinder

- [xyz_spiral.jl](/Users/ren/Codex/MPSToolkit/examples/scarfinder/xyz_spiral.jl)
  Periodic-TEBD ScarFinder benchmark against the exact XYZ spiral family.
