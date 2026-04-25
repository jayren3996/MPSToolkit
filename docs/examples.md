# Examples

The maintained example index lives in [`src/examples.md`](src/examples.md), which is used by
Documenter. This legacy root page is kept as a short pointer so direct links into `docs/` do not
land on stale absolute paths.

## Scripts

Run script examples from the package root:

```bash
julia --project=. examples/daoe/projectors.jl
julia --project=. examples/operator_space/tfim_local_z.jl
julia --project=. examples/operator_space/tfim_total_sz.jl
julia --project=. examples/operator_space/tfim_string.jl
julia --project=. examples/operator_space/tfim_autocorrelator.jl
julia --project=. examples/operator_space/tfim_entanglement.jl
julia --project=. examples/operator_space/custom_hamiltonians.jl
julia --project=. examples/operator_space/xyz_local_z.jl
julia --project=. examples/scarfinder/xyz_spiral.jl
```

## Notebooks

Open notebooks in Jupyter or VS Code, for example:

```bash
jupyter lab examples/tdvp/pbc_tdvp_vs_tebd.ipynb
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
jupyter lab examples/scarfinder/pxp_scarfinder.ipynb
```
