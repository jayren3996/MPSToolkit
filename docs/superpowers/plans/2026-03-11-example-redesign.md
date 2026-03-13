# MPSToolkit Example Redesign Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the existing examples and notebooks with a curated set that teaches the full public MPSToolkit API in a notebook-first format.

**Architecture:** Remove the current mixed example tree and replace it with five focused scripts plus six structured notebooks. Keep model helpers local and lightweight, and ensure the docs point only to the new set. Validate by running the scripts and parsing notebook JSON.

**Tech Stack:** Julia, ITensors.jl, ITensorMPS.jl, Jupyter notebook JSON, MPSToolkit.jl

---

## Chunk 1: Replace the example script set

### Task 1: Remove obsolete example files

**Files:**
- Delete: `/Users/ren/Codex/MPSToolkit/examples/obc_tebd.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/obc_tdvp.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/pxp_common.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/pxp_obc_tdvp.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/pxp_scarfinder_l30.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/pxp_scarfinder_l30_tebd.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/spin1_xy_obc_tdvp.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/spin1_xy_obc_tebd.jl`

- [ ] **Step 1: Delete the obsolete script examples**
- [ ] **Step 2: Confirm only the notebook directory remains under `examples/`**

### Task 2: Add the new curated script set

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/examples/tebd_spin1_xy.jl`
- Create: `/Users/ren/Codex/MPSToolkit/examples/tdvp_spin1_xy.jl`
- Create: `/Users/ren/Codex/MPSToolkit/examples/observables_and_bases.jl`
- Create: `/Users/ren/Codex/MPSToolkit/examples/scarfinder_pxp_tebd.jl`
- Create: `/Users/ren/Codex/MPSToolkit/examples/scarfinder_pxp_tdvp.jl`

- [ ] **Step 1: Write `tebd_spin1_xy.jl` to demonstrate `LocalGateEvolution`, `evolve!`, and TEBD schedule semantics**
- [ ] **Step 2: Write `tdvp_spin1_xy.jl` to demonstrate `TDVPEvolution`, `tdvp_evolve!`, and MPO-based time evolution**
- [ ] **Step 3: Write `observables_and_bases.jl` to demonstrate `energy_density`, `bond_entropy`, `entanglement_spectrum`, `pauli_matrices`, `pauli_basis`, and `pauli_components`**
- [ ] **Step 4: Write `scarfinder_pxp_tebd.jl` to demonstrate TEBD-based ScarFinder with entropy and fidelity selectors**
- [ ] **Step 5: Write `scarfinder_pxp_tdvp.jl` to demonstrate TDVP-based ScarFinder with entropy and fidelity selectors**

## Chunk 2: Replace the notebook set

### Task 3: Remove obsolete notebooks

**Files:**
- Delete: `/Users/ren/Codex/MPSToolkit/examples/notebooks/low_level_usage.ipynb`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/notebooks/pxp_obc_tdvp.ipynb`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/notebooks/pxp_scarfinder_l30.ipynb`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/notebooks/pxp_scarfinder_l30_tebd.ipynb`
- Delete: `/Users/ren/Codex/MPSToolkit/examples/notebooks/spin1_xy_obc_tdvp.ipynb`

- [ ] **Step 1: Delete the obsolete notebooks**

### Task 4: Add the new tutorial notebooks

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/examples/notebooks/00_overview.ipynb`
- Create: `/Users/ren/Codex/MPSToolkit/examples/notebooks/01_evolution_tebd.ipynb`
- Create: `/Users/ren/Codex/MPSToolkit/examples/notebooks/02_evolution_tdvp.ipynb`
- Create: `/Users/ren/Codex/MPSToolkit/examples/notebooks/03_observables_and_bases.ipynb`
- Create: `/Users/ren/Codex/MPSToolkit/examples/notebooks/04_scarfinder_tebd.ipynb`
- Create: `/Users/ren/Codex/MPSToolkit/examples/notebooks/05_scarfinder_tdvp.ipynb`

- [ ] **Step 1: Create `00_overview.ipynb` as the package map and example index**
- [ ] **Step 2: Create `01_evolution_tebd.ipynb` with detailed TEBD argument and keyword explanations**
- [ ] **Step 3: Create `02_evolution_tdvp.ipynb` with detailed TDVP argument and keyword explanations**
- [ ] **Step 4: Create `03_observables_and_bases.ipynb` covering observables, spectrum interpretation, and Pauli-basis helpers**
- [ ] **Step 5: Create `04_scarfinder_tebd.ipynb` covering explicit projection, energy targets, entropy selection, and fidelity selection**
- [ ] **Step 6: Create `05_scarfinder_tdvp.ipynb` covering the same ScarFinder API with MPO TDVP**

## Chunk 3: Update docs and verify

### Task 5: Rewrite example documentation

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/README.md`
- Modify: `/Users/ren/Codex/MPSToolkit/docs/examples.md`

- [ ] **Step 1: Replace the old example list in `README.md` with the new curated script and notebook set**
- [ ] **Step 2: Rewrite `docs/examples.md` as a structured guide to the new notebooks and scripts**

### Task 6: Verification

**Files:**
- Verify: `/Users/ren/Codex/MPSToolkit/examples`
- Verify: `/Users/ren/Codex/MPSToolkit/examples/notebooks`

- [ ] **Step 1: Run the new script examples with `julia --project=.`**
- [ ] **Step 2: Parse every notebook file as JSON to confirm valid notebook structure**
- [ ] **Step 3: Run `julia --project=. -e 'using Pkg; Pkg.test()'` to confirm the package still passes**
