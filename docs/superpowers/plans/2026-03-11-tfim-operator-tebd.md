# TFIM Operator TEBD Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `dev/` notebook that evolves total `S^z` under the TFIM with TEBD in the operator-space basis.

**Architecture:** Keep all new logic notebook-local. The notebook will define Pauli-basis helpers, build bond-resolved TFIM superoperator gates, initialize the operator-space `MPS` for total `S^z`, and provide a small runner that advances the state and records dynamics data.

**Tech Stack:** Julia, LinearAlgebra, ITensors.jl, ITensorMPS.jl, MPSToolkit.jl, JSON notebook format.

---

## Chunk 1: Document the approach

### Task 1: Save the approved design

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/docs/superpowers/specs/2026-03-11-tfim-operator-tebd-design.md`

- [ ] **Step 1: Write the design summary**

Record the approved scope, approach, and verification targets for the notebook.

## Chunk 2: Add the notebook

### Task 2: Create the TFIM operator notebook

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/dev/tfim_operator_tebd.ipynb`

- [ ] **Step 1: Define the notebook helpers**

Include notebook-local helpers for:
- Pauli-basis matrices and basis transforms
- operator-space product-state construction
- total-`S^z` operator-state construction
- TFIM bond Hamiltonians and operator-space gates
- odd-even-odd TEBD schedule assembly

- [ ] **Step 2: Add the notebook runner**

Write the TEBD loop explicitly in the example so it evolves the operator and records time, overlap, and operator-entanglement data without hiding that logic in a helper.

- [ ] **Step 3: Add a smoke example cell**

Instantiate a small system and run a short evolution so the notebook is immediately usable.

## Chunk 3: Verify the notebook

### Task 3: Run lightweight verification

**Files:**
- Test: `/Users/ren/Codex/MPSToolkit/dev/tfim_operator_tebd.ipynb`
- Test: `/Users/ren/Codex/MPSToolkit`

- [ ] **Step 1: Parse the notebook JSON**

Run: `python3 -c 'import json; json.load(open("dev/tfim_operator_tebd.ipynb"))'`
Expected: PASS with no output.

- [ ] **Step 2: Verify Julia can load the package**

Run: `julia --project=. -e 'using MPSToolkit, ITensors, ITensorMPS, LinearAlgebra; println("ok")'`
Expected: PASS with `ok`.

- [ ] **Step 3: Run a one-step smoke test**

Run a Julia command that reconstructs the notebook helpers, builds the evolution object, and advances total `S^z` by one TEBD step.
Expected: PASS with a printed scalar diagnostic and exit code `0`.
