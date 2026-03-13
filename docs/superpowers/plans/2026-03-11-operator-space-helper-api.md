# Operator-Space Helper API Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reusable operator-space helper constructors to `MPSToolkit` and simplify the TFIM development notebook to use them.

**Architecture:** Put the reusable Pauli-basis operator-space logic in `src/operator_space/` and expose it through `MPSToolkit.OperatorSpace`. Keep TFIM-specific assembly in the notebook, but remove notebook-local definitions for product-state construction, total `S^z`, and physical-to-operator-space gate conversion.

**Tech Stack:** Julia, ITensors.jl, ITensorMPS.jl, LinearAlgebra, JSON notebook format.

---

## Chunk 1: Add tests for the new helper API

### Task 1: Write failing operator-space helper tests

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_operator_space.jl`

- [ ] **Step 1: Add a test for `pauli_basis_state`**
- [ ] **Step 2: Run the targeted test and verify it fails because the helper is missing**
- [ ] **Step 3: Add a test for `pauli_total_sz_state`**
- [ ] **Step 4: Add a test for `pauli_gate_from_hamiltonian`**

## Chunk 2: Implement the helper API

### Task 2: Add reusable operator-space helpers

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/src/operator_space/helpers.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/src/MPSToolkit.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/docs/api.md`

- [ ] **Step 1: Implement `pauli_basis_state` with default Pauli ordering**
- [ ] **Step 2: Implement `pauli_total_sz_state`**
- [ ] **Step 3: Implement `pauli_gate` and `pauli_gate_from_hamiltonian`**
- [ ] **Step 4: Export the helpers and document them**

## Chunk 3: Simplify the notebook

### Task 3: Rewrite the TFIM notebook around the helper API

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/dev/tfim_operator_tebd.ipynb`

- [ ] **Step 1: Remove notebook-local general helper code**
- [ ] **Step 2: Update the notebook to use the package helpers directly**
- [ ] **Step 3: Keep the visible usage cell short**

## Chunk 4: Verify

### Task 4: Run verification

**Files:**
- Test: `/Users/ren/Codex/MPSToolkit/test/test_operator_space.jl`
- Test: `/Users/ren/Codex/MPSToolkit/dev/tfim_operator_tebd.ipynb`

- [ ] **Step 1: Run `julia --project=. -e 'using Pkg; Pkg.test(; test_args=[\"test_operator_space\"])'`**
- [ ] **Step 2: Run `julia --project=. -e 'using MPSToolkit, ITensors, ITensorMPS, LinearAlgebra; println(\"ok\")'`**
- [ ] **Step 3: Run `python3 -c 'import json; json.load(open(\"dev/tfim_operator_tebd.ipynb\"))'`**
