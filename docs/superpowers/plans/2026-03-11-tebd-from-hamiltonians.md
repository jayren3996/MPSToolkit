# TEBD From Hamiltonians Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reusable convenience helpers that construct `LocalGateEvolution` from dense local Hamiltonians and use them to simplify the TFIM operator notebook.

**Architecture:** Put the Hamiltonian-to-gate conversion helpers in the evolution layer so they are available to both physical and operator-space workflows. Keep the conversion generic by exposing a `map_hamiltonian` hook, then update the TFIM notebook to use `tebd_evolution_from_hamiltonians` with `pauli_gate_from_hamiltonian`.

**Tech Stack:** Julia, LinearAlgebra, ITensors.jl, ITensorMPS.jl, JSON notebook format.

---

## Chunk 1: Add failing tests

### Task 1: Cover the new TEBD helper API

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_finite_tebd.jl`

- [ ] **Step 1: Add a failing test for `local_gates_from_hamiltonians` default exponentiation**
- [ ] **Step 2: Add a failing test for custom `map_hamiltonian`**
- [ ] **Step 3: Add a failing test for `tebd_evolution_from_hamiltonians`**
- [ ] **Step 4: Run the targeted TEBD test file and verify failure**

## Chunk 2: Implement the helper API

### Task 2: Add the evolution helpers

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/src/evolution/tebd.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/src/MPSToolkit.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/docs/api.md`
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_core.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_docstrings.jl`

- [ ] **Step 1: Implement `local_gates_from_hamiltonians`**
- [ ] **Step 2: Implement `tebd_evolution_from_hamiltonians`**
- [ ] **Step 3: Export and document the helpers**
- [ ] **Step 4: Run the targeted TEBD test file and verify pass**

## Chunk 3: Simplify the TFIM notebook

### Task 3: Use the new helper in the dev notebook

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/dev/tfim_operator_tebd.ipynb`

- [ ] **Step 1: Replace manual gate construction with `tebd_evolution_from_hamiltonians`**
- [ ] **Step 2: Keep the notebook self-contained and short**

## Chunk 4: Verify

### Task 4: Run focused verification

**Files:**
- Test: `/Users/ren/Codex/MPSToolkit/test/test_finite_tebd.jl`
- Test: `/Users/ren/Codex/MPSToolkit/dev/tfim_operator_tebd.ipynb`

- [ ] **Step 1: Run `julia --project=. -e 'using Test, MPSToolkit, ITensors, ITensorMPS, LinearAlgebra; include(\"test/test_finite_tebd.jl\"); include(\"test/test_core.jl\"); include(\"test/test_docstrings.jl\")'`**
- [ ] **Step 2: Run `julia --project=. -e 'using MPSToolkit, ITensors, ITensorMPS, LinearAlgebra; println(\"ok\")'`**
- [ ] **Step 3: Run `python3 -c 'import json; json.load(open(\"dev/tfim_operator_tebd.ipynb\"))'`**
