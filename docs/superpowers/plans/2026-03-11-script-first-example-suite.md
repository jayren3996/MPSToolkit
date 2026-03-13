# Script-First Example Suite Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the old example tree with a compact script-first suite focused on operator space, random-circuit TEBD, and the XYZ ScarFinder workflow.

**Architecture:** Let the new examples drive a final helper layer into the package: generic Strang TEBD builders, spin-1/2 TFIM / XYZ model helpers, and an operator-space dynamics runner. Then rebuild the example tree around short scripts that consume those helpers directly, and update the public docs to match.

**Tech Stack:** Julia, ITensors.jl, ITensorMPS.jl, MPSToolkit.jl, markdown docs.

---

## Chunk 1: Add helper coverage

### Task 1: Test the new helper layer

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_finite_tebd.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_operator_space.jl`
- Create: `/Users/ren/Codex/MPSToolkit/test/test_models.jl`

- [ ] **Step 1: Add failing tests for Strang TEBD helpers**
- [ ] **Step 2: Add failing tests for operator-space dynamics**
- [ ] **Step 3: Add failing tests for spin-1/2 model helpers**
- [ ] **Step 4: Run the targeted tests and verify failure**

## Chunk 2: Implement the helper layer

### Task 2: Add reusable TEBD, model, and operator-dynamics helpers

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/src/evolution/tebd.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/src/operator_space/helpers.jl`
- Create: `/Users/ren/Codex/MPSToolkit/src/models/spinhalf.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/src/MPSToolkit.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/docs/api.md`
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_core.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_docstrings.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/test/runtests.jl`

- [ ] **Step 1: Implement and export the helpers**
- [ ] **Step 2: Update docs and namespace coverage**
- [ ] **Step 3: Run the focused helper tests and verify pass**

## Chunk 3: Replace the example tree

### Task 3: Create the new script-first suite

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/examples/*.jl`
- Delete: overlapping legacy example scripts under `/Users/ren/Codex/MPSToolkit/examples`

- [ ] **Step 1: Add the operator-space scripts**
- [ ] **Step 2: Add the random-circuit TEBD script**
- [ ] **Step 3: Add the XYZ ScarFinder script**
- [ ] **Step 4: Remove legacy overlapping scripts**

## Chunk 4: Update public docs

### Task 4: Rewrite the example index

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/README.md`
- Modify: `/Users/ren/Codex/MPSToolkit/docs/examples.md`

- [ ] **Step 1: Replace notebook-first example references with the new script set**
- [ ] **Step 2: Add the new helper categories to the README/API narrative**

## Chunk 5: Verify representative runs

### Task 5: Run focused verification

**Files:**
- Test: `/Users/ren/Codex/MPSToolkit/examples`
- Test: `/Users/ren/Codex/MPSToolkit`

- [ ] **Step 1: Run the focused test suite covering helpers and namespace docs**
- [ ] **Step 2: Run representative example scripts from each workflow family**
- [ ] **Step 3: Verify package load**
