# MPSToolkit Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Hard rename the package to `MPSToolkit`, reorganize the code into feature modules, and generalize ScarFinder refinement from entropy-only selection to a generic scoring interface with entropy and fidelity selectors.

**Architecture:** Keep the existing finite-MPS backend logic, but move it under an umbrella module with feature-oriented submodules. Introduce a generic selector context and scoring dispatch so ScarFinder can support entropy, fidelity, and future selection rules without changing the loop structure.

**Tech Stack:** Julia, `ITensors.jl`, `ITensorMPS.jl`, Julia `Test`

---

## Chunk 1: Selector Generalization

### Task 1: Add failing tests for generic selector scoring and fidelity selection

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_core.jl`
- Modify: `/Users/ren/Codex/MPSToolkit/test/test_docstrings.jl`

- [ ] **Step 1: Add a failing refinement test that requires selector context**
- [ ] **Step 2: Run the targeted tests and confirm failure**
- [ ] **Step 3: Implement `SelectionContext`, `FidelitySelector`, and generic `score(selector, psi, context)`**
- [ ] **Step 4: Re-run targeted tests**

## Chunk 2: Hard Rename And Module Layout

### Task 2: Rename package/module to `MPSToolkit`

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/Project.toml`
- Create: `/Users/ren/Codex/MPSToolkit/src/MPSToolkit.jl`
- Delete: `/Users/ren/Codex/MPSToolkit/src/ScarFinder.jl`
- Modify: tests, examples, docs, notebooks to use `MPSToolkit`

- [ ] **Step 1: Add failing load tests for `MPSToolkit`**
- [ ] **Step 2: Rename the root module and update imports**
- [ ] **Step 3: Re-run the test suite**

### Task 3: Introduce feature submodules

**Files:**
- Modify: source files under `/Users/ren/Codex/MPSToolkit/src`

- [ ] **Step 1: Reorganize the root module into `Evolution`, `Observables`, `ScarFinder`, `Bases`, and `OperatorSpace` submodules without changing backend behavior**
- [ ] **Step 2: Add Pauli-basis helpers and placeholder `OperatorSpace` namespace**
- [ ] **Step 3: Re-run the test suite**

## Chunk 3: Docs And Examples

### Task 4: Update docs and examples

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/README.md`
- Modify: docs and examples under `/Users/ren/Codex/MPSToolkit/docs` and `/Users/ren/Codex/MPSToolkit/examples`

- [ ] **Step 1: Update docs to describe `MPSToolkit` scope and module layout**
- [ ] **Step 2: Update examples and notebooks to import `MPSToolkit`**
- [ ] **Step 3: Verify notebooks are valid JSON**

## Chunk 4: Verification

### Task 5: Run full verification

**Files:**
- Test: `/Users/ren/Codex/MPSToolkit/test/runtests.jl`

- [ ] **Step 1: Run `julia --project=. -e 'using Pkg; Pkg.test()'`**
- [ ] **Step 2: Smoke test the TDVP and TEBD PXP examples**
