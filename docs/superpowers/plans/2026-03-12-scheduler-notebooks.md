# Scheduler Notebooks Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one TEBD scheduler tutorial notebook and one operator-space DMT scheduler tutorial notebook, then register them in the examples documentation.

**Architecture:** Keep scheduler teaching in notebooks rather than source files. The TEBD notebook lives in `examples/tebd/` and demonstrates three scheduling patterns in one walkthrough. The DMT notebook lives in `examples/operator_space/` and explains the shared schedule concept with the DMT backend. Documentation updates only need to index and describe the new notebooks.

**Tech Stack:** Jupyter notebook JSON, Julia, existing MPSToolkit example APIs

---

## Chunk 1: Create Notebook Content

### Task 1: Add the TEBD scheduler notebook

**Files:**
- Create: `examples/tebd/scheduler_patterns.ipynb`

- [ ] **Step 1: Build notebook sections**
- intro and scheduler concept
- standard sweep example
- 2-site and 3-site brick schedules
- shallow random mixed-span circuit

- [ ] **Step 2: Verify notebook JSON loads**

Run: `python - <<'PY' ...`
Expected: JSON parses and contains the expected markdown/code cells.

### Task 2: Add the DMT scheduler notebook

**Files:**
- Create: `examples/operator_space/dmt_scheduler.ipynb`

- [ ] **Step 1: Build notebook sections**
- operator-space DMT introduction
- Pauli-basis setup
- scheduled `DMTGateEvolution`
- local observable inspection

- [ ] **Step 2: Verify notebook JSON loads**

Run: `python - <<'PY' ...`
Expected: JSON parses and contains the expected markdown/code cells.

## Chunk 2: Document the New Notebooks

### Task 3: Update examples index

**Files:**
- Modify: `docs/examples.md`

- [ ] **Step 1: Add notebook entries**

Document both notebooks in the run instructions and workflow index.

- [ ] **Step 2: Verify docs references**

Run: `rg -n "scheduler_patterns|dmt_scheduler" docs/examples.md`
Expected: both notebook names appear in the run instructions and index.

## Chunk 3: Final Verification

### Task 4: Run lightweight verification

**Files:**
- Create/Modify: none

- [ ] **Step 1: Parse both notebook files**

Run a JSON parse check over both notebooks.

- [ ] **Step 2: Summarize any residual limitations**

Call out that the notebooks are tutorial examples rather than automated tests, and that they rely on an interactive Jupyter environment for execution.
