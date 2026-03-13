# Helper API Notebooks Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one physical-spin TEBD helper notebook and one operator-space TEBD helper notebook, then register both in the examples documentation.

**Architecture:** Keep the examples as tutorial notebooks under their existing workflow folders. The TEBD notebook demonstrates the helper stack for physical-spin evolution. The operator-space notebook mirrors the same progression but highlights `map_hamiltonian=pauli_gate_from_hamiltonian`.

**Tech Stack:** Jupyter notebook JSON, Julia, existing MPSToolkit helper APIs

---

## Chunk 1: Create Notebook Content

### Task 1: Add the physical-spin TEBD helper notebook

**Files:**
- Create: `examples/tebd/tebd_helper_apis.ipynb`

- [ ] **Step 1: Build notebook sections**
- helper overview
- `local_gates_from_hamiltonians`
- `tebd_evolution_from_hamiltonians`
- `tebd_strang_schedule`
- `tebd_strang_evolution`
- manual-vs-helper comparison

- [ ] **Step 2: Verify notebook JSON loads**

Run: `python - <<'PY' ...`
Expected: JSON parses and contains the expected markdown/code cells.

### Task 2: Add the operator-space TEBD helper notebook

**Files:**
- Create: `examples/operator_space/operator_tebd_helper_apis.ipynb`

- [ ] **Step 1: Build notebook sections**
- operator-space helper overview
- `local_gates_from_hamiltonians` with `map_hamiltonian`
- `tebd_evolution_from_hamiltonians`
- `tebd_strang_schedule`
- `tebd_strang_evolution`
- manual-vs-helper comparison

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

Run: `rg -n "tebd_helper_apis|operator_tebd_helper_apis" docs/examples.md`
Expected: both notebook names appear in the run instructions and index.

## Chunk 3: Final Verification

### Task 4: Run lightweight verification

**Files:**
- Create/Modify: none

- [ ] **Step 1: Parse both notebook files**

Run a JSON parse check over both notebooks.

- [ ] **Step 2: Execute notebook code cells outside Jupyter**

Extract code cells into temporary Julia scripts and run them with `julia --project`.

- [ ] **Step 3: Summarize residual limitations**

Call out that the notebooks are tutorial examples rather than automated tests and rely on an interactive notebook environment for the full reading experience.
