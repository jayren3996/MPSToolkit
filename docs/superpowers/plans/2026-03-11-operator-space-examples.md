# Operator-Space Examples Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a dedicated DAOE/FDAOE script and notebook that explain the operator-space projector API in detail and wire them into the public docs.

**Architecture:** Keep the implementation lightweight: one shared example script plus one explanatory notebook, both built around the existing `pauli_siteinds`, `pauli_daoe_projector`, and `fdaoe_projector` APIs. Reuse the same Pauli-basis strings and overlap checks in both assets so the docs and executable example stay aligned.

**Tech Stack:** Julia, ITensors.jl, ITensorMPS.jl, JSON notebook format, existing MPSToolkit docs/examples.

---

## Chunk 1: Add the executable example

### Task 1: Create the DAOE/FDAOE example script

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/examples/operator_space_daoe.jl`
- Test: `/Users/ren/Codex/MPSToolkit/examples/operator_space_daoe.jl`

- [ ] **Step 1: Write the failing test expectation**

Use the script as the behavior target: it should print DAOE/FDAOE overlaps for preserved and damped strings and exit successfully.

- [ ] **Step 2: Run the missing script to verify failure**

Run: `julia --project=. examples/operator_space_daoe.jl`
Expected: FAIL with “No such file” or equivalent because the script does not exist yet.

- [ ] **Step 3: Write the minimal implementation**

Create a script that:
- defines a small helper to build Pauli-basis product `MPS` states
- builds `pauli_siteinds(4)`
- constructs `pauli_daoe_projector(...; lstar=1, gamma=0.3)` and `fdaoe_projector(...; wstar=2, gamma=0.3)`
- applies both projector MPOs to the approved example strings
- prints the measured overlaps next to the expected values

- [ ] **Step 4: Run the script to verify it passes**

Run: `julia --project=. examples/operator_space_daoe.jl`
Expected: PASS with printed overlap summaries for both DAOE and FDAOE.

## Chunk 2: Add the explanatory notebook

### Task 2: Create the notebook tutorial

**Files:**
- Create: `/Users/ren/Codex/MPSToolkit/examples/notebooks/06_operator_space_daoe.ipynb`

- [ ] **Step 1: Define the notebook content**

Include sections for:
- what operator-space sites are
- the `(I, X, Y, Z)` local basis
- DAOE projector construction
- FDAOE projector construction
- why specific strings are preserved or damped
- how `apply` is used

- [ ] **Step 2: Write the notebook JSON**

Build a notebook with markdown explanation cells and runnable code cells that mirror the script.

- [ ] **Step 3: Verify the notebook is valid JSON**

Run: `python3 -c 'import json; json.load(open("examples/notebooks/06_operator_space_daoe.ipynb"))'`
Expected: PASS with no output.

## Chunk 3: Wire the new example into docs

### Task 3: Update the curated example index

**Files:**
- Modify: `/Users/ren/Codex/MPSToolkit/README.md`
- Modify: `/Users/ren/Codex/MPSToolkit/docs/examples.md`

- [ ] **Step 1: Update the run list and notebook order**

Add the new script and notebook alongside the existing curated set.

- [ ] **Step 2: Remove outdated operator-space placeholder wording**

Make the docs describe operator-space projection as supported functionality.

- [ ] **Step 3: Verify the example docs read cleanly**

Run: `sed -n '1,220p' README.md` and `sed -n '1,220p' docs/examples.md`
Expected: the operator-space example is present and placeholder wording is gone.

## Chunk 4: Full verification

### Task 4: Run end-to-end verification

**Files:**
- Test: `/Users/ren/Codex/MPSToolkit/examples/operator_space_daoe.jl`
- Test: `/Users/ren/Codex/MPSToolkit/examples/notebooks/06_operator_space_daoe.ipynb`
- Test: `/Users/ren/Codex/MPSToolkit`

- [ ] **Step 1: Run the new script**

Run: `julia --project=. examples/operator_space_daoe.jl`
Expected: PASS.

- [ ] **Step 2: Parse the notebook JSON**

Run: `python3 -c 'import json; json.load(open("examples/notebooks/06_operator_space_daoe.ipynb"))'`
Expected: PASS.

- [ ] **Step 3: Run the full test suite**

Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS.
