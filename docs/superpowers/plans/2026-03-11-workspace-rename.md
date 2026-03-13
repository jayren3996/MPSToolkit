# Workspace Rename Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rename the workspace directory to `MPSToolkit` and rewrite internal path references so the repo remains self-consistent.

**Architecture:** The change is operational rather than algorithmic: move the repo once, then do a global path rewrite constrained to workspace path strings. Finish with a stale-reference search and a full package test run from the new root.

**Tech Stack:** shell move operations, ripgrep, Perl in-place replacement, Julia package tests.

---

## Chunk 1: Move the workspace

### Task 1: Rename the working directory

**Files:**
- Verify: previous workspace root
- Create: `/Users/ren/Codex/MPSToolkit`

- [ ] **Step 1: Confirm the target path does not already exist**

Run: `test ! -e /Users/ren/Codex/MPSToolkit`
Expected: PASS.

- [ ] **Step 2: Rename the directory**

Run: `mv <old-workspace-root> /Users/ren/Codex/MPSToolkit`
Expected: PASS.

- [ ] **Step 3: Verify the rename**

Run: `test -d /Users/ren/Codex/MPSToolkit && test ! -e <old-workspace-root>`
Expected: PASS.

## Chunk 2: Rewrite repo path references

### Task 2: Update absolute and home-relative workspace paths

**Files:**
- Modify: markdown and notebook files under `/Users/ren/Codex/MPSToolkit`

- [ ] **Step 1: Find stale workspace-path references**

Run: `rg -n "<old-workspace-fragment>" /Users/ren/Codex/MPSToolkit`
Expected: a finite list of docs/specs/plans/notebook matches.

- [ ] **Step 2: Rewrite the path strings**

Use a constrained in-place replacement over matched files only.

- [ ] **Step 3: Re-run the stale-reference search**

Run: `rg -n "<old-workspace-fragment>" /Users/ren/Codex/MPSToolkit`
Expected: no matches.

## Chunk 3: Verify from the new workspace

### Task 3: Run end-to-end verification

**Files:**
- Test: `/Users/ren/Codex/MPSToolkit`

- [ ] **Step 1: Confirm key docs now point at the new root**

Run: `sed -n '1,220p' /Users/ren/Codex/MPSToolkit/README.md`
Expected: links point at `/Users/ren/Codex/MPSToolkit`.

- [ ] **Step 2: Run the full package test suite**

Run: `cd /Users/ren/Codex/MPSToolkit && julia --project=. -e 'using Pkg; Pkg.test()'`
Expected: PASS.
