# Operator-Space DMT Design

**Date:** 2026-03-12

## Goal

Integrate density matrix truncation (DMT) into `MPSToolkit` as an operator-space truncation backend for Pauli-basis `MPS`, starting from faithful 2-site and 3-site implementations based on the legacy `dev/ITensor` code, then generalizing to arbitrary contiguous gate span once correctness is established.

## Scope

In scope:

- operator-space `MPS` in the Pauli basis
- faithful migration of legacy 2-site and 3-site DMT workflows
- a package-quality public API under `MPSToolkit.OperatorSpace`
- tests comparing specialized and generalized implementations
- development-stage benchmarks on a concrete physical model

Out of scope for the first iteration:

- non-Pauli local bases
- infinite systems
- folding DMT into a single opaque evolution routine
- removing the specialized kernels before the generic implementation is validated

## Background

The current package already has:

- dense-gate TEBD evolution for finite `MPS`
- operator-space helpers for Pauli-basis states and gates
- DAOE/FDAOE projector utilities

The legacy DMT code in `dev/ITensor` is operator-space oriented, but it is still script-like and closely tied to specific transport examples. The integration should preserve the useful algorithmic content while broadening the abstraction so DMT is available for generic operator-space dynamics, with transport treated as one benchmark application rather than the defining use case.

## Design Summary

DMT will be integrated as an operator-space truncation backend, not as a new top-level evolution engine.

The flow will remain explicit:

1. apply a local operator-space gate
2. form the active local block
3. build the DMT-preserved subspace from local environments
4. truncate in that constrained basis
5. write tensors back in canonical gauge

This preserves the package's existing separation between evolution and truncation/projection logic.

## Public API

The first iteration should expose a compact API under `MPSToolkit.OperatorSpace`, with names finalized during implementation:

- `dmt2!` for faithful 2-site operator-space DMT
- `dmt3!` for faithful 3-site operator-space DMT
- `dmt_truncate!` or `dmt!` as the eventual generalized contiguous-window entrypoint
- small config/type helpers for DMT options such as `maxdim`, `cutoff`, and preservation-rule settings

The API should remain operator-space focused rather than transport specific.

## Internal Structure

Add `src/operator_space/dmt.jl` and keep responsibilities separated:

- DMT config/types
- left/right environment cache
- faithful 2-site kernel
- faithful 3-site kernel
- generalized contiguous-window kernel
- shared constrained-truncation helpers

Update `src/MPSToolkit.jl` to include and export the new public entrypoints from the `OperatorSpace` namespace.

## Algorithm Staging

### Phase 1: Faithful Specialized Kernels

Port the old 2-site and 3-site DMT logic into package-quality code for Pauli-basis operator-space `MPS`.

Requirements:

- match legacy behavior as closely as possible
- factor out reusable helpers instead of copying large script blocks verbatim
- preserve gauge handling and environment updates carefully

### Phase 2: Generalized Windowed DMT

Design a generalized DMT routine for arbitrary contiguous gate span.

Requirements:

- cover the 2-site and 3-site cases through the same abstraction
- preserve the same local-information constraints as the specialized kernels
- stay readable enough that the generalization does not make maintenance worse

The 2-site and 3-site implementations should remain available until the generalized version is numerically validated against them.

## Preservation Model

The initial preservation model should follow the classic DMT spirit:

- preserve selected low-order local operator information near the active cut/window

The first implementation does not need a fully open-ended user-supplied preservation-rule system unless that abstraction emerges naturally during implementation. However, the internal design should leave room for future custom preservation rules.

## Testing Strategy

Add `test/test_dmt.jl` with three levels of checks:

1. helper-level tests
2. regression tests against legacy-style small examples
3. generalized-vs-specialized equivalence tests

Specific targets:

- environment builder correctness
- constrained matrix truncation helper behavior
- `dmt2!` and `dmt3!` stability on small Pauli-basis operator-space examples
- generalized routine reproducing specialized 2-site and 3-site results to numerical tolerance

## Benchmark Strategy

Keep early benchmarks in `dev/ITensor` or another development-oriented location until the API settles.

Benchmark requirements:

- choose one concrete physical model, such as XXZ or TFIM operator dynamics
- compare the generalized contiguous-window DMT routine against the faithful specialized `dmt2!` and `dmt3!` implementations
- verify that the generalized routine recovers the specialized results for 2-site and 3-site windows to numerical tolerance
- track local observables, profiles, and any relevant truncation diagnostics needed to confirm equivalence

Transport is a valid benchmark application, but the benchmark's main purpose is to validate that the generalized implementation faithfully reproduces the specialized cases.

## Success Criteria

The integration is considered ready when:

- 2-site and 3-site DMT are available in the package for operator-space workflows
- tests confirm the faithful specialized kernels work on small cases
- the generalized routine reproduces specialized results for 2-site and 3-site windows
- at least one benchmark confirms the generalized DMT routine reproduces the specialized `dmt2!` and `dmt3!` behavior on a concrete model

The specialized kernels can be removed only after:

- numerical equivalence with the generalized routine is established
- benchmark quality is not degraded
- the generalized implementation remains understandable and maintainable

## Risks

- the generic contiguous-window formulation may be harder to stabilize than the specialized legacy cases
- subtle gauge/environment mistakes can produce plausible but wrong results
- over-generalizing too early could delay a usable package integration

These risks are why the staged migration is preferred over a generic-first rewrite.

## Open Notes

- The current workspace does not appear to be in a `.git` checkout, so the spec can be written locally but may not be committable from this directory unless a repo root is restored or clarified.
