# Workspace Rename Design

## Goal

Rename the working folder to `/Users/ren/Codex/MPSToolkit` so the workspace name matches the package name.

## Scope

This rename should include:

- moving the workspace directory itself
- rewriting absolute and home-relative path references inside the repo
- keeping docs, notebooks, specs, and plans internally consistent after the move

This change does not alter package code or public Julia APIs.

## Approach

Use a hard workspace rename:

1. write the spec and plan in the current repo
2. move the directory to `/Users/ren/Codex/MPSToolkit`
3. rewrite repo contents from the old workspace path to the new one
4. verify no stale old-workspace references remain
5. run the full test suite from the new location

## Path Rewrite Rules

Rewrite old workspace-path strings inside the repo so they point at the new root.

The main places affected are:

- `README.md`
- `docs/*.md`
- `docs/superpowers/**/*.md`
- notebook JSON files under `examples/notebooks`

## Verification

Verification should include:

- confirming `/Users/ren/Codex/MPSToolkit` exists
- confirming the previous workspace path no longer exists
- `rg` for stale old-workspace references
- `julia --project=. -e 'using Pkg; Pkg.test()'` from the new directory
