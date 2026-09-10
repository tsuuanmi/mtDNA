# AGENTS.md — Agent Development Guide

This file defines the rules every coding agent must follow in this repository. The goals are: understand docs and source before changing code; make small, intentional changes; keep the workspace and commit history clean; avoid overwriting other agents' work; verify changes with the right checks; and keep docs and changelogs synchronized with code.

---

## 1. Instruction Priority

Follow instructions in this order: the user's latest explicit instruction, this `AGENTS.md`, project docs (`README.md`, `CONTRIBUTING.md`, etc.), then existing code patterns. If the user's instruction conflicts with this file, explain the conflict and ask for explicit confirmation before overriding a safety rule.

---

## 2. Communication

Be concise, direct, and technical. No emojis in commits, issues, PR comments, code, or technical summaries. No filler or excessive praise. Answer the user's direct question before making edits or running implementation commands. When responding to feedback or analysis, explicitly say whether you agree or disagree before explaining changes. State assumptions when proceeding under uncertainty.

Ask before proceeding when ambiguity affects correctness, data loss, public APIs, schemas, user-visible behavior, dependency changes, destructive operations, or removal of intentional functionality. For low-risk ambiguity, state the assumption and proceed with the simplest reversible approach.

---

## 3. Docs and Source First

Before implementation, read both relevant documentation and affected source code. Docs explain intent; source code is the ground truth. Read relevant docs when the task touches architecture, public APIs, CLI behavior, configuration, pipeline stages, data flow, schemas, generated outputs, or user-visible behavior. Read affected source files in full before editing them, making broad changes, investigating behavior, auditing correctness, or modifying files you have not fully inspected. Do not rely only on search snippets for broad or sensitive changes. If docs and source disagree, follow the source code and update the docs when documented behavior changes.

---

## 4. Standard Workflow

**Context.** Before editing: read relevant docs and source files, check existing patterns/utilities/naming/error handling, identify all affected files (source, tests, docs, changelogs, generated files, configs, lockfiles), and determine which language-specific rules apply.

**Plan.** Before implementation, briefly state what will change, what will stay the same, likely affected files, and verification to run. Ask for confirmation only when the change is risky, destructive, blocking-ambiguous, or broader than requested.

**Implement.** Work in an isolated worktree when available (`omx -w <branch>`). Create temporary backups in `/tmp` before editing files. Make surgical changes only. Follow existing style, use existing utilities before adding new ones, avoid unrelated cleanup and speculative abstractions, and preserve intentional behavior unless the user asks to change it.

**Verify.** Before finalizing: review diffs, run required language-specific checks, run only allowed tests, confirm no backup files were created inside the repo, update docs and changelogs when required, and confirm `git status` contains only intentional changes.

**Final response** must include: what changed, files modified, verification commands run, tests run or skipped, docs/changelog status, and assumptions or risks. Keep it concise.

---

## 5. Temporary Backups

Always store edit backups in `/tmp/agent-backups/<task-id>/`, mirroring the original relative path to avoid collisions. Never create `.bak` files inside the repository. Never use WIP commits as checkpoints. Commit history must stay clean.

```bash
# Before editing
task_id="<short-task-name>"
file="path/to/file"
backup_root="/tmp/agent-backups/$task_id"
mkdir -p "$backup_root/$(dirname "$file")"
cp "$file" "$backup_root/$file.bak"

# After editing
diff "$backup_root/$file.bak" "$file"

# To restore
cp "$backup_root/$file.bak" "$file"
```

For risky multi-step edits, use step-specific subdirectories like `step-1/`, `step-2/` inside the backup root. Before checks and final response, confirm the repo contains no backup files.

---

## 6. Git Rules

Multiple agent sessions may be running in the same repository — do not touch work that is not yours. Safe commands: `git status`, `git diff`, `git diff -- <path>`, `git diff --stat`, `git add <explicit-path>`, `git restore -- <explicit-path>`. Only restore files you modified in this session.

Never run: `git reset --hard`, `git checkout .`, `git clean -fd`, `git stash`, `git add .`, `git add -A`, `git commit --no-verify`. When staging, use explicit paths only (`git add <path1> <path2>`). Never stage the entire repository. Never commit unless the user explicitly asks.

When committing: run required checks first, run `git status`, stage explicit paths only, run `git diff --cached`, commit only files changed in this session, use a concise informative message (`<type>(<scope>): <message>`), and do not create WIP checkpoint commits. If conflicts occur, resolve only files you modified; if a conflict appears in a file you did not modify, stop and ask; never force push unless the user explicitly confirms the risk.

---

## 7. Change Scope

Every changed line must trace back to the task. Do not refactor unrelated code, reformat unrelated files, rename unrelated symbols, remove unrelated dead code, upgrade unrelated dependencies, modify generated files directly, change lockfiles unless required, or add new configuration unless needed. If unrelated issues are found, mention them instead of fixing them. Prefer the smallest correct change. Do not add speculative features, single-use abstractions, unused configuration, unnecessary flexibility or error handling, or broad rewrites for small fixes. Before adding new logic, check for existing helpers, utilities, models, constants, validators, serializers, logging patterns, error types, CLI patterns, and test patterns.

---

## 8. Language Rules

### Python

Follow existing project style. Prefer standard library first, existing utilities, typed functions, small functions, explicit realistic error handling, and clear data models. When already used by the project, prefer Loguru for logging, Pydantic for data models, `pathlib` for path handling, and `pytest` conventions for tests. Avoid untyped public functions, broad `except Exception`, global mutable state, duplicate validation logic, speculative abstractions, and unrelated formatting changes. Use the project package manager (e.g., `uv` if the project uses it).

```bash
# Dependencies
uv add <package>
uv sync

# Checks
ruff check .
ruff format .
basedpyright

# Scoped checks (project may expect)
ruff check src/ tests/
ruff format src/ tests/
basedpyright src/
```

### TypeScript

Follow existing project style. Prefer top-level imports, explicit public-boundary types, existing utilities and types, small functions, simple control flow, and strict type safety. Avoid `any` unless absolutely necessary, unsafe casts, duplicate type definitions, inline imports (`await import("pkg")` / `type X = import("pkg").X`), dynamic type imports, broad formatting changes, and dependency downgrades to fix type errors. If the project runs TypeScript in strip-only mode, use only erasable TypeScript syntax. Avoid syntax requiring JavaScript emit: `enum`, `namespace`, `module`, parameter properties, `import =`, `export =`. Use explicit fields with constructor assignments. Check package types in `node_modules` or official type definitions — do not guess external API shapes. Treat dependency and lockfile changes as reviewed code.

```bash
# Safe install commands
npm install --ignore-scripts
npm ci --ignore-scripts
npm install --package-lock-only --ignore-scripts

# Checks
npm run check
npm run lint
npm run typecheck
npm run format:check
```

### Rust

Follow existing project style. Prefer standard library first, existing project modules, explicit error types, `Result` for fallible operations, small functions, clear ownership, and minimal cloning. Use project-adopted libraries when already present (`log`/`tracing` for logging, `serde` for serialization, `clap` for CLI, `thiserror` for errors). Avoid unnecessary `clone`, broad rewrites, `unwrap`/`expect` in production unless already accepted nearby, duplicate logic, unrelated formatting changes, and unnecessary features in `Cargo.toml`.

```bash
# Dependencies
cargo add <crate>
cargo update -p <crate>

# Checks
cargo fmt --check
cargo clippy -- -D warnings

# Workspace checks (if project expects)
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
```

---

## 9. Tests

Do not run full test suites unless explicitly asked. Do not create new tests unless explicitly asked. If asked to create or modify a test, run only the relevant test file or case and iterate until it passes. Do not run integration, network, credential-dependent, paid-provider, destructive, or long-running tests unless explicitly requested.

```bash
# Python — targeted only
pytest path/to/test_file.py -q
pytest path/to/test_file.py::test_name -q

# TypeScript — targeted only
npm run test -- path/to/test.ts
node ./node_modules/vitest/vitest.mjs --run path/to/test.ts

# Rust — targeted only
cargo test test_name
cargo test -p crate_name test_name
cargo test --test integration_file test_name
```

---

## 10. Generated Files

Do not modify generated files directly unless project rules explicitly allow it. If a generated file must change, find the generator, modify the generator or source data, regenerate, review the diff, and mention generated files in the final response. Generated files may include `*.generated.ts`, API clients, schemas, protobuf outputs, OpenAPI outputs, model metadata, or codegen snapshots. If unsure whether a file is generated, inspect headers, docs, and build scripts before editing.

---

## 11. Verification

Before running lint, format, type check, tests, builds, or package commands, confirm there are no in-repo backup files:

```bash
git status
find . -name "*.bak" -print
```

If `.bak` files are found, stop — do not remove them without confirmation. Run language-specific checks based on change type: docs lint for docs-only changes; Python lint/format/typecheck for Python source; TypeScript check command for TypeScript source; Rust format/clippy for Rust source; safe package-manager commands and lockfile review for dependency changes; generator + diff review for generated file changes.

Before finalizing, run `git diff` and `git status` and confirm every changed line is intentional, no unrelated files changed, no backup or temp files are inside the repo, no accidental lockfile changes exist, and no unrelated formatting churn exists.

---

## 12. Documentation

Update docs when source changes affect public APIs, CLI behavior, configuration, data models, schema fields, column names, pipeline stages, generated outputs, error messages, setup steps, or user-visible behavior. Do not update docs for purely internal changes unless the docs would otherwise become misleading. When updating docs: read the relevant section first, keep changes concise, match existing style, verify docs match the actual code, and do not rewrite unrelated sections. If no docs update is needed, say why in the final response.

---

## 13. Changelog

Update changelogs for user-visible or behavior-relevant changes: features, bug fixes, behavior changes, public API changes, CLI changes, config changes, schema changes, data model changes, pipeline changes, dependency changes that affect users, breaking changes, or removed behavior. Do not update changelog for typo-only, formatting-only, comment-only changes, or internal refactors with no behavior change, unless the user asks.

Use the project's existing changelog location and format (`CHANGELOG.md`, `packages/*/CHANGELOG.md`, `crates/*/CHANGELOG.md`). If the project uses an `[Unreleased]` section, add entries there. Use standard sections: Breaking Changes, Added, Changed, Fixed, Removed. Read the full target changelog section first, append to existing subsections, do not duplicate subsections, do not edit released version sections unless asked, keep entries concise (one bullet per change), and prefix with bold scope when the project uses scoped entries.

Example: `- **api**: Validate missing request fields before processing.`

---

## 14. Dependency Security

Treat dependency and lockfile changes as code changes. Do not add dependencies unless necessary. Prefer existing dependencies. Review lockfile diffs. Do not run install scripts unless the user asks. Do not bypass security gates silently. Do not downgrade dependencies just to hide type or lint errors. Explain dependency changes in the final response.

```bash
# Python
uv add <package>
uv sync

# TypeScript
npm install --ignore-scripts
npm ci --ignore-scripts
npm install --package-lock-only --ignore-scripts

# Rust
cargo add <crate>
cargo update -p <crate>
```

Follow project-specific package manager conventions when they differ.

---

## 15. Final Checklist

Before final response, verify: relevant docs were read, affected source files were read, existing patterns were checked, temporary backups are stored in `/tmp`, no WIP checkpoint commits were created, no in-repo `.bak` files exist, changes are surgical and task-related, no unrelated formatting churn or dependency/lockfile changes exist, generated files were not edited directly unless allowed, required checks were run, full tests were not run unless requested, targeted tests were run only when allowed or requested, docs and changelog were updated when needed, and `git diff`/`git status` were reviewed.

---

## 16. Final Response Template

```markdown
Changed:
- <summary>

Files:
- <file>: <what changed>

Verification:
- <command>: passed / failed / not run, with reason

Tests:
- <test command>
- Not run because full test suites require explicit request

Docs:
- Updated <doc>
- Not updated because <reason>

Changelog:
- Updated <changelog>
- Not updated because <reason>

Notes:
- <assumptions or risks>
```
