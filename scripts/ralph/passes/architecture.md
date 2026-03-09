# Architecture Pass

You have chosen to do an **architecture pass**. Your goal is to improve the codebase structure, reduce tech debt, and leave the repo better than you found it — without advancing any user stories.

## What You Can Do

- **Refactor** — extract shared logic, reduce duplication, improve module boundaries
- **Split crates** — if a crate is getting too large, split it along natural boundaries
- **Improve types** — make impossible states irrepresentable, add newtypes, tighten enums
- **Improve error handling** — better error types, more specific error variants, clearer messages
- **Fix warnings** — resolve clippy warnings, dead code, unused imports
- **Improve test infrastructure** — shared test helpers, fixtures, assertion utilities
- **Documentation** — add CLAUDE.md files to under-documented modules, update stale docs
- **Dependency cleanup** — remove unused deps, consolidate versions

## Steps

1. Read the progress log at `progress.txt` (Codebase Patterns section)
2. Read `../../CLAUDE.md` for project conventions
3. Survey the codebase — look for:
   - Large files that should be split
   - Duplicated patterns across device models
   - Overly complex functions
   - Missing abstractions that would simplify upcoming stories
   - Stale or misleading comments
4. Pick ONE focused refactoring task (don't try to fix everything)
5. Implement the refactoring
6. Run quality checks: `nix develop --command cargo clippy --workspace -- -D warnings && nix develop --command cargo test --workspace`
7. Verify ALL existing tests still pass (zero regressions)
8. Commit changes
9. Update progress.txt and Codebase Patterns if the refactoring changes conventions

## Progress Report Format

APPEND to progress.txt:
```
## [Date/Time] - Architecture Pass
- What was refactored:
  - [describe the change]
- Why:
  - [what problem it solves, what it unblocks]
- Impact on future work:
  - [how this affects upcoming stories]
- **Codebase pattern updates:**
  - [any patterns that changed as a result]
---
```

## Constraints

- **Do NOT advance user stories** — don't mark anything as `passes: true`
- **Do NOT change prd.json** (that's a project management pass concern)
- **Zero test regressions** — every test that passed before must still pass
- **One focused change** — don't try to refactor the whole codebase at once
- **No feature work** — refactoring only, not new functionality
- Commit with message: `refactor: [brief description of what was improved]`
