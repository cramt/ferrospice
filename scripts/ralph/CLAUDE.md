# Ralph Agent Instructions — Ferrospice

You are an autonomous coding agent working on **ferrospice**, a Rust rewrite of ngspice (circuit simulator). This is a **library-first** project.

## Project Context

- **Reference C source:** `ngspice-upstream/` (read-only, do not modify)
- **Existing netlist parser:** `ferrospice-netlist/` crate (parses SPICE files)
- **Project root CLAUDE.md:** Read `../../CLAUDE.md` for conventions, goals, and dev methodology
- **Test suite:** `ngspice-upstream/tests/` contains 113 `.cir`/`.out` test pairs — the project is done when all pass

## Quality Checks

**Always run commands through `nix develop --command ...`** so that flake.nix stays honest.

```bash
nix develop --command cargo clippy --workspace -- -D warnings
nix develop --command cargo test --workspace
nix develop --command cargo fmt --check
```

If a command fails because a tool is missing, add it to the `devShell` in `flake.nix`.

## Key Rules

- **Correctness over performance.** Get it right first. Optimize later.
- **Make impossible states irrepresentable.** Use enums, newtypes, typestate.
- **Use facet, not serde.** Use unsynn, not syn.
- **Split into subcrates freely.** Don't let things get monolithic.
- **Test-driven.** Write tests first when possible. Port ngspice tests before implementing.
- **Never include Co-Authored-By lines in git commits.**
- **Library-first.** All simulation logic lives in library crates. CLI is a thin wrapper.

## Your Task

1. Read the PRD at `prd.json` (in the same directory as this file)
2. Read the progress log at `progress.txt` (check Codebase Patterns section first)
3. Read `../../CLAUDE.md` for project conventions
4. Check you're on the correct branch from PRD `branchName`. If not, check it out or create from main.
5. Pick the **highest priority** user story where `passes: false`
6. Implement that single user story
7. Run quality checks: `nix develop --command cargo clippy --workspace -- -D warnings && nix develop --command cargo test --workspace`
8. Update CLAUDE.md files if you discover reusable patterns (see below)
9. If checks pass, commit ALL changes with message: `feat: [Story ID] - [Story Title]`
10. Update the PRD to set `passes: true` for the completed story
11. Append your progress to `progress.txt`

## Progress Report Format

APPEND to progress.txt (never replace, always append):
```
## [Date/Time] - [Story ID]
- What was implemented
- Files changed
- **Learnings for future iterations:**
  - Patterns discovered (e.g., "this codebase uses X for Y")
  - Gotchas encountered (e.g., "don't forget to update Z when changing W")
  - Useful context (e.g., "the evaluation panel is in component X")
---
```

The learnings section is critical - it helps future iterations avoid repeating mistakes and understand the codebase better.

## Consolidate Patterns

If you discover a **reusable pattern** that future iterations should know, add it to the `## Codebase Patterns` section at the TOP of progress.txt (create it if it doesn't exist). This section should consolidate the most important learnings:

```
## Codebase Patterns
- Example: Use `faer` for sparse matrix operations
- Example: Device trait has dc_stamp(), ac_stamp(), tran_stamp() methods
- Example: Node "0" is ground and excluded from the MNA matrix
```

Only add patterns that are **general and reusable**, not story-specific details.

## Update CLAUDE.md Files

Before committing, check if any edited files have learnings worth preserving in nearby CLAUDE.md files:

1. **Identify directories with edited files** - Look at which directories you modified
2. **Check for existing CLAUDE.md** - Look for CLAUDE.md in those directories or parent directories
3. **Add valuable learnings** - If you discovered something future developers/agents should know:
   - API patterns or conventions specific to that module
   - Gotchas or non-obvious requirements
   - Dependencies between files
   - Testing approaches for that area

**Do NOT add:**
- Story-specific implementation details
- Temporary debugging notes
- Information already in progress.txt

## Stop Condition

After completing a user story, check if ALL stories have `passes: true`.

If ALL stories are complete and passing, reply with:
<promise>COMPLETE</promise>

If there are still stories with `passes: false`, end your response normally (another iteration will pick up the next story).

## Important

- Work on ONE story per iteration
- Commit frequently
- Keep CI green
- Read the Codebase Patterns section in progress.txt before starting
- Reference ngspice C source in `ngspice-upstream/` when implementing algorithms
