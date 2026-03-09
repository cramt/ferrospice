# Implementation Pass

You have chosen to do an **implementation pass**. Your goal is to pick the next failing user story and implement it.

## Steps

1. Read the PRD at `prd.json` (in the same directory as this file's parent)
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
## [Date/Time] - [Story ID] (implementation)
- What was implemented
- Files changed
- **Learnings for future iterations:**
  - Patterns discovered
  - Gotchas encountered
  - Useful context
---
```

## Consolidate Patterns

If you discover a **reusable pattern** that future iterations should know, add it to the `## Codebase Patterns` section at the TOP of progress.txt (create it if it doesn't exist).

Only add patterns that are **general and reusable**, not story-specific details.

## Update CLAUDE.md Files

Before committing, check if any edited files have learnings worth preserving in nearby CLAUDE.md files:

1. **Identify directories with edited files**
2. **Check for existing CLAUDE.md**
3. **Add valuable learnings** — API patterns, gotchas, dependencies, testing approaches

**Do NOT add:** story-specific details, temporary debugging notes, information already in progress.txt.

## Constraints

- Work on ONE story per iteration
- Commit frequently
- Keep CI green
- Reference ngspice C source in `ngspice-upstream/` when implementing algorithms
