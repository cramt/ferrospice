# Project Management Pass

You have chosen to do a **project management pass**. Your goal is to improve the PRD and project organization — not to write any production code.

## What You Can Do

- **Reprioritize** — reorder story priorities based on dependencies, complexity, or what unblocks the most work
- **Split stories** — break a large story into smaller, more actionable pieces (assign new IDs like US-023a, US-023b or use the next available US-0XX number)
- **Merge stories** — combine stories that are too small or tightly coupled to work on independently
- **Add stories** — identify missing work (gaps in test coverage, missing device models, infrastructure needs) and create new user stories
- **Remove stories** — mark stories as obsolete if the approach has changed, with a note explaining why
- **Update acceptance criteria** — sharpen criteria based on what you've learned from the codebase
- **Update notes** — add implementation hints, reference file paths, dependency information
- **Update progress.txt patterns** — consolidate or correct the Codebase Patterns section

## Steps

1. Read the PRD at `prd.json`
2. Read the progress log at `progress.txt` (especially Codebase Patterns)
3. Read `../../CLAUDE.md` for project conventions
4. Assess the current state:
   - Which stories are done vs remaining?
   - Are priorities still sensible given what's been learned?
   - Are any stories too large to complete in a single iteration?
   - Are there missing stories or gaps?
   - Are acceptance criteria still accurate?
5. Make your changes to `prd.json`
6. Append a summary to `progress.txt`

## Progress Report Format

APPEND to progress.txt:
```
## [Date/Time] - Project Management Pass
- Changes made:
  - [list each change: reprioritized X, split Y, added Z, etc.]
- Rationale:
  - [why these changes improve the project plan]
---
```

## Constraints

- **Do NOT modify any source code** (no Rust files, no Cargo.toml changes)
- **Do NOT mark stories as `passes: true`** unless they were genuinely already complete and mislabeled
- **Do NOT commit code changes** — only commit prd.json and progress.txt changes
- **Preserve completed work** — never set a `passes: true` story back to `false` unless it's genuinely broken
- Commit with message: `chore: project management pass — [brief summary]`
