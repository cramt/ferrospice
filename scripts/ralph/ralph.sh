#!/usr/bin/env bash
# Ralph Wiggum - Long-running AI agent loop
# Usage: ./ralph.sh [-i|--interactive] [--task US-XXX] [max_iterations]

set -e

# Parse arguments
MAX_ITERATIONS=10
TASK=""

INTERACTIVE=false

while [[ $# -gt 0 ]]; do
  case $1 in
    -i|--interactive)
      INTERACTIVE=true
      shift
      ;;
    --task)
      TASK="$2"
      shift 2
      ;;
    --task=*)
      TASK="${1#*=}"
      shift
      ;;
    *)
      # Assume it's max_iterations if it's a number
      if [[ "$1" =~ ^[0-9]+$ ]]; then
        MAX_ITERATIONS="$1"
      fi
      shift
      ;;
  esac
done

# When a specific task is given, force a single iteration
if [ -n "$TASK" ]; then
  MAX_ITERATIONS=1
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRD_FILE="$SCRIPT_DIR/prd.json"
PROGRESS_FILE="$SCRIPT_DIR/progress.txt"
ARCHIVE_DIR="$SCRIPT_DIR/archive"
LAST_BRANCH_FILE="$SCRIPT_DIR/.last-branch"
CLAUDE_MD="$SCRIPT_DIR/CLAUDE.md"

# Archive previous run if branch changed
if [ -f "$PRD_FILE" ] && [ -f "$LAST_BRANCH_FILE" ]; then
  CURRENT_BRANCH=$(jq -r '.branchName // empty' "$PRD_FILE" 2>/dev/null || echo "")
  LAST_BRANCH=$(cat "$LAST_BRANCH_FILE" 2>/dev/null || echo "")

  if [ -n "$CURRENT_BRANCH" ] && [ -n "$LAST_BRANCH" ] && [ "$CURRENT_BRANCH" != "$LAST_BRANCH" ]; then
    # Archive the previous run
    DATE=$(date +%Y-%m-%d)
    # Strip "ralph/" prefix from branch name for folder
    FOLDER_NAME=$(echo "$LAST_BRANCH" | sed 's|^ralph/||')
    ARCHIVE_FOLDER="$ARCHIVE_DIR/$DATE-$FOLDER_NAME"

    echo "Archiving previous run: $LAST_BRANCH"
    mkdir -p "$ARCHIVE_FOLDER"
    [ -f "$PRD_FILE" ] && cp "$PRD_FILE" "$ARCHIVE_FOLDER/"
    [ -f "$PROGRESS_FILE" ] && cp "$PROGRESS_FILE" "$ARCHIVE_FOLDER/"
    echo "   Archived to: $ARCHIVE_FOLDER"

    # Reset progress file for new run
    echo "# Ralph Progress Log" > "$PROGRESS_FILE"
    echo "Started: $(date)" >> "$PROGRESS_FILE"
    echo "---" >> "$PROGRESS_FILE"
  fi
fi

# Track current branch
if [ -f "$PRD_FILE" ]; then
  CURRENT_BRANCH=$(jq -r '.branchName // empty' "$PRD_FILE" 2>/dev/null || echo "")
  if [ -n "$CURRENT_BRANCH" ]; then
    echo "$CURRENT_BRANCH" > "$LAST_BRANCH_FILE"
  fi
fi

# Initialize progress file if it doesn't exist
if [ ! -f "$PROGRESS_FILE" ]; then
  echo "# Ralph Progress Log" > "$PROGRESS_FILE"
  echo "Started: $(date)" >> "$PROGRESS_FILE"
  echo "---" >> "$PROGRESS_FILE"
fi

# Build pass history from prd.json passLog
build_pass_history() {
  if [ ! -f "$PRD_FILE" ]; then
    echo "No pass history yet (first run)."
    return
  fi

  local LOG_LENGTH
  LOG_LENGTH=$(jq '.passLog // [] | length' "$PRD_FILE" 2>/dev/null || echo "0")

  if [ "$LOG_LENGTH" -eq 0 ]; then
    echo "No pass history yet (first run)."
    return
  fi

  # Show last 10 entries
  local OFFSET=0
  if [ "$LOG_LENGTH" -gt 10 ]; then
    OFFSET=$((LOG_LENGTH - 10))
  fi

  echo "Last passes (most recent last):"
  echo ""
  jq -r --argjson offset "$OFFSET" '
    .passLog // [] | .[$offset:] | to_entries[] |
    "  \(.value.iteration). [\(.value.passType)] \(.value.summary) (\(.value.date))\(if .value.storyId then " — " + .value.storyId else "" end)"
  ' "$PRD_FILE" 2>/dev/null || echo "  (could not read pass log)"

  # Summary counts
  echo ""
  local IMPL_COUNT PM_COUNT ARCH_COUNT
  IMPL_COUNT=$(jq '[.passLog // [] | .[] | select(.passType == "implementation")] | length' "$PRD_FILE" 2>/dev/null || echo "0")
  PM_COUNT=$(jq '[.passLog // [] | .[] | select(.passType == "project-management")] | length' "$PRD_FILE" 2>/dev/null || echo "0")
  ARCH_COUNT=$(jq '[.passLog // [] | .[] | select(.passType == "architecture")] | length' "$PRD_FILE" 2>/dev/null || echo "0")
  echo "Totals: $IMPL_COUNT implementation, $PM_COUNT project-management, $ARCH_COUNT architecture"

  # Streak detection — nudge if too many implementation passes in a row
  local RECENT_TYPES
  RECENT_TYPES=$(jq -r '[.passLog // [] | .[-4:] | .[].passType] | join(",")' "$PRD_FILE" 2>/dev/null || echo "")
  if [[ "$RECENT_TYPES" == "implementation,implementation,implementation,implementation" ]]; then
    echo ""
    echo "⚠ The last 4 passes were all implementation. Strongly consider a project-management or architecture pass."
  fi
}

if [ -n "$TASK" ]; then
  echo "Starting Ralph - Task: $TASK - Max iterations: $MAX_ITERATIONS"
else
  echo "Starting Ralph - Max iterations: $MAX_ITERATIONS"
fi

for i in $(seq 1 $MAX_ITERATIONS); do
  echo ""
  echo "==============================================================="
  echo "  Ralph Iteration $i of $MAX_ITERATIONS"
  echo "==============================================================="

  # Build the prompt with pass history injected
  PASS_HISTORY=$(build_pass_history)

  # Use awk for robust multi-line replacement
  PROMPT=$(awk -v history="$PASS_HISTORY" '{
    if ($0 ~ /<!-- PASS_HISTORY_PLACEHOLDER -->/) {
      print history
    } else {
      print
    }
  }' "$CLAUDE_MD")

  # Inject task directive if specified
  if [ -n "$TASK" ]; then
    PROMPT="$PROMPT

## TASK OVERRIDE

You MUST work on story **$TASK** specifically. Do not pick a different story. This is an implementation pass targeting $TASK."
  fi

  # Run claude code with the ralph prompt
  if [ "$INTERACTIVE" = true ]; then
    # Interactive mode: spawn a normal claude code instance you can watch/interact with
    echo "$PROMPT" | claude --dangerously-skip-permissions --verbose
    OUTPUT=""
  else
    OUTPUT=$(echo "$PROMPT" | claude --dangerously-skip-permissions --print --verbose 2>&1 \
      | tee /dev/stderr) || true
  fi

  # Check for completion signal
  if echo "$OUTPUT" | grep -q "<promise>COMPLETE</promise>"; then
    echo ""
    echo "Ralph completed all tasks!"
    echo "Completed at iteration $i of $MAX_ITERATIONS"
    exit 0
  fi

  echo "Iteration $i complete. Continuing..."
  sleep 2
done

echo ""
echo "Ralph reached max iterations ($MAX_ITERATIONS) without completing all tasks."
echo "Check $PROGRESS_FILE for status."
exit 1
