#!/usr/bin/env bash
# Formats Claude Code stream-json output into human-readable log lines
# Reads JSON lines from stdin, outputs formatted text to stdout
# Also captures the final result text for the COMPLETE check

FINAL_TEXT=""

while IFS= read -r line; do
  # Skip empty lines and non-JSON
  [[ -z "$line" ]] && continue
  [[ "$line" != "{"* ]] && echo "$line" && continue

  TYPE=$(echo "$line" | jq -r '.type // empty' 2>/dev/null)
  [[ -z "$TYPE" ]] && continue

  case "$TYPE" in
    assistant)
      # Thinking or text content
      MESSAGE_TYPE=$(echo "$line" | jq -r '.message.type // empty' 2>/dev/null)
      if [[ "$MESSAGE_TYPE" == "assistant" ]]; then
        # Extract content blocks
        CONTENT_LEN=$(echo "$line" | jq '.message.content | length' 2>/dev/null)
        for ((j=0; j<${CONTENT_LEN:-0}; j++)); do
          BLOCK_TYPE=$(echo "$line" | jq -r ".message.content[$j].type // empty" 2>/dev/null)
          case "$BLOCK_TYPE" in
            thinking)
              THINKING=$(echo "$line" | jq -r ".message.content[$j].thinking // empty" 2>/dev/null)
              if [[ -n "$THINKING" ]]; then
                echo ""
                echo "--- thinking ---"
                echo "$THINKING"
                echo "--- /thinking ---"
              fi
              ;;
            text)
              TEXT=$(echo "$line" | jq -r ".message.content[$j].text // empty" 2>/dev/null)
              if [[ -n "$TEXT" ]]; then
                echo "$TEXT"
                FINAL_TEXT+="$TEXT"
              fi
              ;;
            tool_use)
              TOOL_NAME=$(echo "$line" | jq -r ".message.content[$j].name // empty" 2>/dev/null)
              TOOL_INPUT=$(echo "$line" | jq -c ".message.content[$j].input // {}" 2>/dev/null)
              echo ""
              echo ">>> tool: $TOOL_NAME"
              # Show a summary of the input, truncated
              INPUT_PREVIEW=$(echo "$TOOL_INPUT" | head -c 500)
              echo "    input: $INPUT_PREVIEW"
              ;;
          esac
        done
      fi
      ;;
    content_block_start)
      BLOCK_TYPE=$(echo "$line" | jq -r '.content_block.type // empty' 2>/dev/null)
      case "$BLOCK_TYPE" in
        thinking)
          echo ""
          echo "--- thinking ---"
          ;;
        tool_use)
          TOOL_NAME=$(echo "$line" | jq -r '.content_block.name // empty' 2>/dev/null)
          echo ""
          echo ">>> tool: $TOOL_NAME"
          ;;
      esac
      ;;
    content_block_delta)
      DELTA_TYPE=$(echo "$line" | jq -r '.delta.type // empty' 2>/dev/null)
      case "$DELTA_TYPE" in
        thinking_delta)
          THINKING=$(echo "$line" | jq -r '.delta.thinking // empty' 2>/dev/null)
          [[ -n "$THINKING" ]] && printf '%s' "$THINKING"
          ;;
        text_delta)
          TEXT=$(echo "$line" | jq -r '.delta.text // empty' 2>/dev/null)
          if [[ -n "$TEXT" ]]; then
            printf '%s' "$TEXT"
            FINAL_TEXT+="$TEXT"
          fi
          ;;
        input_json_delta)
          PARTIAL=$(echo "$line" | jq -r '.delta.partial_json // empty' 2>/dev/null)
          [[ -n "$PARTIAL" ]] && printf '%s' "$PARTIAL"
          ;;
      esac
      ;;
    content_block_stop)
      echo ""
      ;;
    tool_use)
      TOOL_NAME=$(echo "$line" | jq -r '.tool_name // .name // empty' 2>/dev/null)
      TOOL_INPUT=$(echo "$line" | jq -c '.input // {}' 2>/dev/null)
      echo ""
      echo ">>> tool: $TOOL_NAME"
      INPUT_PREVIEW=$(echo "$TOOL_INPUT" | head -c 500)
      echo "    input: $INPUT_PREVIEW"
      ;;
    tool_result)
      TOOL_NAME=$(echo "$line" | jq -r '.tool_name // .name // empty' 2>/dev/null)
      OUTPUT_PREVIEW=$(echo "$line" | jq -r '.output // .content // empty' 2>/dev/null | head -c 300)
      echo "<<< result ($TOOL_NAME): $OUTPUT_PREVIEW"
      ;;
    result)
      # Final result message
      TEXT=$(echo "$line" | jq -r '.result // .text // empty' 2>/dev/null)
      if [[ -n "$TEXT" ]]; then
        echo ""
        echo "=== FINAL RESULT ==="
        echo "$TEXT"
        FINAL_TEXT+="$TEXT"
      fi
      ;;
  esac
done

# Output the accumulated text as the captured output for the COMPLETE check
echo "$FINAL_TEXT"
