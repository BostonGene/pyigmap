# Format Code

Format and lint code according to project standards.

## Context

```bash
# Check for modified Python files
git diff --name-only --diff-filter=ACMR | grep -E '\.py$' || echo "No Python files modified"

# Check ruff version
uv run ruff --version
```

## Instructions

1. Run the formatting tools in order:
   - `just format` - Auto-format with ruff
   - `just lint` - Check for remaining issues
   - `just type` - Run type checker

2. For any issues that can't be auto-fixed:
   - List the file and line number
   - Explain the issue
   - Suggest a fix

3. After formatting, show:
   - Files that were modified
   - Any remaining issues that need manual attention

## Arguments

- `$ARGUMENTS` - Optional: specific files or directories to format
  - Example: `flows/vibyy_asr/`
  - Example: `common/storage/client.py`

If no arguments, format all modified files.
