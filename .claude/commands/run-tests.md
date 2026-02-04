# Run Tests

Run the project test suite with proper formatting.

## Context

```bash
# Check if test dependencies are installed
uv run python -c "import pytest" 2>&1 || echo "pytest not installed"
```

## Instructions

1. Run the full quality gate suite:
   - `just lint` - Run ruff linter
   - `just format` - Auto-format code
   - `just type` - Run type checker
   - `just tests-unit` - Run unit tests

2. If any step fails:
   - Show the specific errors clearly
   - Suggest fixes if obvious
   - Do NOT proceed to the next step

3. Report summary at the end:
   - Total tests run
   - Tests passed/failed
   - Any warnings or issues

## Arguments

- `$ARGUMENTS` - Optional: specific test path or pattern
  - Example: `tests/unit/test_storage.py`
  - Example: `tests/unit/test_storage.py::test_upload`

If arguments provided, run only those specific tests after linting.
