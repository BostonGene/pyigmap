# Review PR

Review a pull request with parallel subagents for comprehensive analysis.

## Context

```bash
# Get PR info if number provided
if [ -n "$ARGUMENTS" ]; then
  gh pr view "$ARGUMENTS" --json title,body,files,additions,deletions 2>/dev/null || echo "PR not found"
fi

# Current branch diff if no PR number
git diff main...HEAD --stat 2>/dev/null || git diff origin/main...HEAD --stat
```

## Instructions

Review the PR using multiple perspectives:

1. **Style Check**: Verify code follows CONTRIBUTING.md guidelines
   - Naming conventions (flow_vibyy_*, task_*, UPPER_SNAKE_CASE)
   - Import style (absolute imports, isort grouping)
   - Line length (120 chars), quote style (single quotes)
   - Type hints on all functions

2. **Bug Hunt**: Look for common issues
   - Unhandled exceptions
   - Resource leaks (unclosed files/connections)
   - Missing null checks
   - Logic errors in conditionals

3. **Security Review**: Check for vulnerabilities
   - Hardcoded credentials or secrets
   - SQL/command injection risks
   - Unsafe deserialization
   - Missing input validation

4. **Architecture**: Evaluate design
   - Does it follow existing patterns?
   - Is it over-engineered?
   - Are there simpler alternatives?

5. **Tests**: Verify test coverage
   - Are new functions tested?
   - Are edge cases covered?
   - Do tests follow naming conventions?

## Arguments

- `$ARGUMENTS` - Optional: PR number to review
  - Example: `42`

If no PR number, review the current branch's diff against main.

## Output Format

Provide findings grouped by severity:
- **Critical**: Must fix before merge
- **Major**: Should fix, can discuss
- **Minor**: Nice to have improvements
- **Praise**: Things done well
