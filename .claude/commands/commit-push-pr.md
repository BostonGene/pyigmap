# Commit, Push, and Create PR

Automates the commit, push, and PR creation workflow.

## Context

```bash
# Current git status
git status --porcelain

# Current branch
git branch --show-current

# Recent commits for style reference
git log --oneline -5
```

## Instructions

1. Review the git status output above to understand what changed
2. Stage all relevant changes (excluding any sensitive files like .env)
3. Create a conventional commit message following the project style:
   - `feat:` for new features
   - `fix:` for bug fixes
   - `chore:` for maintenance
   - `docs:` for documentation
   - `refactor:` for refactoring
   - Use scope when applicable: `feat(asr):`, `fix(storage):`
4. Push to the remote branch
5. Create a PR with:
   - Clear title matching the commit message
   - Summary section with 1-3 bullet points
   - Test plan section
   - Reference to any related issues using "Fixes #N" or "Relates to #N"

## Arguments

- `$ARGUMENTS` - Optional: PR title override or commit message

If no arguments provided, generate appropriate messages from the changes.
