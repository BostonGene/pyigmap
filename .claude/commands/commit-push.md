# Commit and Push

Automates the commit and push.

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
   - Use scope when applicable: `feat(api):`, `fix(auth):`
4. Push to the remote branch
