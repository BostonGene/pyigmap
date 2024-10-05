# Creating a Release

1. Run the [release make command](Makefile): `make v[X.Y.Z]`
2. Push the changes:
```bash
git push
```
3. Push a tag:
```bash
git tag v[X.Y.Z] && git push origin v[X.Y.Z]
```
4. Create a new release in GitHub
```bash
gh release create $(TAG) --notes-from-tag
```
