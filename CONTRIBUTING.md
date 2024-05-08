# Contributing

Thanks for being willing to contribute!

## Project setup

1. Fork and clone the repo.
2. Install a [Python 3.9](https://www.python.org/downloads/release/python-390/) (or later), [Docker](https://docs.docker.com/engine/install/), Bash 3.2 (or later) and [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html).
3. Execute `pip install pre-commit` and `sudo apt install make`.

## Building and Testing

### Unit (step) tests

1. To build an **EXECUTABLE** docker images, python virtual environment and installs requirements execute:
```bash
make build
```

2. Run unit tests:
```
make test
```

### Integration (workflow) tests

1. To build a NOT EXECUTABLE docker images, nextflow, execute:
```bash
make
```

2. Run integration tests
```bash
make test_wf
```

## Committing and Pushing changes

Please make sure to run the tests before you commit your changes (if you didn't configure `pre-commit`). You can run for it `make test` (for steps) and `make test_wf` (for workflows).  

Also, check that your code meet [PEP8](https://peps.python.org/pep-0008/) requirements (by [ruff](https://github.com/astral-sh/ruff)), dynamic and static typing (by [mypy](https://github.com/python/mypy)). You can automate it using:
```bash
make check # runs ruff linter
make format # runs ruff formatter
make mypy # runs mypy type checker
```

### Pre-commit hooks

Additionally, you can activate [pre-commit](https://pre-commit.com/) hooks. Execute:

```bash
pre-commit install # set up the git hook scripts
```

## Development

* Create an [issue](https://github.com/BostonGene/pyigmap/issues), a branch from `main` (from the task), and then create a pull request from this branch.
* After each minor/major change in pyigmap, make sure to add some notes to [CHANGELOG](CHANGELOG.md).

## Help needed

Please check out [the open issues](https://github.com/BostonGene/pyigmap/issues).

Also, please watch the repo and respond to questions/bug reports/feature requests! Thanks!