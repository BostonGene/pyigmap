# Developing pyigmap

## Building and testing

You need to have a [Python 3.9](https://www.python.org/downloads/release/python-390/) (or later), [Docker](https://docs.docker.com/engine/install/), Bash 3.2 (or later) and [Java 11 (or later, up to 21)](http://www.oracle.com/technetwork/java/javase/downloads/index.html).

After you have installed them, you can do

```bash
cd pyigmap
make build # builds the python virtual environment and installs requirements
make test # runs tests for all steps using the pytest tool
make test_wf # runs tests for all steps using the pytest and pytest-workflow tool
```

Also, you can run [ruff](https://github.com/astral-sh/ruff) and [mypy](https://github.com/python/mypy) type checker (if you need it)

```bash
make check # runs ruff linter
make format # runs ruff formatter
make mypy # runs mypy type checker
```

Additionally, you can employ [pre-commit](https://pre-commit.com/) hooks. To activate them, execute:

```bash
pre-commit install # set up the git hook scripts
```