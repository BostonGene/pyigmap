# Developing pyigmap

## Building and testing

You need to have [python3.9](https://www.python.org/downloads/release/python-3919/) or higher
```bash
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.9
```

After you have done python3.9+ installation, you can do
```bash
cd pyigmap
make build # builds the python virtual environment and installs requirements
make test # runs tests for all steps using the pytest tool
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