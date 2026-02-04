FROM python:3.12-slim-bullseye AS image

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set UV environment variables
ENV UV_SYSTEM_PYTHON=1
ENV UV_COMPILE_BYTECODE=1
ENV UV_PYTHON_DOWNLOADS=never

# Copy dependency files first
WORKDIR /app
COPY pyproject.toml uv.lock ./
COPY src/ ./src/

# Install dependencies and the package itself in editable mode
RUN uv pip install --system -e .

FROM image AS tool

ENTRYPOINT ["python3.12", "-m", "reporter.run"]
