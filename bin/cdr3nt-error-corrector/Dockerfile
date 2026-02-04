FROM python:3.12-slim-bullseye AS image

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set UV environment variables
ENV UV_SYSTEM_PYTHON=1
ENV UV_COMPILE_BYTECODE=1
ENV UV_PYTHON_DOWNLOADS=never

# Copy project files and install
WORKDIR /app
COPY pyproject.toml uv.lock ./
COPY src/ ./src/

# Keep build_ref.py at root for reference building
COPY build_ref.py ./

# Install dependencies and the package itself
RUN uv pip install --system -e .

FROM image AS tool

ENTRYPOINT ["python3.12", "-m", "cdr3nt_error_corrector.run"]
