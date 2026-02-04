FROM python:3.12-slim-bullseye AS image

ENV SOFT_DIR=/usr/local

# Install system dependencies
RUN apt-get update && apt-get install -y \
    pigz \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set UV environment variables
ENV UV_SYSTEM_PYTHON=1
ENV UV_COMPILE_BYTECODE=1
ENV UV_PYTHON_DOWNLOADS=never

# Install seqkit
ENV SEQKIT_VERSION=2.8.1
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=${SOFT_DIR}/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

# Copy project files and install
WORKDIR /app
COPY pyproject.toml uv.lock ./
COPY src/ ./src/

# Install dependencies and the package itself
RUN uv pip install --system -e .

FROM image AS tool

ENTRYPOINT ["python3.12", "-m", "pyumi.run"]
