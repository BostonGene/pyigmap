FROM python:3.12-slim-bullseye AS image

ENV SOFT_DIR=/usr/local

# Install system dependencies
RUN apt-get update && apt-get install -y \
    cmake \
    pigz \
    jq \
    wget \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set UV environment variables
ENV UV_SYSTEM_PYTHON=1
ENV UV_COMPILE_BYTECODE=1
ENV UV_PYTHON_DOWNLOADS=never

# Install calib
ENV CALIB_VERSION=0.3.7
RUN wget -q https://github.com/vpc-ccg/calib/archive/refs/tags/v${CALIB_VERSION}.tar.gz -O CALIB.tar.gz && \
    tar -xzvf CALIB.tar.gz --one-top-level=CALIB --strip-component 1 && \
    cd CALIB && \
    make && mv calib ${SOFT_DIR}/bin && \
    make -C consensus/ && mv consensus/calib_cons ${SOFT_DIR}/bin && cd .. && \
    rm CALIB.tar.gz && rm -r CALIB

# Copy project files and install (no dependencies, but sets up package structure)
WORKDIR /app
COPY pyproject.toml ./
COPY src/ ./src/

# Install the package (no dependencies to install currently)
RUN uv pip install --system --no-deps -e .

FROM image AS tool

ENTRYPOINT ["python3.12", "-m", "calib_dedup.run"]
