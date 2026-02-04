FROM ubuntu:24.04 AS image

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    ca-certificates \
    curl \
    dash \
    g++ \
    gcc \
    pigz \
    automake \
    libtool \
    nasm \
    yasm \
    libbz2-dev \
    liblzma-dev \
    libdeflate0 \
    libpython3-dev \
    libisal-dev \
    make \
    cmake \
    git \
    ncurses-dev \
    wget \
    xxd \
    zlib1g-dev \
    python3.12 \
    python3.12-dev \
    python3.12-venv && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install fastp binary
ENV FASTP_VERSION=1.0.1
ENV LIBDEFLATE_VERSION=1.23
ENV ISAL_VERSION=2.31.1

RUN if [ "$(uname -m)" = "x86_64" ]; then \
      wget http://opengene.org/fastp/fastp.${FASTP_VERSION} -O /usr/local/bin/fastp && \
      chmod a+x /usr/local/bin/fastp; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
      git clone -b "v$ISAL_VERSION" https://github.com/intel/isa-l.git && \
      cd isa-l && ./autogen.sh && ./configure --prefix=/usr --libdir=/usr/lib64 && \
      make -j$(nproc) && make install && cd - && \
      git clone -b "v$LIBDEFLATE_VERSION" https://github.com/ebiggers/libdeflate.git && \
      cd libdeflate && cmake -B build && cmake --build build && cmake --install build && cd - && \
      git clone -b "v$FASTP_VERSION" https://github.com/OpenGene/fastp.git && \
      cd fastp && make -j$(nproc) && make install && cd -; \
    fi

# Install Rust
RUN curl -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:$PATH"

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Create venv for isolation
RUN python3.12 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
ENV VIRTUAL_ENV="/opt/venv"

# Copy project files
WORKDIR /app
COPY pyproject.toml uv.lock ./
COPY mock_merger/ ./mock_merger/
COPY src/ ./src/

# Build with maturin using uv
# Install build dependencies
RUN uv pip install --system maturin patchelf tomli

# Build the Rust extension
RUN maturin build --release --out dist

# Install the built wheel (this installs the fastp package)
RUN uv pip install --system dist/*.whl && rm -rf dist

ENV RUST_LOG=info

FROM image AS tool

ENTRYPOINT ["python3.12", "-m", "fastp.run"]
