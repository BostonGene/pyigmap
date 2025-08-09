FROM ubuntu:24.04

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
    python3 \
    python3-pip \
    python3-venv && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

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

RUN curl -sSf https://sh.rustup.rs | sh -s -- -y \
 && . "/root/.cargo/env" \
 && cargo --version \
 && python3 -m venv /venv

ENV PATH="/root/.cargo/bin:/venv/bin:$PATH"

COPY pyproject.toml mock_merger/Cargo.toml ./
COPY mock_merger/src/ src/
COPY run.py logger.py /usr/local/src/

RUN pip install --upgrade pip && \
    pip install build maturin && \
    maturin build --release -o dist && \
    pip install dist/*.whl && \
    rm -rf dist

ENV RUST_LOG=info
# ENTRYPOINT ["python3", "/usr/local/src/run.py"]