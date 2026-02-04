FROM python:3.12-slim-bookworm AS image

# Install system dependencies
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    dash \
    pigz \
    wget \
    build-essential \
    cmake \
    git \
    gettext \
    pkg-config \
    libssl-dev \
    ca-certificates \
    axel && \
    rm -rf /var/lib/apt/lists/*

# Install SRA toolkit
ENV SRA_VERSION=3.2.1
ENV VDB_VERSION=3.2.1
RUN mkdir /usr/local/bin/sratoolkit && \
    if [ "$(uname -m)" = "x86_64" ]; then \
        wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-ubuntu64.tar.gz -O sratoolkit.tar.gz && \
        tar -xzf sratoolkit.tar.gz --strip-components=1 -C /usr/local/bin/sratoolkit && \
        rm sratoolkit.tar.gz ; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
        git clone -b $VDB_VERSION https://github.com/ncbi/ncbi-vdb.git && \
        git clone -b $SRA_VERSION https://github.com/ncbi/sra-tools.git && \
        mkdir build && \
        cd build && \
        cmake -S "$(cd ../ncbi-vdb; pwd)" -B ncbi-vdb && \
        cmake --build ncbi-vdb && \
        cmake -D VDB_LIBDIR="$PWD/ncbi-vdb/lib" -D CMAKE_INSTALL_PREFIX="$PWD/sratoolkit" -S "$(cd ../sra-tools; pwd)" -B sra-tools && \
        cmake --build sra-tools --target install && \
        cp -r $PWD/sratoolkit/* /usr/local/bin/sratoolkit/ && \
        cd - && rm -rf build ncbi-vdb sra-tools; \
    fi

ENV PATH=${PATH}:/usr/local/bin/sratoolkit/bin

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

# Install dependencies and the package itself
RUN uv pip install --system -e .

FROM image AS tool

ENTRYPOINT ["python3.12", "-m", "fq_downloader.run"]
