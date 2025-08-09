FROM ubuntu:24.04

RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    ca-certificates \
    g++ \
    gcc \
    pigz \
    autoconf \
    automake \
    libtool \
    nasm \
    yasm \
    libbz2-dev \
    liblzma-dev \
    libdeflate0 \
    libpython3-dev \
    libisal-dev \
    libgcc-s1 \
    cpio \
    libc6 \
    make \
    cmake \
    git \
    ncurses-dev \
    python3 \
    python3-pip \
    wget \
    xxd \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    uuid-dev \
    libexpat1-dev \
    libxml2-dev \
    libxslt1-dev \
    libsqlite3-dev \
    libmagic-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    libgif-dev \
    libxpm-dev \
    libfreetype-dev \
    libboost-all-dev \
    libprotobuf-dev \
    protobuf-compiler \
    libsasl2-dev \
    libhdf5-dev \
    liblapack-dev \
    libblas-dev \
    libgsl-dev \
    liburing-dev \
    libuv1-dev \
    libyaml-cpp-dev \
    libreadline-dev \
    libncurses-dev \
    uuid-runtime \
    gawk \
    perl \
    curl \
    pkg-config && \
    apt-get autoremove --purge --yes && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV SEQTK_VERSION=1.4
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz -O seqtk.tar.gz && \
    tar -xvzf seqtk.tar.gz --one-top-level=seqtk --strip-component 1 && \
    cd seqtk && make && mv seqtk /usr/local/bin && \
    cd .. && rm -r seqtk seqtk.tar.gz

ENV SEQKIT_VERSION=2.8.0
RUN if [ "$(uname -m)" = "x86_64" ]; then \
      BIN_ARCH=amd64; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
      BIN_ARCH=arm64; \
    fi && \
    wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_${BIN_ARCH}.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=/usr/local/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

ENV IGBLAST_VERSION=1.22.0
RUN if [ "$(uname -m)" = "x86_64" ]; then \
      wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz -O ncbi-igblast.tar.gz && \
      tar -xvzf ncbi-igblast.tar.gz --one-top-level=/usr/local/bin/ncbi-igblast --strip-component 1 && \
      rm ncbi-igblast.tar.gz ; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
      wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-src.tar.gz -O ncbi-igblast.tar.gz && \
      tar -xvzf ncbi-igblast.tar.gz --one-top-level=ncbi-igblast --strip-component 1 && rm ncbi-igblast.tar.gz && \
      cd ncbi-igblast/c++ && ./configure && make -j$(nproc) && make install && \
      mkdir -p /usr/local/bin/ncbi-igblast/bin/ && cp /ncbi-igblast/c++/ReleaseMT/bin/igblastn /usr/local/bin/ncbi-igblast/bin/ && \
      chmod 755 /usr/local/bin/ncbi-igblast/bin/igblastn && \
      cd / && rm -rf ncbi-igblast* ; \
    fi

ENV IGBLAST_DIR=/usr/local/bin/ncbi-igblast
ENV PATH=${IGBLAST_DIR}:$PATH

COPY run.py /usr/local/src/

ENTRYPOINT ["python3", "/usr/local/src/run.py"]