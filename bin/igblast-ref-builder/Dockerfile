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

ENV IGBLAST_VERSION=1.22.0
ENV IGBLAST_DIR=/usr/local/bin/ncbi-igblast
RUN mkdir -p "${IGBLAST_DIR}" && \
    if [ "$(uname -m)" = "x86_64" ]; then \
      wget -q -O ncbi-igblast.tar.gz "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz" && \
      tar -xvzf ncbi-igblast.tar.gz --one-top-level="${IGBLAST_DIR}" --strip-components=1 && \
      rm -f ncbi-igblast.tar.gz ; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
      wget -q -O ncbi-igblast.tar.gz "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-src.tar.gz" && \
      tar -xvzf ncbi-igblast.tar.gz --one-top-level=ncbi-igblast --strip-components=1 && rm -f ncbi-igblast.tar.gz && \
      cd ncbi-igblast/c++ && ./configure && make -j"$(nproc)" && make install && \
      mkdir -p "${IGBLAST_DIR}/bin" && \
      cp ReleaseMT/bin/makeblastdb "${IGBLAST_DIR}/bin/" && \
      cp ReleaseMT/bin/edit_imgt_file.pl "${IGBLAST_DIR}/bin/" && \
      chmod 755 "${IGBLAST_DIR}/bin/makeblastdb" && \
      cp -r src/app/igblast/internal_data "${IGBLAST_DIR}/" && \
      cd / && rm -rf ncbi-igblast* ; \
    fi && \
    wget -q -O ncbi_human_c_genes.tar "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar" && \
    mkdir -p "${IGBLAST_DIR}/database" && \
    tar -xvf ncbi_human_c_genes.tar -C "${IGBLAST_DIR}/database" && \
    rm -f ncbi_human_c_genes.tar

ENV SEQKIT_VERSION=2.10.0
RUN if [ "$(uname -m)" = "x86_64" ]; then \
      BIN_ARCH=amd64; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
      BIN_ARCH=arm64; \
    fi && \
    wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_${BIN_ARCH}.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=/usr/local/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

ENV PATH=${IGBLAST_DIR}:$PATH

COPY build_ref.py /usr/local/src/

ENTRYPOINT ["python3", "/usr/local/src/build_ref.py"]