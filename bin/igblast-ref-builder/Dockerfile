FROM ubuntu:24.04

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      ca-certificates \
      wget \
      python3 \
      tar \
      xz-utils && \
    rm -rf /var/lib/apt/lists/*

ENV IGBLAST_VERSION=1.22.0
ENV IGBLAST_DIR=/usr/local/bin/ncbi-igblast
ARG TARGETARCH

RUN set -eux; \
    mkdir -p "${IGBLAST_DIR}"; \
    if [ "${TARGETARCH:-amd64}" = "amd64" ]; then \
        wget -q -O /tmp/igblast.tgz "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz" && \
        tar -xzf /tmp/igblast.tgz -C "${IGBLAST_DIR}" --strip-component=1 && \
        rm -f /tmp/igblast.tgz; \
    elif [ "${TARGETARCH}" = "arm64" ]; then \
        apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
            build-essential zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
            libboost-all-dev libssl-dev libexpat1-dev libxml2-dev libprotobuf-dev protobuf-compiler \
            libsasl2-dev libhdf5-dev liblapack-dev libblas-dev libgsl-dev cmake git pkg-config cpio && \
            rm -rf /var/lib/apt/lists/* && \
        wget -q -O /tmp/igblast-src.tgz "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-src.tar.gz" && \
        mkdir -p /tmp/igblast-src && \
        tar -xzf /tmp/igblast-src.tgz -C /tmp/igblast-src --strip-component=1 && rm -f /tmp/igblast-src.tgz && \
        cd /tmp/igblast-src/c++; \
        [ -d /lib64 ] || ln -s /lib /lib64; \
        ./configure && make -j"$(nproc)" && make install && \
        mkdir -p "${IGBLAST_DIR}/bin" && \
        cp /tmp/igblast-src/c++/ReleaseMT/bin/makeblastdb "${IGBLAST_DIR}/bin/" && chmod 755 "${IGBLAST_DIR}/bin/makeblastdb" && \
        cp /tmp/igblast-src/c++/ReleaseMT/bin/edit_imgt_file.pl "${IGBLAST_DIR}/bin/" && \
        cp -r /tmp/igblast-src/c++/src/app/igblast/internal_data "${IGBLAST_DIR}/" && \
        cp -r /tmp/igblast-src/c++/src/app/igblast/optional_file "${IGBLAST_DIR}/" && \
        cd / && rm -rf /tmp/igblast-src; \
    else \
        echo "Unsupported TARGETARCH=${TARGETARCH}"; exit 1; \
    fi && \
    wget -q -O /tmp/ncbi_human_c_genes.tar "https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar" && \
    mkdir -p "${IGBLAST_DIR}/database" && \
    tar -xf /tmp/ncbi_human_c_genes.tar -C "${IGBLAST_DIR}/database" && \
    rm -f /tmp/ncbi_human_c_genes.tar

ENV SEQKIT_VERSION=2.10.0
RUN set -eux; \
    if [ "${TARGETARCH:-amd64}" = "amd64" ]; then BIN_ARCH=amd64; else BIN_ARCH=arm64; fi; \
    wget -q -O /tmp/seqkit.tgz "https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_${BIN_ARCH}.tar.gz" && \
    tar -xzf /tmp/seqkit.tgz -C /usr/local/bin seqkit && \
    chmod +x /usr/local/bin/seqkit && \
    rm -f /tmp/seqkit.tgz

ENV PATH="${IGBLAST_DIR}/bin:${IGBLAST_DIR}:$PATH"

COPY build_ref.py /usr/local/src/
ENTRYPOINT ["python3", "/usr/local/src/build_ref.py"]