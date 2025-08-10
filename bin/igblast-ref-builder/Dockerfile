FROM python:3.13-bullseye

RUN apt-get update \
 && apt-get -y install --no-install-recommends wget pigz ca-certificates \
 && rm -rf /var/lib/apt/lists/*

ENV IGBLAST_VERSION=1.22.0
ENV IGBLAST_DIR=/usr/local/bin/ncbi-igblast
RUN if [ "$(uname -m)" = "x86_64" ]; then \
      wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz -O ncbi-igblast.tar.gz && \
      tar -xvzf ncbi-igblast.tar.gz --one-top-level=${IGBLAST_DIR} --strip-component 1 && \
      rm ncbi-igblast.tar.gz ; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
      wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-src.tar.gz -O ncbi-igblast.tar.gz && \
      tar -xvzf ncbi-igblast.tar.gz --one-top-level=ncbi-igblast --strip-component 1 && rm ncbi-igblast.tar.gz && \
      cd ncbi-igblast/c++ && ./configure && make -j$(nproc) && make install && \
      mkdir -p ${IGBLAST_DIR}/bin/ && cp /ncbi-igblast/c++/ReleaseMT/bin/igblastn ${IGBLAST_DIR}/bin/ && \
      chmod 755 ${IGBLAST_DIR}/bin/igblastn && \
      cd / && rm -rf ncbi-igblast* ; \
    fi && \
    wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar -q -O ncbi_human_c_genes.tar && \
    mkdir ${IGBLAST_DIR}/database && tar -xvf ncbi_human_c_genes.tar -C ${IGBLAST_DIR}/database

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