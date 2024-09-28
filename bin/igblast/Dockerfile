FROM python:3.9-bullseye AS image

RUN apt-get update && apt-get -y install pigz --no-install-recommends

ENV IGBLAST_VERSION=1.22.0
RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz -O ncbi-igblast.tar.gz && \
    tar -xvzf ncbi-igblast.tar.gz --one-top-level=/usr/local/bin/ncbi-igblast --strip-component 1 && \
    rm ncbi-igblast.tar.gz

ENV SEQTK_VERSION=1.4
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz -O seqtk.tar.gz && \
    tar -xvzf seqtk.tar.gz --one-top-level=seqtk --strip-component 1 && \
    cd seqtk && make && mv seqtk /usr/local/bin && \
    cd .. && rm -r seqtk seqtk.tar.gz

ENV SEQKIT_VERSION=2.8.0
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=/usr/local/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

ENV IGBLAST_DIR=/usr/local/bin/ncbi-igblast
ENV PATH=${IGBLAST_DIR}:$PATH

COPY run.py /usr/local/src/

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/src/run.py"]

FROM image as build-ref

RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar -q -O ncbi_human_c_genes.tar && \
    mkdir ${IGBLAST_DIR}/database && tar -xvf ncbi_human_c_genes.tar -C ${IGBLAST_DIR}/database

COPY build_ref.py /usr/local/src/

ENTRYPOINT ["python3.9", "/usr/local/src/build_ref.py"]