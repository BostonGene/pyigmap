FROM python:3.9-bullseye

LABEL image.authors="nikita.syzrantsev@bostongene.com"

RUN apt-get update && apt-get -y install pigz --no-install-recommends

ENV IGBLAST_VERSION=1.22.0
ENV IGBLAST_DIR=/usr/local/bin/ncbi-igblast
RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz -O ncbi-igblast.tar.gz && \
    tar -xvzf ncbi-igblast.tar.gz --one-top-level=/usr/local/bin/ncbi-igblast --strip-component 1 && \
    rm ncbi-igblast.tar.gz && \
    wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar -q -O ncbi_human_c_genes.tar && \
    mkdir ${IGBLAST_DIR}/database && tar -xvf ncbi_human_c_genes.tar -C ${IGBLAST_DIR}/database

ENV SEQKIT_VERSION=2.8.0
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=/usr/local/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

ENV PATH=${IGBLAST_DIR}:$PATH

COPY build_ref.py /usr/local/src/

ENTRYPOINT ["python3.9", "/usr/local/src/build_ref.py"]