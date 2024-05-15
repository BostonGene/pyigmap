FROM python:3.9-bullseye AS image

ENV SOFT_DIR=/usr/local

RUN apt-get update && apt-get -y install wget pigz

COPY requirements.txt .
RUN pip3 install -r requirements.txt

ENV IGBLAST_VERSION=1.22.0
RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz -O ncbi-igblast.tar.gz && \
    tar -xvzf ncbi-igblast.tar.gz --one-top-level=${SOFT_DIR}/bin/ncbi-igblast --strip-component 1 && \
    rm ncbi-igblast.tar.gz

ENV SEQTK_VERSION=1.4
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz -O seqtk.tar.gz && \
    tar -xvzf seqtk.tar.gz --one-top-level=seqtk --strip-component 1 && \
    cd seqtk && make && mv seqtk ${SOFT_DIR}/bin && \
    cd .. && rm -r seqtk seqtk.tar.gz

ENV SEQKIT_VERSION=2.8.0
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=${SOFT_DIR}/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

ENV IGBLAST_DIR=${SOFT_DIR}/bin/ncbi-igblast
ENV PATH=${IGBLAST_DIR}:$PATH

COPY run.py ${SOFT_DIR}/run.py
COPY utils.py ${SOFT_DIR}/utils.py
COPY logger.py ${SOFT_DIR}/logger.py

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/run.py"]

FROM image as ref

RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/ncbi_human_c_genes.tar -q -O ncbi_human_c_genes.tar && \
    tar -xvf ncbi_human_c_genes.tar -C ${IGBLAST_DIR}/database \

COPY build_ref.py ${SOFT_DIR}/build_ref.py

ENTRYPOINT ["python3.9", "/usr/local/build_ref.py"]