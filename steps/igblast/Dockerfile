FROM python:3.9-bullseye

RUN apt-get update && apt-get -y install wget pigz

COPY requirements.txt .
RUN pip3 install -r requirements.txt

ENV IGBLAST_VERSION=1.22.0
RUN wget https://ftp.ncbi.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz -O ncbi-igblast.tar.gz && \
    tar -xvzf ncbi-igblast.tar.gz --one-top-level=ncbi-igblast --strip-component 1 && \
    rm ncbi-igblast.tar.gz

ENV SEQTK_VERSION=1.4
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v${SEQTK_VERSION}.tar.gz -O seqtk.tar.gz && \
    tar -xvzf seqtk.tar.gz --one-top-level=seqtk --strip-component 1 && \
    cd seqtk && make && mv seqtk /usr/local/bin/ && \
    cd .. && rm -r seqtk seqtk.tar.gz

ENV SEQKIT_VERSION=2.8.0
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=seqkit --strip-component 1 && \
    mv seqkit /usr/local/bin/ && \
    rm seqkit.tar.gz

WORKDIR /

COPY run.py .
COPY logger.py .

ENTRYPOINT ["python3.9", "run.py"]