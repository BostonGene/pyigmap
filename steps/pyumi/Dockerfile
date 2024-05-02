FROM python:3.9-bullseye

ENV SOFT_DIR=/usr/local

RUN apt-get update && apt-get install -y pigz

COPY requirements.txt .
RUN pip3 install -r requirements.txt

ENV SEQKIT_VERSION=2.8.1
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=${SOFT_DIR}/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

COPY utils.py ${SOFT_DIR}/utils.py
COPY run.py ${SOFT_DIR}/run.py
COPY extract.py ${SOFT_DIR}/extract.py
COPY pattern.py ${SOFT_DIR}/pattern.py
COPY logger.py ${SOFT_DIR}/logger.py

# Uncomment it if you need to run an executable container
ENTRYPOINT ["python3.9", "/usr/local/run.py"]