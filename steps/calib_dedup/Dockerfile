FROM python:3.9-bullseye AS image

ENV SOFT_DIR=/usr/local

RUN apt-get update && apt-get install -y cmake pigz

COPY requirements.txt .
RUN pip3 install -r requirements.txt

ENV CALIB_VERSION=0.3.7
RUN wget -q https://github.com/vpc-ccg/calib/archive/refs/tags/v${CALIB_VERSION}.tar.gz -O CALIB.tar.gz && \
    tar -xzvf CALIB.tar.gz --one-top-level=CALIB --strip-component 1 && \
    cd CALIB && \
    make && mv calib ${SOFT_DIR}/bin && \
    make -C consensus/ && mv consensus/calib_cons ${SOFT_DIR}/bin && cd .. && \
    rm CALIB.tar.gz && rm -r CALIB

ENV SEQKIT_VERSION=2.8.0
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=${SOFT_DIR}/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

COPY run.py ${SOFT_DIR}/run.py
COPY barcode.py ${SOFT_DIR}/barcode.py
COPY logger.py ${SOFT_DIR}/logger.py

FROM image AS tool
# Uncomment it if you need to run an executable container
ENTRYPOINT ["python3.9", "/usr/local/run.py"]