FROM python:3.9-bullseye

RUN apt-get update && apt-get install -y cmake pigz

COPY requirements.txt .
RUN pip3 install -r requirements.txt

ENV CALIB_VERSION=0.3.7
RUN wget -q https://github.com/vpc-ccg/calib/archive/refs/tags/v${CALIB_VERSION}.tar.gz -O CALIB.tar.gz && \
    tar -xzvf CALIB.tar.gz --one-top-level=CALIB --strip-component 1 && \
    cd CALIB && \
    make && mv calib /usr/local/bin/ && \
    make -C consensus/ && mv consensus/calib_cons /usr/local/bin/ && cd .. && \
    rm CALIB.tar.gz && rm -r CALIB

ENV SEQKIT_VERSION=2.8.0
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=seqkit --strip-component 1 && \
    mv seqkit /usr/local/bin/ && \
    rm seqkit.tar.gz

WORKDIR /

COPY run.py .
COPY barcode.py .
COPY logger.py .

ENTRYPOINT ["python3.9", "run.py"]
