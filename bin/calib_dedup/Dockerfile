FROM python:3.9-bullseye AS image

ENV SOFT_DIR=/usr/local

RUN apt-get update && apt-get install -y cmake pigz jq

ENV CALIB_VERSION=0.3.7
RUN wget -q https://github.com/vpc-ccg/calib/archive/refs/tags/v${CALIB_VERSION}.tar.gz -O CALIB.tar.gz && \
    tar -xzvf CALIB.tar.gz --one-top-level=CALIB --strip-component 1 && \
    cd CALIB && \
    make && mv calib ${SOFT_DIR}/bin && \
    make -C consensus/ && mv consensus/calib_cons ${SOFT_DIR}/bin && cd .. && \
    rm CALIB.tar.gz && rm -r CALIB

COPY utils.py ${SOFT_DIR}/utils.py
COPY run.py ${SOFT_DIR}/run.py
COPY logger.py ${SOFT_DIR}/logger.py

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/run.py"]
