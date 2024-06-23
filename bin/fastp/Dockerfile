FROM python:3.9-bullseye AS image

ENV SOFT_DIR=/usr/local

RUN apt-get update && apt-get install -y pigz

COPY requirements.txt .
RUN pip3 install -r requirements.txt

ENV FASTP_VERSION=0.23.4
RUN wget -q http://opengene.org/fastp/fastp.${FASTP_VERSION} -O ${SOFT_DIR}/bin/fastp && \
    chmod a+x ${SOFT_DIR}/bin/fastp

COPY run.py ${SOFT_DIR}/run.py
COPY mock_merge.py ${SOFT_DIR}/mock_merge.py
COPY logger.py ${SOFT_DIR}/logger.py

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/run.py"]