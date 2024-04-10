FROM python:3.9-bullseye

ENV SOFT_DIR=/usr/local/bin

RUN apt-get update && apt-get install -y jq pigz

ENV FFQ_VERSION=0.3.1
RUN pip install ffq==$FFQ_VERSION

COPY run.py ${SOFT_DIR}/run.py