FROM python:3.9-bullseye as image

RUN apt-get update && apt-get install -y jq pigz

ENV FFQ_VERSION=0.3.1
RUN pip install ffq==$FFQ_VERSION
