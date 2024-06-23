FROM python:3.9-bullseye AS image

ENV SOFT_DIR=/usr/local

RUN apt-get update

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY utils.py ${SOFT_DIR}/utils.py
COPY run.py ${SOFT_DIR}/run.py
COPY airr.py ${SOFT_DIR}/airr.py
COPY pgen.py ${SOFT_DIR}/pgen.py
COPY correct.py ${SOFT_DIR}/correct.py
COPY logger.py ${SOFT_DIR}/logger.py
COPY filter.py ${SOFT_DIR}/filter.py

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/run.py"]