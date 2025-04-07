FROM python:3.9-bullseye AS image

COPY requirements.txt .
RUN pip3 install --no-deps -r requirements.txt

COPY run.py airr.py pgen.py correct.py logger.py filter.py /usr/local/src/

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/src/run.py"]

FROM image AS build-ref

RUN apt-get update && apt-get install git

COPY build_ref.py /usr/local/src/

ENTRYPOINT ["python3.9", "/usr/local/src/build_ref.py"]