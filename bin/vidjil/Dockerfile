FROM python:3.9-bullseye AS image

RUN apt-get update && \
    apt-get -y install parallel && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV VIDJIL_VERSION=2024.02
RUN wget https://www.vidjil.org/releases/vidjil-algo-${VIDJIL_VERSION}_x86_64 -O /usr/local/bin/vidjil-algo && \
    chmod 755 /usr/local/bin/vidjil-algo

COPY run.py /usr/local/src/

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/src/run.py"]