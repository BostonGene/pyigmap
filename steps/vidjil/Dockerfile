FROM python:3.9-bullseye

RUN apt-get update && apt-get -y install wget pigz parallel

ENV VIDJIL_VERSION=2024.02
RUN wget https://www.vidjil.org/releases/vidjil-algo-${VIDJIL_VERSION}_x86_64 -O vidjil-algo && \
    chmod 755 vidjil-algo && \
    mv vidjil-algo /usr/local/bin/

COPY run.py .
COPY logger.py .

ENTRYPOINT ["python3.9", "run.py"]