FROM python:3.9-bullseye as not_exec

ENV SOFT_DIR=/usr/local

RUN apt-get update && apt-get -y install wget pigz parallel

ENV VIDJIL_VERSION=2024.02
RUN wget https://www.vidjil.org/releases/vidjil-algo-${VIDJIL_VERSION}_x86_64 -O ${SOFT_DIR}/bin/vidjil-algo && \
    chmod 755 ${SOFT_DIR}/bin/vidjil-algo

COPY run.py ${SOFT_DIR}/run.py
COPY logger.py ${SOFT_DIR}/logger.py

FROM not_exec as exec
# Uncomment it if you need to run an executable container
ENTRYPOINT ["python3.9", "/usr/local/run.py"]