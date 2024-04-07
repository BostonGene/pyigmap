FROM python:3.9-bullseye

ENV SOFT_DIR=/usr/local/bin

RUN apt-get update && apt-get -y install wget pigz parallel

ENV VIDJIL_VERSION=2024.02
RUN wget https://www.vidjil.org/releases/vidjil-algo-${VIDJIL_VERSION}_x86_64 -O vidjil-algo && \
    chmod 755 vidjil-algo && \
    mv vidjil-algo $SOFT_DIR

COPY run.py ${SOFT_DIR}/run.py
COPY logger.py ${SOFT_DIR}/logger.py

# Uncomment it if you need to run an executable container
# ENTRYPOINT ["python3.9", "${SOFT_DIR}/run.py"]