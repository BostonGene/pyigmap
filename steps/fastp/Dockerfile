FROM python:3.9-bullseye AS image

ENV SOFT_DIR=/usr/local

ENV FASTP_VERSION=0.23.4
RUN wget -q http://opengene.org/fastp/fastp.${FASTP_VERSION} -O ${SOFT_DIR}/bin/fastp && \
    chmod a+x ${SOFT_DIR}/bin/fastp

COPY run.py ${SOFT_DIR}/run.py

FROM image AS tool
# Uncomment it if you need to run an executable container
ENTRYPOINT ["python3.9", "/usr/local/run.py"]