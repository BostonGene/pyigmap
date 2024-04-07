FROM python:3.9-bullseye

ENV SOFT_DIR=/usr/local/bin

ENV FASTP_VERSION 0.23.4
RUN wget -q http://opengene.org/fastp/fastp.${FASTP_VERSION} -O ${SOFT_DIR}/fastp && \
    chmod a+x ${SOFT_DIR}/fastp

COPY run.py ${SOFT_DIR}/run.py

# Uncomment it if you need to run an executable container
# ENTRYPOINT ["python3.9", "${SOFT_DIR}/run.py"]