FROM python:3.9-bullseye AS image

ENV FASTP_VERSION=0.23.4
RUN wget -q http://opengene.org/fastp/fastp.${FASTP_VERSION} -O /usr/local/bin/fastp && \
    chmod a+x /usr/local/bin/fastp

COPY run.py mock_merge.py logger.py /usr/local/src/

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/src/run.py"]