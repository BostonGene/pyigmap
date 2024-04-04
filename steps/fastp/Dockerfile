FROM python:3.9-bullseye

ENV FASTP_VERSION 0.23.4
RUN wget -q http://opengene.org/fastp/fastp.${FASTP_VERSION} -O /usr/local/bin/fastp && \
    chmod a+x /usr/local/bin/fastp

WORKDIR /

COPY run.py .

ENTRYPOINT ["python3.9", "run.py"]