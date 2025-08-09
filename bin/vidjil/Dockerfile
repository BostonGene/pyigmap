FROM python:3.9-bullseye

RUN apt-get update && \
    apt-get -y install make pigz && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV SEQKIT_VERSION=2.10.0
RUN if [ "$(uname -m)" = "x86_64" ]; then \
      BIN_ARCH=amd64; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
      BIN_ARCH=arm64; \
    fi && \
    wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_${BIN_ARCH}.tar.gz -O seqkit.tar.gz && \
    tar -xvzf seqkit.tar.gz --one-top-level=/usr/local/bin/seqkit --strip-component 1 && \
    rm seqkit.tar.gz

ENV VIDJIL_VERSION=2024.02
RUN if [ "$(uname -m)" = "x86_64" ]; then \
        wget https://www.vidjil.org/releases/vidjil-algo-${VIDJIL_VERSION}_x86_64 -O /usr/local/bin/vidjil-algo && \
        chmod 755 /usr/local/bin/vidjil-algo; \
    elif [ "$(uname -m)" = "aarch64" ]; then \
        wget https://www.vidjil.org/releases/vidjil-algo-${VIDJIL_VERSION}.tar.gz && \
        tar -xvzf vidjil-algo-${VIDJIL_VERSION}.tar.gz && \
        cd vidjil-algo-${VIDJIL_VERSION} && make && cp vidjil-algo /usr/local/bin/vidjil-algo && \
        chmod 755 /usr/local/bin/vidjil-algo && \
        cd .. && rm -rf vidjil-algo-${VIDJIL_VERSION}*; \
    fi

COPY run.py /usr/local/src/

ENTRYPOINT ["python3.9", "/usr/local/src/run.py"]