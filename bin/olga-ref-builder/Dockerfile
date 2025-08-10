FROM python:3.13-bullseye

RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
    git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY build_ref.py /usr/local/src/

ENTRYPOINT ["python3", "/usr/local/src/build_ref.py"]