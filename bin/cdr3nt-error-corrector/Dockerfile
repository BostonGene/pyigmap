FROM python:3.13-bullseye AS image

COPY requirements.txt .
RUN python3 -m venv venv && \
    venv/bin/pip3 install --no-deps -r requirements.txt

COPY run.py airr.py pgen.py correct.py logger.py filter.py /usr/local/src/

# ENTRYPOINT ["python3", "/usr/local/src/run.py"]