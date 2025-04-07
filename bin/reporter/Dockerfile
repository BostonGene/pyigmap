FROM python:3.10-bullseye AS image

COPY requirements.txt .
RUN pip3 install --no-deps -r requirements.txt

COPY run.py logger.py utils.py viz.py /usr/local/src/

FROM image AS tool

ENTRYPOINT ["python3.10", "/usr/local/src/run.py"]