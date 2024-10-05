FROM python:3.9-bullseye AS image

COPY requirements.txt .
RUN pip3 install --no-deps -r requirements.txt

COPY run.py airr.py pgen.py correct.py logger.py filter.py /usr/local/src/

FROM image AS tool

ENTRYPOINT ["python3.9", "/usr/local/src/run.py"]