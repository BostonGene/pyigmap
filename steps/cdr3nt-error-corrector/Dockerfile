FROM python:3.9-bullseye

RUN apt-get update

COPY requirements.txt .
RUN pip3 install -r requirements.txt

COPY utils.py .
COPY run.py .
COPY airr.py .
COPY pgen.py .
COPY correct.py .
COPY logger.py .
COPY filter.py .

ENTRYPOINT ["python3.9", "run.py"]