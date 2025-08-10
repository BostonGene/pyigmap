FROM python:3.13-bullseye

ENV VENV=/opt/venv
RUN python3 -m venv "$VENV"
ENV PATH="$VENV/bin:$PATH"

COPY requirements.txt .
RUN pip install --upgrade pip && \
    pip install --no-deps --no-cache-dir -r requirements.txt

COPY run.py airr.py pgen.py correct.py logger.py filter.py /usr/local/src/

# ENTRYPOINT ["python3", "/usr/local/src/run.py"]