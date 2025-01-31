FROM python:3.8-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends --no-install-suggests\
    build-essential \
    libzmq3-dev \
    poppler-utils  \
    tesseract-ocr  \
    qpdf  \
    ffmpeg  \
    libsm6  \
    libxext6 \
    && \
    apt-get clean -qy && rm -rf /var/lib/apt/lists/*

COPY requirements-docker.txt /app/requirements.txt

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && pip install --no-cache-dir -r requirements.txt -v

COPY ./ /app/drug_server

ENV PYTHONPATH=/app/drug_server \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN useradd --create-home appuser
USER appuser
