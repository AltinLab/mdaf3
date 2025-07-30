# Dockerfile
FROM python:3.13-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
       ca-certificates \
       git \
       build-essential \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY wheelhouse/mdaf3-*.whl /tmp/

RUN pip install --no-cache-dir /tmp/mdaf3-*.whl
