# Dockerfile
FROM python:3.13-slim

WORKDIR /app

COPY wheelhouse/mdaf3-*.whl /tmp/

RUN pip install --no-cache-dir /tmp/mdaf3-*.whl
