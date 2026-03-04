FROM python:3.11-slim-bookworm

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PPPRPP_SERVER_NAME=0.0.0.0 \
    PPPRPP_SERVER_PORT=7860 \
    PPPRPP_INBROWSER=0 \
    GRADIO_ANALYTICS_ENABLED=False

RUN apt-get update && apt-get install -y --no-install-recommends \
    hmmer \
    openjdk-17-jre-headless \
    perl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt ./
RUN python -m pip install --upgrade pip && pip install -r requirements.txt

COPY . .

EXPOSE 7860

CMD ["python", "app.py"]
