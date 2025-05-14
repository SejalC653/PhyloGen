FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender-dev \
    cmake \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN wget http://www.microbesonline.org/fasttree/FastTree -O /usr/local/bin/FastTree && \
    chmod +x /usr/local/bin/FastTree

RUN apt-get update && \
    apt-get install -y mafft

RUN apt-get update && \
    apt-get install -y muscle

RUN apt-get update && \
    apt-get install -y clustalo

WORKDIR /app

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "main.py"]
