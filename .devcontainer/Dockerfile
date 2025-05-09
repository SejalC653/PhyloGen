FROM python:3.9-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app

# Install system dependencies and bioinformatics tools
RUN apt-get update && apt-get install -y \
    mafft \
    muscle \
    clustalw \
    fasttree \
    build-essential \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --upgrade pip && pip install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

# Expose Streamlit port
EXPOSE 8501

# Default command to run the Streamlit app
CMD ["streamlit", "run", "PG.py", "--server.port=8501", "--server.address=0.0.0.0"]
