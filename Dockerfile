# Dockerfile for SCT Cosmology Series
# One-command reproducibility for CAR analysis

FROM python:3.9-slim

LABEL maintainer="DR JM NIPOK"
LABEL description="SCT Cosmology Series - Codified Acoustic Relation (CAR) Framework"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    gfortran \
    make \
    cmake \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /sct

# Copy requirements first for better caching
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the entire repository
COPY . .

# Set environment variables
ENV PYTHONPATH="/sct:${PYTHONPATH}"

# Run the CAR calculator to verify installation
RUN python sct_core.py

# Default command
CMD ["bash"]