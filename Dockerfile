FROM python:3.11-slim

LABEL maintainer="DR JM NIPOK <nipok@njit.edu>"
LABEL description="SCT Cosmology Series — CAR framework environment"
LABEL version="1.7"

# System dependencies for CAMB (Fortran), scientific computing
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        gcc \
        g++ \
        make \
        cmake \
        git \
        wget \
        libopenmpi-dev \
        openmpi-bin \
        liblapack-dev \
        libblas-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /sct

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir \
        numpy scipy matplotlib astropy h5py pyyaml getdist \
    && pip install --no-cache-dir camb \
    && rm -rf /root/.cache/pip

# Copy repository
COPY . .

# Verify the core calculator runs correctly
RUN python sct_core.py && echo "sct_core.py verification: PASS"

# Create output directory
RUN mkdir -p /output/chains /output/figures /output/data

VOLUME ["/output"]
ENV OUTPUT_DIR=/output

# Default: run the CAR calculator; override with docker run ... bash chains/run_chains.sh
CMD ["python", "sct_core.py"]
