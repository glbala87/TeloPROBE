FROM python:3.12-slim AS base

LABEL maintainer="TeloPROBE Developers"
LABEL description="TeloPROBE: Probabilistic Robust Observation of Boundary Endpoints for telomere length estimation"
LABEL version="2.0.0"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    wget \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install minimap2 and samtools from source (specific versions for reproducibility)
ARG MINIMAP2_VERSION=2.28
ARG SAMTOOLS_VERSION=1.21

RUN wget -q https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
    && tar -xjf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
    && mv minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 /usr/local/bin/ \
    && rm -rf minimap2-${MINIMAP2_VERSION}*

RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local --without-curses \
    && make -j$(nproc) \
    && make install \
    && cd .. && rm -rf samtools-${SAMTOOLS_VERSION}*

# Install TeloPROBE
WORKDIR /opt/teloprobe

COPY pyproject.toml .
COPY src/ src/
COPY data/ data/
COPY config/ config/
COPY workflow/ workflow/

RUN pip install --no-cache-dir -e ".[dev]"

# Verify installation
RUN teloprobe --version \
    && minimap2 --version \
    && samtools --version | head -1 \
    && python -c "import teloprobe; print(f'TeloPROBE {teloprobe.__version__} ready')"

# Default working directory for user data
WORKDIR /data

ENTRYPOINT ["teloprobe"]
CMD ["--help"]
