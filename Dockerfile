FROM ghcr.io/astral-sh/uv:0.5-python3.12-bookworm-slim

# Set environment variables
ENV UV_COMPILE_BYTECODE=1


RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y git libcurl4-openssl-dev procps

RUN mkdir -p /bgdatacache/ /home/user/.config/bbglab/

ENV BGDATA_LOCAL="/bgdatacache"
ENV BBGLAB_HOME="/home/user/.config/bbglab/"

RUN pip install bgdata bgreference

RUN echo "# Version of the bgdata config file\n\
version = 2\n\
\n\
# The default local folder where you want to store the data packages\n\
local_repository = \"/bgdatacache\"\n\
\n\
# The remote URL from where do you want to download the data packages\n\
remote_repository = \"https://bbglab.irbbarcelona.org/bgdata\"\n\
\n\
# If you want to force bgdata to work only locally\n\
# offline = True\n\
\n\
# Cache repositories\n\
[cache_repositories]" > /home/user/.config/bbglab/bgdatav2.conf

RUN bgdata get datasets/genomereference/hg38 && bgdata get datasets/genomereference/hg19
RUN bgdata get datasets/genomereference/mm39 && bgdata get datasets/genomereference/mm10

RUN python3 -c "from bgreference import hg38, hg19; hg38(1, 1300000, 3000); hg19(1, 1300000, 3000)"
RUN python3 -c "from bgreference import mm39, mm10; mm39(1, 1300000, 3000); mm10(1, 1300000, 3000)"

RUN chmod -R 777 /bgdatacache/*

# Install dependencies and build the project
RUN mkdir -p /root/.cache/uv \
    && git clone https://github.com/bbglab/omega.git \
    && cd omega \
    && git checkout main-dev \
    && uv sync --frozen --no-dev \
    && uv build \
    && pip install dist/*.tar.gz \
    # Clean up unnecessary files
    && rm -rf /root/.cache/uv

RUN omega --help

