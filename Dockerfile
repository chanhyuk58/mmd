FROM rocker/verse:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    libnlopt-dev \
    cmake \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "options(Ncpus = parallel::detectCores()); \
    install.packages(c('future', 'furrr'), repos='https://cloud.r-project.org/')"
