FROM rocker/tidyverse:latest

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libnlopt-dev \
    libgsl-dev \
    libomp-dev \
    libblas-dev \
    liblapack-dev \
    cmake \
    pkg-config \
    zlib1g-dev \
    openssh-client \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "pkgs <- c('ranger', 'np', 'Rcpp', 'nloptr', 'doFuture', 'foreach'); \
    options(Ncpus = 2); \
    install.packages(pkgs, repos='https://cloud.r-project.org/'); \
    inst <- installed.packages()[, 'Package']; \
    if (!all(pkgs %in% inst)) { \
        cat('\n\nFAILED INSTALLATION:\n'); \
        print(setdiff(pkgs, inst)); \
        quit(status = 1); \
    }"
