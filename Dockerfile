FROM rocker/tidyverse:latest

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libnlopt-dev \
    libgsl-dev \
    libomp-dev \
    cmake \
    pkg-config \
    zlib1g-dev \
    openssh-client \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "pkgs <- c('ranger', 'np', 'Rcpp', 'nloptr', 'doFuture', 'foreach'); \
    options(Ncpus = parallel::detectCores()); \
    install.packages(pkgs, \
    repos='https://cloud.r-project.org/')" \
    && R -e "inst <- installed.packages()[, 'Package']; \
    if (!all(pkgs %in% inst)) stop(paste('Missing:', pkgs[!pkgs %in% inst]))"
