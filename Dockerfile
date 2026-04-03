FROM rocker/tidyverse:latest

# Install system dependencies for nloptr, np, and OpenMP
RUN apt-get update && apt-get install -y --no-install-recommends \
    libnlopt-dev \
    libgsl-dev \
    libomp-dev \
    cmake \
    pkg-config \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "options(Ncpus = parallel::detectCores()); \
    install.packages(c('Rcpp', 'nloptr', 'doParallel', 'foreach'), \
    repos='https://cloud.r-project.org/')" \
    && R -e "pkgs <- c('Rcpp', 'nloptr', 'doParallel', 'foreach'); \
    inst <- installed.packages()[, 'Package']; \
    if (!all(pkgs %in% inst)) stop(paste('Missing:', pkgs[!pkgs %in% inst]))"
