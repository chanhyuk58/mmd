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
RUN R -e "options(Ncpus = parallel::detectCores()); \
    install.packages(c('np', 'Rcpp', 'nloptr', 'doFuture', 'foreach'), \
    repos='https://cloud.r-project.org/')" \
    && R -e "pkgs <- c('np', 'Rcpp', 'nloptr', 'doFuture', 'foreach'); \
    inst <- installed.packages()[, 'Package']; \
    if (!all(pkgs %in% inst)) stop(paste('Missing:', pkgs[!pkgs %in% inst]))"
