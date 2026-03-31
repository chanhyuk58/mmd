FROM rocker/tidyverse:latest

RUN apt-get update && apt-get install -y \
    libnlopt-dev \
    cmake \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages( \
    c('Rcpp', 'np', 'nloptr', 'txtplot'), repos='https://cloud.r-project.org/')" \
    && R -e "if (!all(c('Rcpp', 'np', 'nloptr', 'txtplot') %in% installed.packages())) \
    stop('Package installation failed!')"

