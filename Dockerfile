FROM rocker/tidyverse:latest

RUN apt-get update && apt-get install -y \
    libnlopt-dev \
    cmake \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages( \
        c('Rcpp', 'np', 'nloptr'), \
        repos='https://cloud.r-project.org/', \
        clean = TRUE)" \
    && rm -rf /tmp/downloaded_packages

