FROM rocker/tidyverse:latest

RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*
RUN R -e "install.packages( \
        c('Rcpp', 'np', 'nloptr', 'stats'), \
        repos='https://cloud.r-project.org/', \
        clean = TRUE)" \
    && rm -rf /tmp/downloaded_packages

