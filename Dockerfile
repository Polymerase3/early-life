FROM rocker/r-ver:4.3.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff5-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libxt-dev \
    libzip-dev \
    fonts-dejavu \
    fonts-montserrat \
  && rm -rf /var/lib/apt/lists/*

COPY phiper_0.2.4.tar.gz /tmp/phiper_0.2.4.tar.gz

RUN R -q -e "install.packages(c('rlang','ggplot2','Cairo','openxlsx','dplyr','purrr','locfdr','withr','future','plotly','htmlwidgets','tibble','duckdb'), repos='https://cloud.r-project.org')" \
  && R -q -e "install.packages('/tmp/phiper_0.2.4.tar.gz', repos=NULL, type='source')"

WORKDIR /work
ENTRYPOINT ["Rscript"]
CMD ["R/02-analysis.R"]
