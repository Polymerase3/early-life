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
    gfortran \
    libgsl-dev \
    fonts-dejavu \
    xz-utils \
    fontconfig \
    wget \
    unzip \
  && rm -rf /var/lib/apt/lists/*

# Install Montserrat font (not available via apt in this base image)
RUN mkdir -p /usr/local/share/fonts/truetype/montserrat \
  && wget -O /usr/local/share/fonts/truetype/montserrat/Montserrat-VF.ttf "https://github.com/google/fonts/raw/main/ofl/montserrat/Montserrat%5Bwght%5D.ttf" \
  && fc-cache -f

COPY phiper_0.2.6.tar.gz /tmp/phiper_0.2.6.tar.gz

RUN R -q -e "install.packages(c('rlang','ggplot2','Cairo','openxlsx','dplyr','purrr','locfdr','withr','future','future.apply','plotly','htmlwidgets','tibble','duckdb','chk','dbplyr','vegan','RcppParallel','Rtsne'), repos='https://cloud.r-project.org')" \
  && R -q -e "install.packages('/tmp/phiper_0.2.6.tar.gz', repos=NULL, type='source')"

WORKDIR /work
ENTRYPOINT ["Rscript"]
CMD ["R/03-analysis.R"]
