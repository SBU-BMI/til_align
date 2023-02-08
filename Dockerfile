FROM rocker/r-ver:4.2.1

RUN apt-get update \
    && apt-get install -yq --no-install-recommends \
        libcurl4-openssl-dev \
        libproj15 \
        libgdal-dev \
        pandoc \
        texlive \
        lmodern \
        texlive-latex-extra \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled \
        abind \
        curl \
        plyr \
        dplyr \
        ggplot2 \
        ggpubr \
        survival \
        survminer \
        readr \
        raster \
        terra \
        rmarkdown \
        rjson \
    && rm -rf /tmp/downloaded_packages \
    && strip /usr/local/lib/R/site-library/*/libs/*.so

WORKDIR /code
COPY scripts_and_docker/ .
RUN chmod +x *.sh *.R *.rmd
ENV PATH="/code/":$PATH
RUN mkdir -p /data

# Copy this to a writable directory due to
# https://github.com/rstudio/rmarkdown/issues/1975
RUN mkdir /tmp/rmarkdowndir/ \
    && chmod a+rwx /tmp/rmarkdowndir \
    && cp /code/Descriptive_Statistics.rmd /tmp/rmarkdowndir/

CMD ["/bin/bash"]
