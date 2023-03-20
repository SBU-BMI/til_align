# Dockerfile containing system libraries and R packages for TIL-align.
#
# This is separated into a separate Dockerfile because it can take 5 hours to build
# for ARM64 (using qemu emulation on an AMD64 build platform). This part of the Docker
# image remains constant, so it is in its own Dockerfile and is used as a base image
# for the TIL-align Docker image.

FROM rocker/r-ver:4.2.1

RUN apt-get update \
    && apt-get install -yq --no-install-recommends \
        cmake \
        libcurl4-openssl-dev \
        libproj15 \
        libgdal-dev \
        pandoc \
        texlive \
        lmodern \
        texlive-latex-extra \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus 4 \
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
