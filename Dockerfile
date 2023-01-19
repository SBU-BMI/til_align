FROM rocker/r-ver:4.2.1

#RUN mkdir /TILalignment
#RUN mkdir /TILalignment/data
#RUN mkdir /TILalignment/data/tilPreds
#RUN mkdir /TILalignment/data/cancPreds
#RUN mkdir /TILalignment/data/outputs

#RUN mkdir /Downstream
#RUN mkdir /Downstream/data

RUN mkdir /code
RUN mkdir /data

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
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

COPY callAlign_withDocker.sh /code/callAlign.sh
COPY commandLineAlign.R /code/commandLineAlign.R
COPY callAnalytics_withDocker.sh  /code/callAnalytics.sh
COPY renderWrapper.R /code/renderWrapper.R
COPY Descriptive_Statistics.rmd /code/Descriptive_Statistics.rmd

#RUN chmod 0755 /TILalignment/*
#RUN chmod 0755 /Downstream/*
RUN chmod 0755 /code/*

ENV PATH="./":$PATH
#RUN echo ${PATH}
WORKDIR "/code/"
#RUN ls -l .

CMD ["/bin/bash"]
