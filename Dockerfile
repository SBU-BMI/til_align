# Bootstrap a previous built container with all dependencies.
# Build this container with 'docker build -f base.Dockerfile .'
FROM kaczmarj/tilalign:9fd8a17bd13a36a8d59a1ff25b9f36d4cb75beb9

# Add dependency from https://github.com/SBU-BMI/til_align/pull/30
RUN install2.r --error --skipinstalled reshape2 \
    && rm -rf /tmp/downloaded_packages \
    && strip /usr/local/lib/R/site-library/*/libs/*.so

WORKDIR /code
COPY scripts_and_docker/ .
RUN chmod +x *.sh *.R *.rmd
ENV PATH="/code/":$PATH
RUN mkdir -p /data

# Copy this to a writable directory due to
# https://github.com/rstudio/rmarkdown/issues/1975
RUN mkdir -p /tmp/rmarkdowndir/ \
    && chmod a+rwx /tmp/rmarkdowndir \
    && cp /code/Descriptive_Statistics.rmd /tmp/rmarkdowndir/

CMD ["/bin/bash"]
