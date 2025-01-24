FROM rocker/tidyverse:3.6.3

LABEL maintainer="Toni Hermoso Pulido <toni.hermoso@crg.eu>"

ARG BOWTIE_VERSION=1.2.1.1
ARG GETOPT_VERSION=1.20.3
ARG OPTPARSE_VERSION=1.7.3
ARG PSIPLOT_VERSION=2.3.0

# Install external dependencies 
RUN apt-get update --allow-releaseinfo-change && apt-get upgrade -y && apt-get install -y --no-install-recommends python curl libcurl4-openssl-dev libssl-dev libsqlite3-dev libxml2-dev qpdf git
RUN apt-get install -y build-essential libharfbuzz-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfribidi-dev

# Install bowtie 
RUN cd /usr/local; curl --fail --silent --show-error --location --remote-name https://github.com/BenLangmead/bowtie/releases/download/v$BOWTIE_VERSION/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip
RUN cd /usr/local; unzip -d /usr/local bowtie-${BOWTIE_VERSION}-linux-x86_64.zip

RUN cd /usr/local; rm bowtie-${BOWTIE_VERSION}-linux-x86_64.zip

# Let's put in PATH
RUN cd /usr/local/bin; ln -s ../bowtie-${BOWTIE_VERSION}/bowtie* .

# Not used anymore
# COPY deps.R /usr/local
# RUN Rscript /usr/local/deps.R > /tmp/deps.log

# Github packages
RUN cd /usr/local/; curl --fail --silent --show-error --location --remote-name https://github.com/trevorld/r-getopt/archive/v${GETOPT_VERSION}.tar.gz
RUN cd /usr/local/; curl --fail --silent --show-error --location --remote-name https://github.com/trevorld/r-optparse/archive/v${OPTPARSE_VERSION}.tar.gz
RUN cd /usr/local/; curl --fail --silent --show-error --location --remote-name https://github.com/kcha/psiplot/archive/v${PSIPLOT_VERSION}.tar.gz
RUN Rscript -e "install.packages( \"/usr/local/v${GETOPT_VERSION}.tar.gz\", repos = NULL )"
RUN Rscript -e "install.packages( \"/usr/local/v${OPTPARSE_VERSION}.tar.gz\", repos = NULL )"
RUN Rscript -e "install.packages( \"/usr/local/v${PSIPLOT_VERSION}.tar.gz\", repos = NULL )"
RUN rm /usr/local/v${PSIPLOT_VERSION}.tar.gz; rm /usr/local/v${OPTPARSE_VERSION}.tar.gz; rm /usr/local/v${GETOPT_VERSION}.tar.gz

# Install Vast-tools
RUN mkdir -p /usr/local/vast-tools
COPY vast-tools /usr/local/vast-tools
COPY install.R /usr/local/vast-tools
COPY lib /usr/local/vast-tools/lib
COPY bin /usr/local/vast-tools/bin
COPY R /usr/local/vast-tools/R
COPY VERSION /usr/local/vast-tools/VERSION

VOLUME /VASTDB

RUN cd /usr/local/vast-tools; ln -s /VASTDB .

RUN cd /usr/local/vast-tools; ./install.R --quiet

# Let's put in PATH
RUN cd /usr/local/bin; ln -s ../vast-tools/vast-tools .

# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*

# Shared mounting
VOLUME /share


