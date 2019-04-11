FROM rocker/r-ver:3.3.3

# File Author / Maintainer
MAINTAINER Toni Hermoso Pulido <toni.hermoso@crg.eu>

ARG BOWTIE_VERSION=1.2.1.1
ARG PSIPLOT_VERSION=2.3.0

# Install external dependencies 
RUN apt-get update -qq && apt-get install -y --no-install-recommends python curl libcurl4-openssl-dev libssl-dev libsqlite3-dev libxml2-dev qpdf git

# Install bowtie 
RUN cd /usr/local; curl --fail --silent --show-error --location --remote-name https://github.com/BenLangmead/bowtie/releases/download/v$BOWTIE_VERSION/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip
RUN cd /usr/local; unzip -d /usr/local bowtie-${BOWTIE_VERSION}-linux-x86_64.zip

RUN cd /usr/local; rm bowtie-${BOWTIE_VERSION}-linux-x86_64.zip

# Let's put in PATH
RUN cd /usr/local/bin; ln -s ../bowtie-${BOWTIE_VERSION}/bowtie* .

COPY deps.R /usr/local

RUN Rscript /usr/local/deps.R

# Psiplot
RUN cd /usr/local/; curl --fail --silent --show-error --location --remote-name https://github.com/kcha/psiplot/archive/v${PSIPLOT_VERSION}.tar.gz
RUN Rscript -e "install.packages( \"/usr/local/v${PSIPLOT_VERSION}.tar.gz\", repos = NULL )"
RUN rm /usr/local/v${PSIPLOT_VERSION}.tar.gz

# Install Vast-tools
RUN mkdir -p /usr/local/vast-tools
COPY vast-tools /usr/local/vast-tools
COPY install.R /usr/local/vast-tools
COPY lib /usr/local/vast-tools/lib
COPY bin /usr/local/vast-tools/bin
COPY R /usr/local/vast-tools/R

VOLUME /VASTDB

RUN cd /usr/local/vast-tools; ln -s /VASTDB .

RUN cd /usr/local/vast-tools; ./install.R

# Let's put in PATH
RUN cd /usr/local/bin; ln -s ../vast-tools/vast-tools .

# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*

# Shared mounting
VOLUME /share


