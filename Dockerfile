# Use specific R version for reproducibility
FROM rocker/r-ver:4.0.3

# Maintainer information
LABEL maintainer="Toni Hermoso Pulido <toni.hermoso@crg.eu>"

# Define arguments for tool versions
ARG BOWTIE_VERSION=1.2.1.1
ARG PSIPLOT_VERSION=2.3.0

# Install dependencies in a single RUN command for efficiency
RUN apt-get update -qq && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    python curl libcurl4-openssl-dev libssl-dev libsqlite3-dev \
    libxml2-dev qpdf git build-essential libharfbuzz-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev libfribidi-dev wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install bowtie in a single step to reduce layers
RUN curl -fsSL https://github.com/BenLangmead/bowtie/releases/download/v${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip -o /tmp/bowtie.zip && \
    unzip /tmp/bowtie.zip -d /usr/local && \
    ln -s /usr/local/bowtie-${BOWTIE_VERSION}/bowtie* /usr/local/bin && \
    rm /tmp/bowtie.zip

# Copy and install R dependencies
COPY deps.R /usr/local
RUN Rscript /usr/local/deps.R > /tmp/deps.log && rm /tmp/deps.log

# Install Psiplot
RUN curl -fsSL https://github.com/kcha/psiplot/archive/v${PSIPLOT_VERSION}.tar.gz -o /tmp/psiplot.tar.gz && \
    Rscript -e "install.packages('/tmp/psiplot.tar.gz', repos = NULL)" && \
    rm /tmp/psiplot.tar.gz

# Copy and install Vast-tools
COPY . /usr/local/vast-tools/
RUN chmod +x /usr/local/vast-tools/*.R && \
    ln -s /usr/local/vast-tools/vast-tools /usr/local/bin/vast-tools

# Ensure VASTDB directory exists
RUN mkdir -p /usr/local/vast-tools/VASTDB

# Run VAST-TOOLS setup script
WORKDIR /usr/local/vast-tools
RUN /usr/local/vast-tools/automatic_Hsa_Mmus_install.R --quiet

# Define shared volume and set default command
VOLUME /share
CMD ["bash"]
