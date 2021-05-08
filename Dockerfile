FROM r-base

MAINTAINER Lars Skjærven "larsss@gmail.com"

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libssl-dev \
    libcurl4-openssl-dev \
    libssh2-1-dev \
    libxml2-dev \
    netcdf-bin \ 
    libnetcdf-dev \ 
    wget \
    muscle \
    pymol

# Download and install dssp
RUN wget ftp://ftp.cmbi.umcn.nl/pub/software/dssp/dssp-2.0.4-linux-amd64 -O /usr/bin/dssp && \
    chmod +x /usr/bin/dssp

# Install required R packages
RUN R -e "install.packages('BH')"
RUN R -e "install.packages('Rcpp')"
RUN R -e "install.packages('RCurl')"
RUN R -e "install.packages('XML')"
RUN R -e "install.packages('ncdf4')"
RUN R -e "install.packages('igraph')"
RUN R -e "install.packages('bigmemory')"
RUN R -e "install.packages('knitr')"
RUN R -e "install.packages('testthat')"
RUN R -e "install.packages('httr')"

# Install Bio3D from source
COPY ./ver_devel/bio3d /tmp/bio3d
RUN R CMD INSTALL /tmp/bio3d && \
    rm -rf /tmp/bio3d

CMD ["R"]
