FROM rhub/fedora-clang-devel 

RUN yum -y install openssl-devel libxml2-devel libcurl-devel

# Install required R packages
RUN /tmp/R-devel/bin/R -e "install.packages('bio3d', dep=T, repos = 'https://cran.uib.no')"
RUN /tmp/R-devel/bin/R -e "install.packages('devtools', dep=T, repos = 'https://cran.uib.no')"


# Install Bio3D from source
COPY ./bio3d /tmp/bio3d
COPY ./util/check.R /tmp/check.R
RUN /tmp/R-devel/bin/R -e "devtools::install('/tmp/bio3d')"

#CMD ["/tmp/R-devel/bin/R"]
CMD ["/tmp/R-devel/bin/Rscript", "/tmp/check.R"]
