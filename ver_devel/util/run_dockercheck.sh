#!/bin/bash

cran_repo="https://cran.uib.no"
cran="FALSE"
run_dont_test="FALSE"

# red hat linux
rhub_image=centos6-epel
rhub_image=centos6-epel-rdt
rhub_image=fedora-gcc-devel
rhub_image=fedora-clang-devel

# debian
rhub_image=ubuntu-gcc-devel
#rhub_image=ubuntu-gcc-release
#rhub_image=debian-gcc-devel
#rhub_image=debian-gcc-patched
#rhub_image=debian-gcc-release

dir=`pwd`
basename="$(basename $dir)"

if [ $basename != "util" ]; then
    echo "Please run script from bio3d/ver_devel/util";
    exit 1;
fi

basedir=`pwd`
dockerfile=Dockerfile

image=bio3d_$rhub_image
pkgs_rhel="RUN yum -y install openssl-devel libxml2-devel libcurl-devel"
pkgs_deb="RUN apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev netcdf-bin libnetcdf-dev"

devel_path='ENV PATH="/tmp/R-devel/bin:${PATH}"'

rm $dockerfile
touch $dockerfile

echo "FROM rhub/${rhub_image}" >> $dockerfile

if [[ $rhub_image == *"fedora"* || $rhub_image == *"centos"* ]]; then
  echo $pkgs_rhel >> $dockerfile
fi

if [[ $rhub_image == *"ubuntu"* || $rhub_image == *"debian"* ]]; then
  echo $pkgs_deb >> $dockerfile
fi

if [[ $rhub_image == *"devel"* ]]; then
  echo $devel_path >> $dockerfile
fi

devtools_command="devtools::check('/bio3d', cran=$cran, run_dont_test=$run_dont_test)"

cat <<EOF>>$dockerfile

# Install required R packages
RUN R -e "install.packages('bio3d', dep=T, repos = '$cran_repo')"
RUN R -e "install.packages('devtools', dep=T, repos = '$cran_repo')"

# Install Bio3D from source
COPY ./bio3d /bio3d
RUN R CMD INSTALL /bio3d

#CMD ["bash"]
CMD ["R", "-e", "$devtools_command"]
EOF

cd $basedir
cd ..
docker build . -t $image -f $basedir/$dockerfile
docker run -ti $image
cd $basedir
