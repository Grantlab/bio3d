# docker build -f ../Dockerfile_html -t xinqyao/bio3d-html:latest .
# docker push xinqyao/bio3d-html:latest

if test ! -r sandbox; then
  mkdir sandbox
fi

rm -rf /tmp/bio3d-html
mkdir /tmp/bio3d-html

docker run -v /tmp/.X11-unix:/tmp/.X11-unix -v /tmp/.docker.xauth:/tmp/.docker.xauth -e XAUTHORITY=/tmp/.docker.xauth -e DISPLAY=:0 --rm -it --volume=/tmp/bio3d-html:/tmp/bio3d --workdir=/bio3d xinqyao/bio3d-html

fnm=$( basename /tmp/bio3d-html/* )
cp -r /tmp/bio3d-html/$fnm sandbox/
rm -f html
ln -s sandbox/$fnm/html html
rm -rf /tmp/bio3d-html
