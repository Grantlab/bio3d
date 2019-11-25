# docker build -f ../Dockerfile_html -t xinqyao/bio3d-html:latest .
# docker push xinqyao/bio3d-html:latest
rm -rf /tmp/bio3d-html
mkdir /tmp/bio3d-html

docker run -v /tmp/.X11-unix:/tmp/.X11-unix -v /tmp/.docker.xauth:/tmp/.docker.xauth -e XAUTHORITY=/tmp/.docker.xauth -e DISPLAY=:0 --rm -it --volume=/tmp/bio3d-html:/tmp/bio3d --workdir=/bio3d xinqyao/bio3d-html

mv /tmp/bio3d-html/* sandbox/
rm -rf /tmp/bio3d-html
