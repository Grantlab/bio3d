#!/bin/bash


docker build . -t bio3dbuild -f util/Dockerfile
docker run -ti bio3dbuild

