#!/bin/bash

R CMD build ../bio3d
file=(`ls -t ../bio3d_*.tar.gz`)
R CMD check --as-cran $file

