#!/bin/bash

R CMD build ../bio3d
R CMD check --as-cran bio3d_2.0.tar.gz

