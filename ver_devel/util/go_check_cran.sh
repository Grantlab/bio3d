#!/bin/bash

R CMD build bio3d
R CMD check --as-cran bio3d_1.2.tar.gz

