options(rgl.useNULL=TRUE)

## load configuration
source("config.r")

## source auxiliary functions
source("server/hmmer.R")
source("server/db.R")
source("server/utils.R")

## load packages
library(bio3d)
#library(devtools)
#load_all("~/workspace/bio3d/ver_devel/bio3d")


library(DT)
library(lattice)
library(shiny)
library(rCharts)
library(rgl)
library(shinyRGL)
library(reshape2)
library(maptools)
library(threejs)
library(abind)

if(configuration$db$use)
  library(RMySQL)


shinyServer(function(input, output, session) {

  source("server/blast.R", local=TRUE)$value
  source("server/blast_render.R", local=TRUE)$value
  source("server/align.R", local=TRUE)$value
  source("server/fit.R", local=TRUE)$value
  source("server/pca.R", local=TRUE)$value
  source("server/nma.R", local=TRUE)$value
  source("server/pca-vs-nma.R", local=TRUE)$value


})
