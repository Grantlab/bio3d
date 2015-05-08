## load configuration
source("config.r")

## source auxiliary functions
source("server/hmmer.R")
source("server/db.R")
source("server/utils.R")

## load packages
library(bio3d)
library(lattice)
library(shiny)
library(rCharts)
library(reshape2)
library(maptools)

if(configuration$db$use)
  library(RMySQL)


shinyServer(function(input, output, session) {

  source("server/blast.R", local=TRUE)$value
  source("server/blast_render.R", local=TRUE)$value
  source("server/align.R", local=TRUE)$value
  source("server/fit.R", local=TRUE)$value

  source("server/pca.R", local=TRUE)$value
  source("server/nma.R", local=TRUE)$value
  
})
