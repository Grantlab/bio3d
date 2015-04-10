## load configuration
source("config.r")

## source auxiliary functions
source("server/hmmer.R")
source("server/db.R")

## load packages
library(bio3d)
library(lattice)
library(shiny)
library(rCharts)

if(configuration$db$use)
  library(RMySQL)


shinyServer(function(input, output, session) {

  source("server/input.R", local=TRUE)$value
  source("server/blast.R", local=TRUE)$value
  source("server/align.R", local=TRUE)$value
  source("server/fit.R", local=TRUE)$value

  source("server/pca.R", local=TRUE)$value
  source("server/nma.R", local=TRUE)$value
  
  
})
