options(rgl.useNULL=TRUE)

## load configuration
source("config.r")

## source auxiliary functions
source("server/hmmer.R")
source("server/db.R")
source("server/utils.R")

## load packages
library(bio3d)

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
library(shinyBS)

##library(dendextend)
##library(dynamicTreeCut)

if(configuration$db$use)
  library(RMySQL)

## system check
source("server/syscheck.R")
system_check()


shinyServer(function(input, output, session) {

  ## assign reactive values
  source("server/rv.R", local=TRUE)$value
  
  source("server/blast.R", local=TRUE)$value
  source("server/pdb.R", local=TRUE)$value
  source("server/seq.R", local=TRUE)$value
  source("server/pfam.R", local=TRUE)$value

  source("server/align.R", local=TRUE)$value
  source("server/fit.R", local=TRUE)$value
  source("server/pca.R", local=TRUE)$value
  source("server/nma.R", local=TRUE)$value
  source("server/pca-vs-nma.R", local=TRUE)$value


})
