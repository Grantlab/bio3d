options(rgl.useNULL=TRUE)

## load configuration
source("config.r")

library(bio3d)
library(lattice)
library(shiny)
library(rgl)
library(shinyRGL)

shinyServer(function(input, output, session) {

  source("server/pdb.R", local=TRUE)$value
  source("server/nma.R", local=TRUE)$value
  source("server/gl.R", local=TRUE)$value
  source("server/dccm.R", local=TRUE)$value
  source("server/overlap.R", local=TRUE)$value
  source("server/geostas.R", local=TRUE)$value
  
})
