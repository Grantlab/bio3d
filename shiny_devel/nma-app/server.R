options(rgl.useNULL=TRUE)

## load configuration
source("config.r")

source("server/hmmer.R")
source("server/utils.R")
source("server/db.R")
source("server/utils.R")


library(bio3d)
library(lattice)
library(shiny)
library(rgl)
library(shinyRGL)
library(DT)
library(threejs)


if(configuration$db$use)
  library(RMySQL)

## system check
source("server/syscheck.R")
system_check()


shinyServer(function(input, output, session) {

  source("server/pdb.R", local=TRUE)$value
  source("server/nma.R", local=TRUE)$value
  source("server/gl.R", local=TRUE)$value
  source("server/dccm.R", local=TRUE)$value
  source("server/overlap.R", local=TRUE)$value
  source("server/geostas.R", local=TRUE)$value

  source("server/blast.R", local=TRUE)$value
  source("server/blast_render.R", local=TRUE)$value
  
})
