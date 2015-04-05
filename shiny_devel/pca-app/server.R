library(bio3d)
library(lattice)
library(shiny)
library(rCharts)
library(RMySQL)
source("server/hmmer.R")
source("server/db.R")

## wget http://www.uniprot.org/docs/pdbtosp.txt
## wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz


## Define server logic for PCA shiny demo

shinyServer(function(input, output, session) {

  source("server/input.R", local=TRUE)$value
  source("server/blast.R", local=TRUE)$value
  source("server/align.R", local=TRUE)$value
  source("server/fit.R", local=TRUE)$value

  source("server/pca.R", local=TRUE)$value
  source("server/nma.R", local=TRUE)$value
  
  
})
