
###########################
##-- Sets reactive values for the app
###########################

### set default vaules
rv <- reactiveValues()
rv$input_type <- "pdb"

## input = pdb
rv$pdbid <- "2LUM"
rv$chainid <- "A"
rv$blast <- readRDS("2LUM_blast.RDS")
rv$limit_hits <- 5
rv$cutoff <- 41

## input = sequence
rv$sequence <- "MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE"

## input = multi pdb
rv$pdb_codes <- "1TND, 1KJY_A"
rv$pdbids <- c("1TND_A", "1TND_B", "1TND_C", "1KJY_A")

## for subsequent tabs
rv$aligned <- FALSE
rv$fitted <- FALSE
rv$modes <- NULL

## used for hiding stuff on the blast tab 
rv$hits_found <- TRUE

## selected accession ids
rv$selacc <- NULL

## pdbs object for nma calculations
rv$pdbs4nma <- NULL
