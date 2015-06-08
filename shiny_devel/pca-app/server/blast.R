###################################
##-- User path for storing stuff  #
################################### 
data_path <- reactive({
  dir <- paste0(format(Sys.time(), "%Y-%m-%d"), "_", randstr())
  path <- paste0(configuration$user_data, "/", dir)
  dir.create(path)
  return(path)
})


###########################
##-- PDB AND BLAST INPUT  #
###########################

### set default vaules
rv <- reactiveValues()
rv$pdbid <- "2LUM"
rv$chainid <- "A"
rv$limit_hits <- 5
rv$cutoff <- 41

observeEvent(input$pdbid, {
  rv$pdbid <- input$pdbid
})

observeEvent(input$chainId, {
  rv$chainid <- input$chainId
})

observeEvent(input$limit_hits, {
  rv$limit_hits <- input$limit_hits
})

observeEvent(input$cutoff, {
  rv$cutoff <- input$cutoff
})

observeEvent(input$reset_pdbid, {
  updateTextInput(session, "pdbid", value = "2LUM")
})


## returns PDB code (4-characters)
get_pdbid <- reactive({
  if (is.null(rv$pdbid))
    return()
  else {
    pdbid <- rv$pdbid
    
    if(!nchar(pdbid)==4) {
      stop("Provide a PDB code of 4 characters")
    }
    
    return(pdbid)
  }
})

## returns the selected chain ID
get_chainid <- reactive({
  if (is.null(rv$chainid))
    return()
  else
    return(rv$chainid)
})

## returns all chain IDs in PDB
get_chainids <- reactive({
  pdbid <- get_pdbid()
  
  if(is.null(pdbid))
    return()
  
  anno <- input_pdb_annotation()
  return(anno$chainId)
})

## returns the PDB object (trimmed to chain ID)
get_pdb <- reactive({
  message("get_pdb called")
  pdbid <- get_pdbid()
  chainid <- get_chainid()
  
  anno <- input_pdb_annotation()

  if(is.vector(input$chainId)) {
    ind <- which(anno$chainId==input$chainId[1])
    anno <- anno[ind,]
  }
  else {
    anno <- anno[1,]
  }

  id <- anno$acc
  if(configuration$pdbdir$archive) {
    ids <- paste0(tolower(substr(id, 1, 4)), "_", substr(id, 6, 6))
    raw.files <- paste0(configuration$pdbdir$splitfiles, "/", substr(ids, 2, 3), "/pdb", ids, ".ent.gz")
    
    if(!file.exists(raw.files))
      stop("PDB not found")
    
    pdb <- read.pdb(raw.files)
  }
  else {
    raw.files <- get.pdb(unique(substr(id, 1,4)), path=configuration$pdbdir$rawfiles, gzip=TRUE)
    pdb <- read.pdb(raw.files)
    pdb <- trim.pdb(pdb, chain=substr(id, 6,6))
  }
  
  return(pdb)
})

## returns the sequence (fasta object)
get_sequence <- reactive({
  message("get_sequence called")
  pdbid <- get_pdbid()
  chainid <- get_chainid()
  
  if(is.null(input$input_type))
    stop("no input provided")

  ## option 1 - PDB code provided
  if(input$input_type == "pdb") {
    anno <- input_pdb_annotation()

    if(is.vector(input$chainId)) {
      ind <- which(anno$chainId==chainid[1])
      seq <- unlist(strsplit(anno$sequence[ind], ""))
    }
    else {
      seq <- unlist(strsplit(anno$sequence[1], ""))
    }
    seq <- as.fasta(seq)
  }

  ## option 2 - sequence provided
  if(input$input_type == "sequence") {
    if(nchar(input$sequence)==0)
      stop("sequence is of length 0")

    inp <- unlist(strsplit(input$sequence, "\n"))
    inds <- grep("^>", inp, invert=TRUE)

    if(!length(inds)>0)
      stop("Error reading input sequence")

    inp <- toupper(paste(inp[inds], collapse=""))
    seq <- as.fasta(unlist(strsplit(inp, "")))
  }

  return(seq)
})

## returns annotation data for input PDB
input_pdb_annotation <- reactive({
  message("input_pdb_annotation called")
  
  pdbid <- get_pdbid()
  
  if(nchar(pdbid)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())

    progress$set(message = 'Fetching PDB data',
                 detail = 'Please wait')
    progress$set(value = 3)

    anno <- get_annotation(pdbid, use_chain=FALSE)
    progress$set(value = 5)

    return(anno)
  }
  else {
    stop("Provide a 4 character PDB code")
  }
})

## prints a short summary of the input PDB
output$input_pdb_summary <- renderPrint({
  message("input_pdb_summary called")
  pdbid <- get_pdbid()
  chainid <- get_chainid()
  if(is.null(pdbid)) {
    return()
  }
  
  input$input_type
  anno <- input_pdb_annotation()
  
  if(is.vector(chainid)) {
    ind <- which(anno$chainId==chainid[1])
    anno <- anno[ind,]
  }
  else {
    anno <- anno[1,]
  }

  cat(" Protein: ", anno$compound[1], "\n",
      "Species: ", anno$source[1], "\n")

})

## returns BLAST results
run_blast <- reactive({
  message("run_blast called")
  
  if(is.null(rv$chainid)) {
    stop("no chainId provided")
  }

  if(input$input_type == "multipdb") {
    stop(" no need to blast when input is multipdb")
  }
  else {
    message("blasting")
    input_sequence <- get_sequence()

    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())

    progress$set(message = 'Blasting',
                 detail = 'Please wait')
    progress$set(value = 2)

    ## use local version if possible
    if(configuration$hmmer$local)
      hmm <- phmmer_local(input_sequence)
    else
      hmm <- hmmer(input_sequence)

    if(!nrow(hmm) > 0)
      stop("No BLAST hits found")

    progress$set(value = 5)
    return(hmm)
  }
})


## filters the blast results based on the cutoff
## returns logical vectors of all and limited hits
## as well as accession ids
filter_hits <- reactive({
  message("filter_hits called")

  blast <- run_blast()
  cutoff <- set_cutoff(blast, cutoff=NULL)
  grps <- cutoff$grps
  cutoff <- cutoff$cutoff

  hits <- blast$score > cutoff
  acc <- blast$acc[hits]
  hits2 <- rep(FALSE, length(hits))
  
  if(!is.null(rv$limit_hits)) {
    limit <- as.numeric(rv$limit_hits)
    if(limit>0)
      hits2[1:limit] <- TRUE
  }
  hits2 <- hits & hits2

  ## hits_all: logical vector of all hits
  ## hits: logical vector of limited hits
  ## acc: character vector of PDB ids
  out <- list(hits=hits2, hits_all=hits,
              acc=acc[hits2], acc_all=acc[hits], grps=grps)

  return(out)
})


set_cutoff <- function(blast, cutoff=NULL) {

  x <- blast
  cluster <- TRUE
  cut.seed <- NULL

  if(is.null(x$evalue))
    stop("missing evalues")

  x$mlog.evalue <- x$score

  ##- Find the point pair with largest diff evalue
  dx <- abs(diff(x$mlog.evalue))
  dx.cut = which.max(dx)


  if(!is.null(cutoff)) {
    ##- Use suplied cutoff
    gps = rep(2, length(x$mlog.evalue))
    gps[ (x$mlog.evalue >= cutoff) ] = 1

  }
  else {
    if(cluster) {
      ## avoid clustering with many hits
      nhit <- length(x$mlog.evalue)
      if(nhit > 2000) {
        cluster <- FALSE
      }
    }

    if(is.null(cut.seed)) {
      ## Use mid-point of largest diff pair as seed for
      ##  cluster grps (typical PDB values are ~110)
      cut.seed = mean( x$mlog.evalue[dx.cut:(dx.cut+1)] )
    }

    if(cluster){
      ##- Partition into groups via clustering
      ##  In future could use changepoint::cpt.var
      hc <- hclust( dist(x$mlog.evalue) )
      if(!is.null(cutoff)) { cut.seed=cutoff }
      gps <- cutree(hc, h=cut.seed)
    }

    if(!cluster || (length(unique(gps))==1)) {
      ##- Either we don't want to run hclust or hclust/cutree
      ##   has returned only one grp so here we will divide
      ##   into two grps at point of largest diff
      gps = rep(2, length(x$mlog.evalue))
      gps[1:dx.cut]=1
    }
  }

  gp.inds <- na.omit(rle2(gps)$inds)
  gp.nums <- x$mlog.evalue[gp.inds]

  #cat("  * Possible cutoff values:   ", floor(gp.nums), "\n",
  #    "           Yielding Nhits:   ", gp.inds, "\n\n")

  if( is.null(cutoff) ) {
    ## Pick a cutoff close to cut.seed
    i <- which.min(abs(gp.nums - cut.seed))
    cutoff <- floor( gp.nums[ i ] )
  }

  inds <- x$mlog.evalue >= cutoff
  cat("  * Chosen cutoff value of:   ", cutoff, "\n",
      "           Yielding Nhits:   ", sum(inds), "\n")

  out <- list(inds=which(inds), gp.inds=gp.inds, grps=gps, cutoff=cutoff)
  return(out)

}
