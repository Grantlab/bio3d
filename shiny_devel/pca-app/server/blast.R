###########################
##-- BLAST and PDB INPUT  #
###########################

## set user data to store stuff
data_path <- reactive({
  dir <- paste0(format(Sys.time(), "%Y-%m-%d"), "_", randstr())
  path <- paste0(configuration$user_data, "/", dir)
  dir.create(path)
  return(path)
})


##- PDB input UI that is responsive to 'reset button' below
output$resetable_pdb_input <- renderUI({
  ## 'input$reset_pdb_input' is just used as a trigger for reset
  reset <- input$reset_pdb_input
  textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21") #)
})

### Input PDB object ### 
get_pdb <- reactive({
  anno <- input_pdb_annotation()
  if(is.vector(input$chainId)) {
    ind <- which(anno$chainId==input$chainId[1])
    anno <- anno[ind,]
  }
  else {
    anno <- anno[1,]
  }

  id <- anno$acc
  raw.files <- get.pdb(substr(id, 1,4), path=configuration$pdbdir$rawfiles, gzip=TRUE)
  pdb <- read.pdb(raw.files)
  pdb <- trim.pdb(pdb, chain=substr(id, 6,6))
  return(pdb)
})

### Input sequence ###
get_sequence <- reactive({

  ## option 1 - PDB code provided
  if(input$input_type == "pdb") {
    anno <- input_pdb_annotation()

    if(is.vector(input$chainId)) {
      ind <- which(anno$chainId==input$chainId[1])
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
      stop()

    inp <- unlist(strsplit(input$sequence, "\n"))
    inds <- grep("^>", inp, invert=TRUE)

    if(!length(inds)>0)
      stop("Error reading input sequence")

    inp <- toupper(paste(inp[inds], collapse=""))
    seq <- as.fasta(unlist(strsplit(inp, "")))
  }

  return(seq)
})

### Annotate input PDB
input_pdb_annotation <- reactive({
  if(is.null(input$pdbid)) {
    return()
  }

  if(nchar(input$pdbid)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())

    progress$set(message = 'Fetching PDB data',
                 detail = 'Please wait')
    progress$set(value = 3)

    anno <- get_annotation(input$pdbid, use_chain=FALSE)
    for(i in 4:5)
      progress$set(value = i)

    return(anno)
  }
  else {
    stop("Provide a 4 character PDB code")
  }
})

## short summary of the input PDB
output$input_pdb_summary <- renderPrint({
  input$input_type
  anno <- input_pdb_annotation()

  if(is.vector(input$chainId)) {
    ind <- which(anno$chainId==input$chainId[1])
    anno <- anno[ind,]
  }
  else {
    anno <- anno[1,]
  }

  cat(" Protein: ", anno$compound[1], "\n",
      "Species: ", anno$source[1], "\n")

})

### Run BLAST
run_blast <- reactive({
  if(is.null(input$chainId)) {
    stop()
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

    print(head(hmm))

    progress$set(value = 5)
    return(hmm)
  }
})


## filters the blast results based on the cutoff
## returns logical vectors of all and limited hits
## as well as accession ids
filter_hits <- reactive({
  input$input_type

  blast <- run_blast()
  cutoff <- set_cutoff(blast, cutoff=NULL)
  grps <- cutoff$grps
  cutoff <- cutoff$cutoff

  hits <- blast$score > cutoff
  acc <- blast$acc[hits]

  limit <- as.numeric(input$limit_hits)
  hits2 <- rep(FALSE, length(hits))
  hits2[1:limit] <- TRUE
  hits2 <- hits & hits2

  ## hits_all: logical vector of all hits
  ## hits: logical vector of limited hits
  ## acc: character vector of PDB ids
  out <- list(hits=hits2, hits_all=hits,
              acc=acc[hits2], acc_all=acc[hits], grps=grps)

  print(acc[hits2])

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

  cat("  * Possible cutoff values:   ", floor(gp.nums), "\n",
      "           Yielding Nhits:   ", gp.inds, "\n\n")

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
