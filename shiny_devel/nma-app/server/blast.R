###########################
##-- BLAST (For overlap)  #
###########################

observeEvent(input$run_blast, {
  rv$blast <- run_blast2()
})

### Input sequence ###
get_sequence <- reactive({
  pdb <- final_pdb()
  seq <- as.fasta(pdbseq(pdb))
  return(seq)
})

### Annotate input PDB
input_pdb_annotation <- reactive({
  message("annotating")
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

### Run BLAST
run_blast2 <- reactive({
  message("blasting")
  input_sequence <- get_sequence()
  print(input_sequence)
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Blasting',
               detail = 'Please wait')
  progress$set(value = 3)

  ## use local version if possible
  if(configuration$hmmer$local)
    hmm <- phmmer_local(input_sequence)
  else
    hmm <- hmmer(input_sequence)

  if(!nrow(hmm) > 0)
    stop("No BLAST hits found")
  
  progress$set(value = 5)

  cutoff <- set_cutoff(hmm, cutoff=NULL)
  cutoff <- cutoff$cutoff

  hits <- hmm$score > cutoff
  acc <- hmm$acc[hits]

  rv$blast <- acc
  return(acc)

  ##out <- list(blast=hmm, pdbid=get_pdbid6())
  ##return(hmm)
})


## filters the blast results based on the cutoff
## returns logical vectors of all and limited hits
## as well as accession ids
#filter_hits <- reactive({
  ##blast <- run_blast()
#  blast <- rv$blast
#  cutoff <- set_cutoff(blast, cutoff=NULL)
#  cutoff <- cutoff$cutoff
#
#  hits <- blast$score > cutoff
#  acc <- blast$acc[hits]
#  return(acc)
#})


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
