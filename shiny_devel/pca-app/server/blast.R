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
rv$blast <- readRDS("2LUM_blast.RDS")
rv$pfam <- readRDS("2LUM_pfam.RDS")
rv$limit_hits <- 5
rv$cutoff <- 41
rv$sequence <- "MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE"
rv$pdb_codes <- "1TND, 1KJY"

observeEvent(input$pdbid, {
  if(nchar(input$pdbid)>3) {
    if(input$input_type == "pdb") {
      rv$pdbid <- substr(input$pdbid, 1, 4)
      
      if(rv$pdbid == "2LUM") {
        blast <- rv$blast
      }
      else {
        blast <- run_blast()
        rv$blast <- blast
        rv$pfam <- pdb.pfam(rv$pdbid)
      }
      
      cut <- set_cutoff(blast, cutoff=NULL)
      rv$cutoff <- cut$cutoff
      rv$limit_hits <- 5
    }
  }
})

observeEvent(input$sequence, {
  if(nchar(input$sequence)>10) {
    if(input$input_type == "sequence") {
      rv$sequence <- input$sequence
      blast <- run_blast()
      rv$blast <- blast
      cut <- set_cutoff(blast, cutoff=NULL)
      rv$cutoff <- cut$cutoff
      rv$limit_hits <- 5
    }
  }
})



observeEvent(input$chainId, {
  rv$chainid <- input$chainId
})

observeEvent(input$limit_hits, {
  rv$limit_hits <- as.numeric(input$limit_hits)
})

observeEvent(input$cutoff, {
  rv$cutoff <- as.numeric(input$cutoff)
})

observeEvent(input$reset_cutoff, {
  if(input$input_type != "multipdb") {
    blast <- rv$blast
    cut <- set_cutoff(blast, cutoff=NULL)

    updateSliderInput(session, "cutoff", value = cut$cutoff)
    updateSliderInput(session, "limit_hits", value = 5)
  }
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
      ##stop("Provide a PDB code of 4 characters")
      return()
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

  #if(is.null(input$input_type))
  #  return()

  ## option 1 - PDB code provided
  if(input$input_type == "pdb") {
    pdbid <- get_pdbid()
    chainid <- get_chainid()

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
    message(rv$sequence)

    if(nchar(rv$sequence)==0)
      return()
      ##stop("sequence is of length 0")

    inp <- unlist(strsplit(rv$sequence, "\n"))
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

  cat(anno$compound[1], "\n",
      "(", anno$source[1], ")")

})

## prints a longer log of the input PDB
output$pdb_log <- renderPrint({
  pdbid <- get_pdbid()
  chainid <- get_chainid()
  invisible(capture.output( pdb <- get_pdb() ))

  pdbsum(pdb, pdbid=pdbid, chainid=chainid)
    
})

##-- PFAM annotation of single or multiple PDBs
output$pfam_table <- DT::renderDataTable({
  message("input_pdb_pfam called")
  pdbid <- get_pdbid()
  #chainid <- get_chainid()
  if(is.null(pdbid)) {
    return()
  }

  ##pfam <- pdb.pfam(pdbid)
  pfam <- rv$pfam
  pfam['PFAM'] <- lapply(pfam['PFAM'], function(x) {
      pfid <- regmatches(x, gregexpr("(?<=\\().*?(?=\\))", x, perl=T))[[1]]
      link <- gsub(pfid,
                   tags$a(href=paste0("http://pfam.xfam.org/family/",pfid), target="_blank", pfid),
                   x)
  })
  colnames(pfam)=c("ID", "Name", "Description","eValue")
  DT::datatable(pfam, escape = FALSE, selection = "none",
    rownames = FALSE,
    options = list(
        dom = "t",
        autoWidth = TRUE,
        columnDefs = list(
            list( orderable = 'false', targets = c(0,1,2) )
        ),
    initComplete = JS(
    "function(settings, json) {",
    '$(this.api().table().header()).find("th").removeClass("sorting");',
    '$(this.api().table().header()).find("th").prop("onclick",null).off("click");',
    "}")
    )
  )
})
##To Do:
##       Rm '1 of 1' at bottom of table
##       Add PFAM URL link to Name. paste0("http://pfam.sanger.ac.uk/family/",acc)
##       Add this to multiple PDB IDs well/div also.


## Returns HMMER results
## dataframe with columns: acc, evalue, score, desc
run_blast <- reactive({
  message("run_blast called")

  if(input$input_type != "multipdb") {
    message("blasting")

    if(input$input_type == "pdb")
      if(is.null(rv$chainid))
        return(NULL)

    if(input$input_type == "sequence")
      if(!nchar(rv$sequence) > 10)
        return(NULL)

    input_sequence <- get_sequence()

    if(is.null(input_sequence))
      return(NULL)

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

    ## to uppercase_untouched
    hmm$acc <- format_pdbids(hmm$acc)
    ##blast <- hmm
    ##saveRDS(blast, file="2LUM_blast.RDS")

    progress$set(value = 5)
    return(hmm)
  }
  else {
    return()
  }
})


## filters the blast results based on the cutoff
## returns logical vectors of all and limited hits
## as well as accession ids
filter_hits <- reactive({
  message("filter_hits called")

  ##blast <- run_blast()
  blast <- rv$blast
  cutoff <- rv$cutoff

  ## logical vector
  inds <- blast$score > cutoff

  ## limited by input$limit_hits
  limit <- rv$limit_hits
  inds2 <- inds
  if(limit > 0 & limit < sum(inds)) {
    inds2[(limit+1):length(inds2)] <- FALSE
  }

  ## accession ids above cutoff
  hits <- blast$acc[inds]
  hits2 <- blast$acc[inds2]

  out <- list(hits=hits, inds=inds,
              hits_limited=hits2, inds_limited=inds2,
              cutoff=cutoff)

  return(out)
})


set_cutoff <- function(blast, cutoff=NULL) {

  if(is.null(blast))
    return(NULL)

  x <- blast
  cluster <- FALSE
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

  out <- list(inds=inds, gp.inds=gp.inds, grps=gps, cutoff=cutoff)
  return(out)
}
