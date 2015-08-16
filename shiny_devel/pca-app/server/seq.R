
## returns the sequence (fasta object)
get_sequence <- reactive({
  message("get_sequence called")

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

    if(!length(seq) > 10)
      warning(paste("sequence is of length", nchar(seq)))
    
    seq <- as.fasta(seq)
  }

  ## option 2 - sequence provided
  if(input$input_type == "sequence") {
    seq <- parse_inseq()
    inds <- check_inseq()$inds
    seq <- as.fasta(seq[inds])
  }

  return(seq)
})

## splits the input string, and returns
## a char vector 
parse_inseq <- reactive({
  inp <- unlist(strsplit(input$sequence, "\n"))
  inds <- grep("^>", inp, invert=TRUE)
  
  inp <- toupper(paste(inp[inds], collapse=""))
  seq <- unlist(strsplit(inp, ""))
  return(seq)
})

check_inseq <- reactive({
  seq <- parse_inseq()
  
  ## allowed values
  aa1 <- c("A","C","D","E","F","G",
           "H","I","K","L","M","N","P","Q",
           "R","S","T","V","W","Y")
  
  inds <- seq %in% aa1
  return(list(inds=inds, nallowed=unique(seq[!inds])))
  
})

output$input_seq_summary <- renderPrint({
  message("input_seq_summary called")

  blast <- rv$blast

  if(!is.null(blast)) {
    tophit <- format_pdbids(blast$acc[1])
    anno <- get_annotation(tophit)
  
    cat(paste("", tophit, "\n  ",
              anno$compound[1], "\n",
              "  (", anno$source[1], ")"))
  }
  else {
    cat("No hits found")
  }
})


output$input_seq_summary2 <- renderPrint({
  message("input_seq_summary2 called")

  ## blast hits
  blast <- rv$blast
  if(!is.null(blast))
    nhits <- length(blast$acc)
  else
    nhits <- 0

  cat(nhits)
})

output$input_seq_summary3 <- renderPrint({
  message("input_seq_summary3 called")

  ## trouble
  trbl <- check_inseq()

  if(sum(!trbl$inds)>0) {
    str <- paste("Non-protein residue codes detected:\n",
                 paste(trbl$nallowed, collapse=", "))
  }
  else {
    str <- "Seqence is OK"
  }

  if(sum(trbl$inds)<11)
    str <- "Sequence must be of length > 10"
    
  cat(str)
})

