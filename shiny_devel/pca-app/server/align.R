#####################
##-- Align
####################

fetch_pdbs <- reactive({
  
  if(length(input$pdb_ids) > 0) 
    ids <- input$pdb_ids
  else
    ids <- hits1()
  
  unq <- unique(substr(ids, 1,4))
  
  progress <- shiny::Progress$new(session, min=1, max=length(unq))
  
  
  progress$set(message = 'Fetching PDBs',
               detail = 'This may take some time...')
  
  ##raw.files <- get.pdb(ids, gzip=TRUE)
  raw.files <- vector("character", length(unq))
  for(i in 1:length(unq)) {
    raw.files[i] <- get.pdb(unq[i], path="raw_files", gzip=TRUE)
    progress$set(value = i)
  }
  progress$close()
  
  progress <- shiny::Progress$new(session, min=1, max=length(ids))
  
  
  progress$set(message = 'Splitting PDBs',
               detail = 'This may take some time...')
  
  ##files <- pdbsplit(raw.files, ids)
  
  ## this is possibly error prone
  files <- vector("character", length(ids))
  for(i in 1:length(unq)) {
    inds <- grep(unq[i], ids)
    
    files[inds] <- pdbsplit(raw.files[i], ids[inds], overwrite=FALSE)
    progress$set(value = i)
  }
  
  progress$close()
  return(files)
})


align <- reactive({
  files <- fetch_pdbs()
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Aligning PDBs',
               detail = 'This may take some time...')
  progress$set(value = 2)

  if(!input$reset_fasta)
    inc1 <<- 0
  
  reset <- (input$reset_fasta > inc1)
  inc1 <<- input$reset_fasta
  
  if(!is.null(input$fastafile_upload) & !reset) {
    infile <- input$fastafile_upload
    aln <- read.fasta(infile$datapath)
    pdbs <- read.fasta.pdb(aln)
  }
  else {
    pdbs <- pdbaln(files)
  }

  if(input$omit_missing) {
    conn <- inspect.connectivity(pdbs, cut=4.05)
    pdbs <- trim.pdbs(pdbs, row.inds=which(conn))
  }
  
  progress$set(value = 5)
  rownames(pdbs$ali) <- basename.pdb(rownames(pdbs$ali))
  pdbs <<- pdbs
  return(pdbs)
})


####   Alignment output   ####
output$alignment_summary <- renderPrint({
  invisible(capture.output( aln <- align() ))
  id <- aln$id
  ali <- aln$ali
  
  nstruct <- length(id)
  dims <- dim(ali)
  gaps <- gap.inspect(ali)
  dims.nongap <- dim(ali[, gaps$f.inds, drop = FALSE])
  dims.gap <- dim(ali[, gaps$t.inds, drop = FALSE])
  
  cat("Alignment dimensions:\n",
      "  ", dims[1L], " sequence rows \n",
      "  ", dims[2L], " position columns ", "(", dims.nongap[2L], " non-gap, ", dims.gap[2L], " gap) ", "\n", sep = "")
  cat("\n")
   
})

output$alignment <- renderUI({
  invisible(capture.output( aln <- align() ))

  ali <- aln$ali
  ids <- basename.pdb(aln$id)
  
  tmp1 <- conserv(ali, method="entropy10")
  tmp2 <- conserv(ali, method="identity")
  cons <- rep(" ", ncol(ali))
  cons[ tmp1==1 ] <- "^"
  cons[ tmp2==1 ] <- "*"
  
  width <- 80
  nseq <- length(ids)
  nres <- ncol(ali)
  
  block.start <- seq(1, nres, by=width)
  if(nres < width) {
    block.end <- nres
  } else {
    block.end <- unique(c( seq(width, nres, by=width), nres))
  }
  nblocks <- length(block.start)

  buff <- matrix("", ncol=3, nrow=nseq)
  block.annot  <- rep("", width)
  block.annot[ c(1,seq(10, width, by=10)) ] = "."

  out <- list()
  for(i in 1:nblocks) {
    positions <- block.start[i]:block.end[i]
    x <- ali[, positions, drop=FALSE]
    x <- cbind(buff, x, buff)
    
    annot <- block.annot[1:length(positions)]
    annot[length(annot)] <- block.end[i]
    annot[1] <- block.start[i]
    annot <- c(buff[1,], annot, buff[1,])
    
    annot2 <- cons[positions]
    annot2 <- c(buff[1,], annot2, buff[1,])

    spl <- unlist(strsplit(annot[4], ""))
    for(k in 1:length(spl))
      annot[4-k+1] <- spl[length(spl)-k+1]
    
    n <- length(annot)
    spl <- unlist(strsplit(annot[n-3], ""))
    for(k in 1:length(spl))
      annot[n-3+k-1] <- spl[k]


    ## annotation row: numbering
    tmp <- list()
    tmp[[1]] <- span(
      span(" ", class="aln_id"),
      lapply(annot[1:3],     function(x) span(x, class="aln_buff")),
      lapply(annot[4:(n-4)], function(x) span(x, class="aln_aminoacid")),
      lapply(annot[(n-3):n], function(x) span(x, class="aln_buff")),
      class="aln_row"
      )

    ## sequence row
    for(j in 1:nrow(x)) {
      tmp[[j+1]] <- span(
        span(ids[j], class="aln_id"),
        lapply(x[j, 1:3],     function(x) span(x, class="aln_buff", class=x)),
        lapply(x[j, 4:(n-4)], function(x) span(x, class="aln_aminoacid", class=x)),
        lapply(x[j, (n-3):n], function(x) span(x, class="aln_buff", class=x)),
        class="aln_row"
        )
    }

    ## annotation row: conservation
    tmp[[j+2]] <- span(
      span(" ", class="aln_id"),
      lapply(annot2[1:3],     function(x) span(x, class="aln_buff")),
      lapply(annot2[4:(n-4)], function(x) span(x, class="aln_aminoacid")),
      lapply(annot2[(n-3):n], function(x) span(x, class="aln_buff")),
      class="aln_row"
      )

    out[[i]] <- span(tmp, class="aln_block")
  }
    
  pre(class="alignment", 
      out
      )
    
})



####    Download functions ####
output$pdbsRData = downloadHandler(
  filename = 'pdbs.RData',
  content = function(file) {
    save(pdbs, file=file)
  })


output$fastafile = downloadHandler(
  filename = 'aln.fasta',
  content = function(file) {
    write.fasta(align(), file=file)
  })
