#####################
##-- Align
####################

fetch_pdbs <- reactive({
  
  print(input$pdb_ids)
  if(length(input$pdb_ids) > 0) 
    ids <- input$pdb_ids
  else
    ids <- get_hit_ids()
  
  unq <- unique(substr(ids, 1,4))

  print(ids)
  print(unq)
  
  progress <- shiny::Progress$new(session, min=1, max=length(unq))
  
  
  progress$set(message = 'Fetching PDBs',
               detail = 'This may take some time...')
  
  ##raw.files <- get.pdb(ids, gzip=TRUE)
  raw.files <- vector("character", length(unq))
  for(i in 1:length(unq)) {
    raw.files[i] <- get.pdb(unq[i], path=configuration$pdbdir$rawfiles, gzip=TRUE)
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
    
    files[inds] <- pdbsplit(raw.files[i], ids[inds], overwrite=FALSE,
                            path=configuration$pdbdir$splitfiles)
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
    pdbs <- pdbaln(files, verbose=TRUE,
                   exefile=configuration$muscle$exefile,
                   outfile=tempfile(pattern="aln", fileext=".fasta"))
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
  
  nseq <- length(ids)
  nres <- ncol(ali)
  width <- ifelse(nres > 80, 80, nres)
  
  block.start <- seq(1, nres, by=width)
  if(nres < width) {
    block.end <- nres
  } else {
    block.end <- unique(c( seq(width, nres, by=width), nres))
  }
  nblocks <- length(block.start)
  bufsize <- nchar(nres)-1
  

  x <- matrix("", ncol=bufsize + width + bufsize, nrow=nseq)
  block.annot <- x[1,]
  block.annot[c(bufsize + 1, seq(10+bufsize, width+bufsize, by=10)) ] <- "."
  cons.annot <- x[1,]

  #buff <- matrix("", ncol=bufsize, nrow=nseq)
  #block.annot  <- rep(" ", width)
  #block.annot[ c(1,seq(10, width, by=10)) ] = "."

  out <- list()
  for(i in 1:nblocks) {

    positions <- block.start[i]:block.end[i]
    n <- length(positions)
    aln.inds <- (bufsize + 1):(n + bufsize)
    buf.inds1 <- 1:bufsize
    buf.inds2 <- (1 + bufsize + n):(2 * bufsize + n)

    if(n < (ncol(x) - 2 * bufsize)) {
      x <- x[, 1:(n + 2 * bufsize)]
      block.annot  <- block.annot[1:(n + 2 * bufsize)]
      cons.annot  <- block.annot[1:(n + 2 * bufsize)]
    }
    x[, aln.inds] <- ali[, positions, drop=FALSE]
   
    annot <-  block.annot
    annot[aln.inds[1]] <- block.start[i]
    annot[aln.inds[length(aln.inds)]] <- block.end[i]
    
    annot2 <- cons.annot
    annot2[aln.inds] <- cons[positions]

    m <- buf.inds1[bufsize]+1
    spl <- unlist(strsplit(annot[m], ""))
    for(k in 1:length(spl))
      annot[m-k+1] <- spl[length(spl)-k+1]   
    
    m <- buf.inds2[1]-1
    spl <- unlist(strsplit(annot[m], ""))
    for(k in 1:length(spl))
      annot[m+k-1] <- spl[k]


    ## annotation row: numbering
    tmp <- list()
    tmp[[1]] <- span(
      span(" ", class="aln_id"),
      lapply(annot[buf.inds1],     function(x) span(x, class="aln_buff")),
      lapply(annot[aln.inds], function(x) span(x, class="aln_aminoacid")),
      lapply(annot[buf.inds2], function(x) span(x, class="aln_buff")),
      class="aln_row"
      )

    ## sequence row
    for(j in 1:nrow(x)) {
      tmp[[j+1]] <- span(
        span(ids[j], class="aln_id"),
        lapply(x[j, buf.inds1],     function(x) span(x, class="aln_buff", class=x)),
        lapply(x[j, aln.inds], function(x) span(x, class="aln_aminoacid", class=x)),
        lapply(x[j, buf.inds2], function(x) span(x, class="aln_buff", class=x)),
        class="aln_row"
        )
    }

    ## annotation row: conservation
    tmp[[j+2]] <- span(
      span(" ", class="aln_id"),
      lapply(annot2[buf.inds1],     function(x) span(x, class="aln_buff")),
      lapply(annot2[aln.inds], function(x) span(x, class="aln_aminoacid")),
      lapply(annot2[buf.inds2], function(x) span(x, class="aln_buff")),
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
