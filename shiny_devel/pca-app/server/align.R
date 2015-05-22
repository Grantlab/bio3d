#####################
##-- Align
####################

get_acc <- reactive({

  ## check first input from checkboxes
  if(!is.null(input$pdb_ids)) {
    ids <- input$pdb_ids

    if(length(ids) > 0)
      return(toupper(ids))
  }

  if(input$input_type != "multipdb") {
    blast <- run_blast()
    hits <- filter_hits()
    acc <- hits$acc
    return(toupper(acc))
  }
  else {
    acc <- toupper(unique(trim(unlist(strsplit(input$pdb_codes, ",")))))
    acc <- acc[acc!=""]

    anno <- get_annotation(acc, use_chain=FALSE)
    inds <- unlist(sapply(acc, grep, anno$acc))
    anno <- anno[inds, ]
    acc <- anno$acc
    return(toupper(acc))
   }
})

fetch_pdbs <- reactive({
  ids <- get_acc()
  unq <- unique(substr(ids, 1,4))

  progress <- shiny::Progress$new(session, min=1, max=length(unq))
  progress$set(message = 'Fetching PDBs',
               detail = 'Please wait ...')

  tryfiles <- paste0(configuration$pdbdir$rawfiles, "/", unq, ".pdb")
  if(all(file.exists(tryfiles))) {
    raw.files <- tryfiles
  }
  else {
    raw.files <- vector("character", length(unq))
    for(i in 1:length(unq)) {
      raw.files[i] <- get.pdb(unq[i], path=configuration$pdbdir$rawfiles, gzip=TRUE)
      progress$set(value = i)
    }
  }
  gc()
  progress$close()

  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = 'Splitting PDBs',
               detail = 'Please wait ...',
               value = 0)

  tryfiles <- paste0(configuration$pdbdir$splitfiles, "/", ids, ".pdb")
  if(all(file.exists(tryfiles))) {
    files <- tryfiles
  }
  else {
    files <- pdbsplit(pdb.files=raw.files, ids=ids, overwrite=FALSE,
                      path=configuration$pdbdir$splitfiles,
                      progress=progress)
  }
  gc()
  progress$close()
  return(files)
})


align <- reactive({
  files <- fetch_pdbs()

  progress <- shiny::Progress$new()
  on.exit(progress$close())

  progress$set(message = 'Aligning PDBs',
               detail = 'Please wait ...',
               value = 0)

  if(!input$reset_fasta)
    inc1 <<- 0

  reset <- (input$reset_fasta > inc1)
  inc1 <<- input$reset_fasta

  if(!is.null(input$fastafile_upload) & !reset) {
    infile <- input$fastafile_upload
    aln <- read.fasta(infile$datapath)
    pdbs <- read.fasta.pdb(aln, progress=progress)
  }
  else {
    pdbs <- pdbaln(files, verbose=TRUE,
                   exefile=configuration$muscle$exefile,
                   outfile=tempfile(pattern="aln", fileext=".fasta"),
                   progress=progress)
  }

  if(input$omit_missing) {
    conn <- inspect.connectivity(pdbs, cut=4.05)
    pdbs <- trim.pdbs(pdbs, row.inds=which(conn))
  }

  rownames(pdbs$ali) <- basename.pdb(rownames(pdbs$ali))
  progress$close()
  gc()
  ##save(pdbs, file="pdbs.RData")
  return(pdbs)
})

####################################
####   Alignment output         ####
####################################

check_aln <- reactive({
  aln <- align()
  gaps <- gap.inspect(aln$ali)

  ## hard-coded limit here -- matches stop.at argument in core.find()
  if(!length(gaps$f.inds) > 15) {
    stop("Insufficent non-gap regions in alignment to proceed")
  }
})

output$alignment_summary <- renderPrint({
  invisible(capture.output( aln <- align() ))
  check_aln()

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

  progress <- shiny::Progress$new(session, min=1, max=nblocks)
  progress$set(message = 'Generating HTML',
               detail = 'Please wait')


  out <- list()
  for(i in 1:nblocks) {
    x[x!=""] <- ""
    positions <- block.start[i]:block.end[i]
    n <- length(positions)
    aln.inds <- (bufsize + 1):(n + bufsize)
    buf.inds1 <- 1:bufsize
    buf.inds2 <- (1 + bufsize + n):(2 * bufsize + n)

    if(n < (ncol(x) - 2 * bufsize)) {
      x <- x[, 1:(n + 2 * bufsize), drop=FALSE]
      block.annot  <- block.annot[1:(n + 2 * bufsize)]
      cons.annot  <- block.annot[1:(n + 2 * bufsize)]
    }
    #x[, aln.inds] <- ali[, positions, drop=FALSE]
    require(abind)
    x.mat <- abind(x, x, along=3)
    x.mat[, aln.inds, 1] <- ali[, positions, drop=FALSE]
    x.mat[, aln.inds, 2] <- aln$resno[, positions, drop=FALSE]

    annot <-  block.annot
    annot[aln.inds[1]] <- block.start[i]
    annot[aln.inds[length(aln.inds)]] <- block.end[i]

    annot2 <- cons.annot
    annot2[aln.inds] <- cons[positions]
    annot2[buf.inds1] <- ""
    annot2[buf.inds2] <- ""

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
    print(aln)
    ## sequence row
    for(j in 1:nrow(x)) {
      tmp[[j+1]] <- span(
        span(ids[j], class="aln_id"),
        lapply(x[j, buf.inds1], function(x) span(x, class="aln_buff", class=x)),
        lapply(seq_along(1:length(x.mat[j, aln.inds,1])), function(x) span(x.mat[j,aln.inds,1][x], class="aln_aminoacid", class=x.mat[j,aln.inds,1][x]
            , "title"="", "data-placement"="bottom",
            "data-original-title"="Residue info","data-toggle"="popover",
            "data-content"=paste0("<table class=\"tb_pop\"><tbody><tr><td>PDB residue number:</td><td>",
               x.mat[j,aln.inds,2][x],"</td></tr><tr><td>Other</td><td>info</td></tr></tbody></table>"
            )
                           )),
        lapply(x[j, buf.inds2], function(x) span(x, class="aln_buff", class=x)),
        class="aln_row"
        )
    }

    ## annotation row: conservation
    tmp[[j+2]] <- span(
      span(" ", class="aln_id"),
      lapply(annot2[buf.inds1],     function(x) span(x, class="aln_buff")),
      lapply(annot2[aln.inds], function(x) span(x, class="aln_aminoacid", title='boo')),
      lapply(annot2[buf.inds2], function(x) span(x, class="aln_buff")),
      class="aln_row"
      )

    out[[i]] <- span(tmp, class="aln_block")
    progress$set(value = i)
  }
  progress$close()
  gc()
  #cat(as.character(out))
  pre(class="alignment",
      out
      )
})



####################################
####     Download functions     ####
####################################

aln2file <- reactive({
  fn <- paste0(data_path(), '/aln.fasta')
  write.fasta(align(), file=fn)
  return(fn)
})

output$fastafile = downloadHandler(
  filename = 'aln.fasta',
  content = function(file) {
    aln2file()
  })
