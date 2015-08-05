#####################
##-- Align
####################

## selected accession ids
rv$selacc <- NULL

observeEvent(input$selected_pdbids, {

  if(!length(input$selected_pdbids) > 0) {
    rv$selacc <- get_acc()
  }
  
  if(!all(rv$selacc %in% input$selected_pdbids))
    rv$selacc <- input$selected_pdbids
  
  if(!all(input$selected_pdbids %in% rv$selacc))
    rv$selacc <- input$selected_pdbids
  
})


observeEvent(input$omit_missing, {

  allacc <- get_acc()

  if(input$omit_missing) {
    pdbs <- align_pdbs()
    conn <- inspect.connectivity(pdbs, cut=4.05)

    if(!all(conn)) {

      selacc <- rv$selacc[ rv$selacc %in% pdbs$lab[conn] ]
      rv$selacc <- selacc
      
      updateSelectInput(session, "selected_pdbids", 
                        choices = allacc, selected = selacc)
    }
  }
  else {
    rv$selacc <- allacc
    updateSelectInput(session, "selected_pdbids", 
                      choices = allacc, selected = allacc)
      
  }
  
})
  

observeEvent(input$reset_fasta, {

  rv$selacc <- get_acc()
  acc <- rv$selacc
  updateSelectInput(session, "selected_pdbids", 
                    choices = acc, selected = acc)
  
  updateCheckboxInput(session, "omit_missing", value = FALSE)
})

output$include_hits <- renderUI({

  acc <- get_acc()
  names(acc) <- acc
  
  selectInput("selected_pdbids", "Include / exclude hits",
              choices = acc, selected = acc, multiple=TRUE)

})

## returns the accession IDs of selected hits (from the SEARCH tab)
get_acc <- reactive({

  row.inds <- sort(as.numeric(input$blast_table_rows_selected))
  
  message( paste(row.inds, collapse=" ") ) 
  
  if(input$input_type != "multipdb") {
    blast <- rv$blast
    acc <- blast$acc[ row.inds ]
  }
  else {
    acc <- get_multipdbids()
  }
  
  acc <- format_pdbids(acc)
  rv$selacc <- acc
  
  return(acc)
})



## returns the filenames of splitted PDBs
split_pdbs <- reactive({
  ids <- get_acc()
  unq <- unique(substr(ids, 1,4))

  ## archive mode - pre splitted PDBs
  if(configuration$pdbdir$archive) {
    
    ## file names are in lower case except from chain ID
    ## e.g. pdb2lum_A.ent.gz
    ids <- format_pdbids(ids, casefmt=tolower)
    
    files <- paste0(configuration$pdbdir$splitfiles, "/",
                    substr(ids, 2, 3), "/pdb", ids, ".ent.gz")

    if(any(!file.exists(files))) {
      ## pdb2lum_A.ent.gz --> 2lum_A
      missing_ids <- pdbfilename2label(
        files[ !file.exists(files) ]
        )
      missing <- paste(missing_ids, collapse=", ")
      ##stop(paste("PDB file(s) missing:", missing))
      warning(paste("PDB file(s) missing:", missing))

      files <- files[ file.exists(files) ]

      if(!length(files) > 0)
        stop("PDB files not found")
    }
  }

  ## download mode - downloads and splits PDBs
  else {
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
  }

  return(files)
})

## performs the actual alignment of the PDBs
align_pdbs <- reactive({
  files <- split_pdbs()

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

  
  if(configuration$pdbdir$archive) {
    ## pdb2lum_A.ent.gz --> 2lum_A
    pdbs$lab <- pdbfilename2label(pdbs$id)
  }
  else {
    pdbs$lab <- basename.pdb(pdbs$id)
  }
  pdbs$lab <- format_pdbids(pdbs$lab)
  
 rownames(pdbs$ali) <- pdbs$lab
 progress$close()
 
 gc()
 return(pdbs)
})

## filters the existing alignment obtained by align_pdbs()
## misleading name, but provides the final alignment 
align <- reactive({
  pdbs <- align_pdbs()

  #if(input$omit_missing) {
  #  conn <- inspect.connectivity(pdbs, cut=4.05)
  #  if(!all(conn)) {
  #    labs <- pdbs$lab
  #    pdbs <- trim(pdbs, row.inds=which(conn))
  #    pdbs$lab <- labs[conn]
  #  }
  #}
  
  if(length(rv$selacc) > 0) {
    if(length(rv$selacc) != length(pdbs$lab)) {
      ids <- rv$selacc
      labs <- pdbs$lab
      
      inds <- which(pdbs$lab %in% ids)
      if(length(inds) > 0) {
       pdbs <- trim(pdbs, row.inds = inds)
        pdbs$lab <- labs[inds]
      }
    }
  }

  gc()
  return(pdbs)
})

seqide <- reactive({
  pdbs <- align()
  ide <- seqidentity(pdbs)
  rownames(ide) <- pdbs$lab
  colnames(ide) <- pdbs$lab
  return(ide)
})

seqconserv <- reactive({
  pdbs <- align()
  return(conserv(pdbs$ali, method=input$conserv_method))
})

####################################
####   Alignment output         ####
####################################

check_aln <- reactive({
  aln <- align()
  gaps <- gap.inspect(aln$ali)

  ## hard-coded limit here -- matches stop.at argument in core.find()
  if(!length(gaps$f.inds) > 5) {
    stop("Insufficient non-gap regions in alignment to proceed")
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


  wanted <- get_acc()
  
  missing <- !wanted %in% aln$lab
  if(any(missing))
    cat("Omitted PDBs:", paste(sort(wanted[missing]), collapse=", "))
        
  
})


output$missres_summary <- renderPrint({
  invisible(capture.output( pdbs <- align() ))
  check_aln()

  ids <- pdbs$lab
  nstructs <- length(ids)

  conn <- inspect.connectivity(pdbs, cut=4.05)
  if(sum(!conn) > 0) {
    inds <- which(!conn)

    cat(sum(!conn), "PDBs with missing in-structure residues:\n")
    cat("  ", paste(ids[!conn], collapse=", "))
    cat("\n")
  }
  else {
    cat("No PDBs with missing in-structure residues\n")
  }
})




output$alignment <- renderUI({
  message("generating HTML")
  invisible(capture.output( aln <- align() ))

  ali <- aln$ali
  ids <- aln$lab

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
    ##x[, aln.inds] <- ali[, positions, drop=FALSE]

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
      lapply(annot[aln.inds], function(x) span(x, class="aln_number")),
      lapply(annot[buf.inds2], function(x) span(x, class="aln_buff")),
      class="aln_row"
      )

    ## sequence row
    tmp.j <- prep.seq.row(tmp, i, x, x.mat, ids, buf.inds1, buf.inds2, aln.inds, width)
    tmp <- tmp.j$tmp
    j <- tmp.j$j
    rm(tmp.j)

    ## annotation row: conservation
    tmp[[j+2]] <- span(
      span(" ", class="aln_id"),
      lapply(annot2[buf.inds1],     function(x) span(x, class="aln_buff")),
      lapply(annot2[aln.inds], function(x) span(x, class="aln_conserv")),
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
####  Non-reactive functions    ####
####################################

`prep.seq.row` <- function(tmp, i, x, x.mat, ids, buf.inds1, buf.inds2, aln.inds, width){
    if( length(ids)<25 ) {
        for(j in 1:nrow(x)) {
          tmp[[j+1]] <- span(
            span(ids[j], class="aln_id"),
            lapply(x[j, buf.inds1], function(x) span(x, class="aln_buff", class=x)),
            lapply(seq_along(1:length(x.mat[j, aln.inds,1])), function(x) span(x.mat[j,aln.inds,1][x], class="aln_aminoacid", class=x.mat[j,aln.inds,1][x]
                , "title"="",
                "data-original-title"="Residue info","data-toggle"="popover",
                "data-content"=paste0(
                    "<table class=\"tb_pop\">
                    <tbody>
                        <tr><td>PDB ID:</td><td>",
                            ids[j],"</td></tr>",
                        "<tr><td>PDB residue number:</td><td>",
                            x.mat[j,aln.inds,2][x],"</td></tr>
                        <tr><td>Position in alignment:</td><td>",
                        j,", ",((i-1)*width)+x,"</td></tr>",
                    "</tbody>
                    </table>"
                )
                               )),
            lapply(x[j, buf.inds2], function(x) span(x, class="aln_buff", class=x)),
            class="aln_row"
            )
        }
    } else {
        for(j in 1:nrow(x)) {
            tmp[[j+1]] <- span(
                span(ids[j], class="aln_id"),
                lapply(x[j, buf.inds1], function(x) span(x, class="aln_buff", class=x)),
                lapply(seq_along(1:length(x.mat[j, aln.inds,1])), function(x) span(x.mat[j, aln.inds,1][x], class="aln_aminoacid", class=x.mat[j,aln.inds,1][x])),
                lapply(x[j, buf.inds2], function(x) span(x, class="aln_buff", class=x)),
                class="aln_row"
            )
        }
    }
    return(list(tmp=tmp, j=j))
}



####################################
####     Plotting functions     ####
####################################

make.plot.seqide.heatmap <- function() {
  pdbs <- align()
  ide <- seqide()
  rownames(ide) <- pdbs$lab
  colnames(ide) <- pdbs$lab
  
  hc <- hclust(as.dist(1-ide))
  grps <- cutree(hc, k=input$clusters_seq)
  mar <- as.numeric(c(input$margins0, input$margins0))
  plot1 <- heatmap(1-ide, distfun=as.dist, symm=TRUE,
          ColSideColors=as.character(grps),
          RowSideColors=as.character(grps),
          cexRow=input$cex0, cexCol=input$cex0,
          margins=mar
          )
  return(plot1)
}

output$seqide_heatmap <- renderPlot({
  print(make.plot.seqide.heatmap())
})

make.plot.seqide.dendrogram <- function() {
  pdbs <- align()
  ide <- seqide()
  
  hc <- hclust(as.dist(1-ide))
  mar <- c(input$margins0, 5, 3, 1)
    
  plot3 <- hclustplot(hc, k=input$clusters_seq, labels=pdbs$lab, cex=input$cex0,
                      ylab="Identity distance", main="Sequence identity clustering",
                      fillbox=FALSE, mar = mar)
  return(plot3)
}

output$seqide_dendrogram <- renderPlot({
  print(make.plot.seqide.dendrogram())
})


make.plot.conservation <- function() {
  pdbs <- align()
  sse <- pdbs2sse(pdbs, ind=1, rm.gaps=FALSE, exefile=configuration$dssp$exefile)
  x <- seqconserv()
  ylab <- paste("Sequence", input$conserv_method)

  ##mar <- as.numeric(c(input$margins0, input$margins0))
  p1 <- plot.bio3d(x, sse = sse,
                   ylab=ylab, xlab="Alignment Position",
                   cex=input$cex0)
  
  return(p1)
}

output$conservation <- renderPlot({
  print(make.plot.conservation())
})


####################################
####     Download functions     ####
####################################

aln2file <- reactive({
  path <- data_path()
  aln <- align()
  fn <- paste0(path, '/aln.fasta')
  write.fasta(aln, file=fn)
  return(fn)
})

output$fastafile = downloadHandler(
  filename = 'aln.fasta.zip',
  content = function(file) {
    zip(file, files=aln2file(), flags = "-9Xj")
  })

seqide2txt <- reactive({
  path <- data_path()
  pdbs <- align()
  ide <- round(seqide(), 2)

  file <- paste0(path, "/", "seqide.dat")
  write.table(ide, file=file, quote=FALSE)
  return(file)
})

output$seqideZIP = downloadHandler(
  filename = 'seqide.zip',
  content = function(file) {
    zip(file, files=seqide2txt(), flags = "-9Xj")
  })

output$seqide_heatmap2pdf = downloadHandler(
  filename = "ide_heatmap.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, onefile=T, width=input$width0, height=input$height0)
    print(make.plot.seqide.heatmap())
    dev.off()
})

output$seqide_dendrogram2pdf = downloadHandler(
  filename = "seqide_dendrogram.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, onefile=T, width=input$width0, height=input$height0)
    print(make.plot.seqide.dendrogram())
    dev.off()
  })

output$conservation2pdf = downloadHandler(
  filename = "seqconserv.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, onefile=T, width=input$width0, height=input$height0)
    print(make.plot.conservation())
    dev.off()
})
