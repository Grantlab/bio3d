#####################
##-- Align
####################

## selected accession ids
rv$selacc <- NULL

observeEvent(input$selected_pdbids, {
  rv$fitted <- FALSE
  if(!length(input$selected_pdbids) > 0) {
    rv$selacc <- get_acc()
  }
  
  if(!all(rv$selacc %in% input$selected_pdbids))
    rv$selacc <- input$selected_pdbids
  
  if(!all(input$selected_pdbids %in% rv$selacc))
    rv$selacc <- input$selected_pdbids
  
})


observeEvent(input$omit_missing, {
  rv$fitted <- FALSE
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
  rv$fitted <- FALSE
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
  message("split__pdbs called")
  
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
 rv$aligned <- TRUE
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

  if(!length(gaps$f.inds) > 5) {
    stop("Insufficient non-gap regions in alignment to proceed")
  }

  if(!nrow(aln$ali)>1)
    stop("Insufficient PDBs selected")
})

output$alignment_summary <- renderUI({
  invisible(capture.output( aln <- align() ))
  check_aln()

  id <- aln$id
  ali <- aln$ali

  nstruct <- length(id)
  dims <- dim(ali)
  gaps <- gap.inspect(ali)
  dims.nongap <- dim(ali[, gaps$f.inds, drop = FALSE])
  dims.gap <- dim(ali[, gaps$t.inds, drop = FALSE])

  str <- paste("<strong>Alignment dimensions:</strong><br>",
               paste0("<ul><li>", dims[1L], " sequence rows</li>"),
               paste0("<li>", dims[2L], " position columns</li>"),
               paste0("<li>(", dims.nongap[2L], " non-gap, ", dims.gap[2L], " gap)</li></ul>"),
               sep = "")
  if(dims[1L] <= 15) {
      updateRadioButtons(session, inputId = 'show_alignment', label = 'Show alignment',
                        choices = c('Show' = 'yes', 'Hide' = 'no'),
                        selected = 'yes', inline = TRUE)
  }
  shinyjs::runjs("$(window).scrollTop(0);")
  HTML(str)
})


output$omitted_pdbs_summary <- renderUI({
  invisible(capture.output( pdbs <- align() ))
  check_aln()
  
  wanted <- get_acc()
  
  missing <- !wanted %in% pdbs$lab
  if(any(missing))
    HTML("<strong>Omitted PDB(s):</strong>", paste(sort(wanted[missing]), collapse=", "))
})

pdbs_connectivity <- reactive({
  invisible(capture.output( pdbs <- align() ))
  check_aln()

  ids <- pdbs$lab
  nstructs <- length(ids)
  
  conn <- inspect.connectivity(pdbs, cut=4.05)
  return(conn)
})

output$missres_summary <- renderUI({
  invisible(capture.output( pdbs <- align() ))
  check_aln()

  ids <- pdbs$lab
  nstructs <- length(ids)

  conn <- pdbs_connectivity()
  if(sum(!conn) > 0) {
    inds <- which(!conn)
    
    str <- paste(paste0("<strong>", sum(!conn), " PDB(s) with missing in-structure residues:</strong><br>"),
                 #"<ul><li>", 
                 paste0(ids[!conn], collapse=", "),
                 #"</li></ul>", 
                 sep="")
    HTML(str)
  }
  #else {
  #  str <- "No PDBs with missing in-structure residues"
  #}
})

output$npdbs_with_missres <- reactive({
  if(!rv$aligned)
    return(0)
  
  conn <- pdbs_connectivity()
  return(sum(!conn))
})
outputOptions(output, 'npdbs_with_missres', suspendWhenHidden=FALSE)




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
  shinyjs::runjs('$("input[name=\'show_alignment\']").change(function(){
                 $("html, body").animate({scrollTop:$("#alignment_row").offset().top}, "smooth");
                 });')
  pre(class="alignment",
      out
      )
})

####################################
####  Non-reactive functions    ####
####################################

`prep.seq.row` <- function(tmp, i, x, x.mat, ids, buf.inds1, buf.inds2, aln.inds, width){
    if( length(ids)<50 ) {
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
####   Clustering functions     ####
####################################

hclust_seqide <- reactive({
  pdbs <- align()
  ide <- seqide()
  rownames(ide) <- pdbs$lab
  colnames(ide) <- pdbs$lab
  
  d <- as.dist(1-ide)
  h <- hclust(d, method=input$hclustMethod)
  
  h$labels <- pdbs$lab
  return(h)
})

cutree_seqide <- reactive({
  hc <- hclust_seqide()

  if(is.null(input$splitTreeK))
    k <- NA
  else
    k <- as.numeric(input$splitTreeK)
  
  cut <- cutreeBio3d(hc, minDistance=input$minDistance, k=k)
  return(cut)
})

observeEvent(input$setk, {
  hc <- hclust_seqide()
  cut <- cutreeBio3d(hc, minDistance=input$minDistance, k=NA)
  updateSliderInput(session, "splitTreeK", value = cut$autok)
})


output$kslider <- renderUI({
  cut <- cutree_seqide()
  sliderInput("splitTreeK", "Cluster/partition into K groups:",
              min = 1, max = 10, value = cut$k, step=1)
})
                    

####################################
####     Plotting functions     ####
####################################

output$schematic_alignment <- renderPlot({
  aln <- align()
  hc <- hclust_seqide()
  grps <- cutree_seqide()$grps

  if(input$cluster_alignment) {
    layout(matrix(c(4,2,3,1), ncol=2),
           heights = c(.1, 1),
           widths = c(0.3, 1))
    par(mar=c(4, 0.1, 0.1, 4))
  }
  else {
    layout(matrix(c(2, 1), nrow=2),
           heights = c(.1, 1))
    par(mar=c(4, 2, 0.1, 4))
  }
  
  ## 1: gap, 0: non-gap
  gaps <- gap.inspect(aln$ali)
  
  mat <- gaps$bin
  if(any(mat==1)) {
    mat[ mat == 1 ] <- -1
    mat[ mat == 0 ] <- 1
    mat[ mat == -1 ] <- 0
  }
  else {
    mat <- mat+1
  }
  
  ## re-order matrix
  if(input$cluster_alignment) 
    mat <- mat[ hc$order, ]
  else 
    mat <- mat[ seq(nrow(mat), 1), ]

  if(any(mat==0))
    col <- c("#FFFFFF", "#9F9F9F")
  else
    col <- c("#9F9F9F")
  
  image(t(mat), col = col, axes=FALSE)

  by <- pretty(0:ncol(aln$ali), n = 6)
  by <- by[2]
    
  labs <- seq(0, ncol(aln$ali), by=by)
  labs[1] <- 1

  at <- labs / labs[length(labs)]
  at[1] <- 0
  
  axis(1, at=at, labels=labs)
  mtext("Alignment index", 1, line=2, cex=1.0)
  at <- seq(0, 1, length.out=length(aln$lab))
  
  labs <- aln$lab
  if(input$cluster_alignment)
    labs <- labs[hc$order]
  
  mtext(labs, side=4, line=2-1.25, at=at, cex=0.8, las=2)



  ## cluster dendrogram
  if(input$cluster_alignment) {
    par(mar=c(4, 0.1, 0.1, 0.1))
    ddr <- as.dendrogram(hc)
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }


  ## conservation bar on top of alignment
  cons <- conserv(aln$ali)
  ng <- rep(0, length(cons))
  ng[gaps$f.inds] <- 1
 

  if(input$cluster_alignment)
    par(mar=c(.1, 0.1, 2, 4))
  else
    par(mar=c(.1, 2, 2, 4))

  #if(all(c(0,1) %in% cons))
  #  col <- c("#FFFFFF", "#FF0000")
  #else 
  col <- colorRampPalette(c("white", "red"))( 10 )

  image(as.matrix(cons), col = col, axes=FALSE)

  if(input$cluster_alignment)
    at <- 0.4
  else
    at <- NA
  
  mtext(quote(bold("Sequence Alignment Overview")), side = 3, line = 0.5, cex = 1.25, at=at)
})


make.plot.seqide.heatmap <- function() {
  pdbs <- align()
  ide <- seqide()
  rownames(ide) <- pdbs$lab
  colnames(ide) <- pdbs$lab
  
  hc <- hclust_seqide()
  grps <- cutree_seqide()$grps
  
  mar <- as.numeric(c(input$margins0, input$margins0))
  plot1 <- heatmap(1-ide,
                   hclustfun = function(x) {
                     hclust(x, method=input$hclustMethod)
                   }, 
                   distfun=as.dist, symm=TRUE,
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

  mar <- c(input$margins0, 5, 3, 1)
  hc <- hclust_seqide()
  cut <- cutree_seqide()
  
  if(vec_is_sorted(hc$height)) {
    max_branch_gap <- max(diff(hc$height))
    k <- length(unique(cut$grps))
    
    hclustplot(hc, k=k, labels=pdbs$lab, cex=input$cex0,
               ylab="Identity distance", main="Sequence identity clustering",
               fillbox=FALSE, mar = mar)

    #par("lty" = 2); par("lwd" = 2); 
    #rect.hclust(hc, k=k, cluster = cut$grps, border = "grey50")

    ## red line only if auto k 
    if(max_branch_gap >= input$minDistance & cut$autok == cut$k ) {
      abline(h = cut$h, col="red", lty=2)
    }
  } else {
    plot(hc, main="", xlab="", sub="")
  }
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
