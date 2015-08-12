#####################
##-- Superimpose and RMSD stuff
####################

init_show_pdbs <- TRUE
find_core <- reactive({
  pdbs <- align()
  check_aln()
  gaps <- gap.inspect(pdbs$ali)

  if(length(gaps$f.inds) <= 15)
    stop.at <- length(gaps$f.inds)-1
  else
    stop.at <- 15

  progress <- shiny::Progress$new()
  on.exit(progress$close())

  progress$set(message = 'Finding core',
               detail = 'Please wait ...',
               value = 0)

  core <- core.find(pdbs, stop.at=stop.at, progress=progress)
  progress$close()
  return(core)
})

fit <- reactive({
  pdbs <- align()
  check_aln()
  core <- find_core()

  if(input$fit_type == "full") {
    pdbs$xyz <- pdbfit(pdbs)
  }
  else {
    core <- find_core()
    pdbs$xyz <- pdbfit(pdbs, core)
  }
  
  if(init_show_pdbs) {
    updateCheckboxInput(session, 'show_pdbs', 'Show PDBs', value=TRUE)
    init_show_pdbs <<- FALSE
  }

  if(configuration$pdbdir$archive) {
    pdbs$lab <- pdbfilename2label(pdbs$id)
  }
  else {
    pdbs$lab <- basename.pdb(pdbs$id)
  }
  pdbs$lab <- format_pdbids(pdbs$lab)
  rv$fitted <- TRUE
  return(pdbs)
})

rmsd1 <- reactive({
  pdbs <- fit()

  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  progress$set(message = 'Calculating RMSD',
               detail = 'Please wait ...',
               value = 0)
  rd <- rmsd(pdbs, fit=TRUE)
  progress$close()
  
  return(rd)
})


####################################
####   Clustering functions     ####
####################################

hclust_rmsd <- reactive({
  pdbs <- fit()
  rd <- rmsd1()
  d <- as.dist(rd)
  h <- hclust(d, method=input$hclustMethod_rmsd)
  
  h$labels <- pdbs$lab
  return(h)
})

cutree_rmsd <- reactive({
  hc <- hclust_rmsd()

  if(is.null(input$splitTreeK_rmsd))
    k <- NA
  else
    k <- as.numeric(input$splitTreeK_rmsd)
  
  cut <- cutreeBio3d(hc, minDistance=input$minDistance_rmsd, k=k)
  return(cut)
})

observeEvent(input$setk_rmsd, {
  hc <- hclust_rmsd()
  cut <- cutreeBio3d(hc, minDistance=input$minDistance_rmsd, k=NA)
  updateSliderInput(session, "splitTreeK_rmsd", value = cut$autok)
})


output$kslider_rmsd <- renderUI({
  cut <- cutree_rmsd()
  sliderInput("splitTreeK_rmsd", "Cluster/partition into K groups:",
              min = 1, max = 10, value = cut$k, step=1)
})


representatives <- reactive({
  pdbs <- fit()
  rd <- rmsd1()
  hc <- hclust_rmsd()
  grps <- cutree_rmsd()$grps

  unq <- unique(grps)
  all.inds <- seq(1, length(grps))
  rep.inds <- rep(NA, length(unq))

  for(i in 1:length(unq)) {
    grp <- unq[i]
    inds <- which(grps==grp)

    if(length(inds) > 1) {
      tmp <- rd[inds, inds]

      m <- rowMeans(tmp)
      rep.ind <- which.min(m)
      rep.ind <- all.inds[inds][rep.ind]
    }
    else {
      rep.ind <- inds[1]
    }

    rep.inds[i] <- rep.ind
  }
  out <- data.frame(pdbid=pdbs$lab[rep.inds], clusterid=grps[rep.inds])
  rownames(out) <- NULL
  return(out)
})


####################################
####     Plotting functions     ####
####################################


#output$qtl_rmsd_heatmap <- iplotCorr_render({
#  pdbs <- fit()
#  rd <- rmsd1()
#  rownames(rd) <- pdbs$lab
#  colnames(rd) <- pdbs$lab
#  hc <- hclust(as.dist(rd))
#  grps <- cutree(hc, k=input$clusters)
#  iplotCorr(rd)
#})


make.plot.heatmap <- function() {
  pdbs <- fit()
  rd <- rmsd1()
  rownames(rd) <- pdbs$lab
  colnames(rd) <- pdbs$lab

  hc <- hclust_rmsd()
  cut <- cutree_rmsd()

  grps1 <- cut$grps
  grps2 <- grps1

  ## cluster by seq ide
  if(input$rowcol_seqide) {
    grps2 <- cutree_seqide()$grps
  }

  ## plot it
  mar <- as.numeric(c(input$margins, input$margins))
  plot1 <- heatmap(rd,
                   hclustfun = function(x) {
                     hclust(x, method=input$hclustMethod_rmsd)
                   }, 
                   distfun=as.dist, symm=TRUE,
                   ColSideColors=as.character(grps1),
                   RowSideColors=as.character(grps2),
                   cexRow=input$cex, cexCol=input$cex,
                   margins=mar
                   )
  return(plot1)
}

output$rmsd_heatmap <- renderPlot({
  print(make.plot.heatmap())
})

make.plot.rmsd.hist <- function() {
  pdbs <- fit()
  rd <- rmsd1()
  x <- rd[upper.tri(rd)]
  plot2 <- hist(x, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")
  ##, freq=FALSE)
  ##lines(density(x), col="gray", lwd=3)
  ##return(plot2)
}

output$rmsd_hist <- renderPlot({
  print(make.plot.rmsd.hist())
})

make.plot.rmsd.dendogram <- function() {
  pdbs <- fit()
  rd <- rmsd1()

  mar <- c(input$margins, 5, 3, 1)
  hc <- hclust_rmsd()
  cut <- cutree_rmsd()
  
  if(vec_is_sorted(hc$height)) {
    max_branch_gap <- max(diff(hc$height))
    k <- cut$k
    
    hclustplot(hc, k=k, labels=pdbs$lab, cex=input$cex,
               ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE,
               mar = mar)
    
    if(max_branch_gap >= input$minDistance_rmsd && cut$k == cut$autok) {
      abline(h = cut$h, col="red", lty=2)
    }
  } else {
    plot(hc, main="", xlab="", sub="")
  }
}

output$rmsd_dendrogram <- renderPlot({
  print(make.plot.rmsd.dendogram())
})

make.plot.rmsf <- function(){
  pdbs <- fit()
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  resno <- pdbs$resno[1, gaps.res$f.inds]
  sse <- pdbs2sse(pdbs, ind=1, rm.gaps=TRUE, exefile=configuration$dssp$exefile)

  rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])
  plot4 <- plot.bio3d(rf, resno=resno, sse=sse,
             ylab="RMSF (Å)", xlab="Residue No.")

  return(plot4)
}

output$rmsf_plot <- renderPlot({
  print(make.plot.rmsf())
})

####################################
####     Data summary           ####
####################################


output$representatives <- DT::renderDataTable({
  pdbs <- fit()
  reps <- representatives()

  DT::datatable(reps, extensions = 'Scroller', escape = TRUE,
                selection = 'none',
                colnames = c("ID", "Cluster ID"),
                
                options = list(
                  deferRender = FALSE,
                  dom = "frtiS",
                  scrollY = 200,
                  scrollCollapse = TRUE)
                )
})




output$print_core <- DT::renderDataTable({
  pdbs <- fit()
  core <- find_core()
  gaps <- gap.inspect(pdbs$ali)

  if(!is.null(input$reference_id))
    ind <- grep(input$reference_id, pdbs$lab)
  else
    ind <- 1
  
  resid <- pdbs$resid[ind, core$atom]
  resno <- pdbs$resno[ind, core$atom]
  bds <- bounds(resno)

  resid.inds <- which(resno %in% bds[,"start"])
  start <- paste0(resid[resid.inds], bds[, "start"])

  resid.inds <- which(resno %in% bds[,"end"])
  end <- paste0(resid[resid.inds], bds[, "end"])

  out <- data.frame(start, end, length=bds[, "length"])
  DT::datatable(out, extensions = 'Scroller', escape = TRUE,
                selection = 'none',
                colnames = c("Start", "End", "Length"),
                
      options = list(
                deferRender = FALSE,
                dom = "frtiS",
                scrollY = 200,
                scrollCollapse = TRUE)
  )
})

output$reference_selector <- renderUI({
  pdbs <- fit()
  ids <- pdbs$lab
  names(ids) <- ids
  selectInput("reference_id", "Reference PDB id:", ids)
})



output$rmsd_table <- DT::renderDataTable({
  pdbs <- fit()
  rd <- rmsd1()
  
  hc <- hclust_rmsd()
  grps <- cutree_rmsd()$grps

  if(!is.null(input$reference_id))
    ind <- grep(input$reference_id, pdbs$lab)
  else
    ind <- 1

  gaps <- gap.inspect(pdbs$xyz)
  rd <- rmsd(pdbs$xyz[ind, gaps$f.inds],
             pdbs$xyz[, gaps$f.inds],
             fit = TRUE)
  names(rd) <- pdbs$lab

  x <- data.frame(ids=names(rd), rmsd=round(rd,1), "cluster"=grps)
  DT::datatable(x, extensions = 'Scroller', escape = TRUE,
                selection = 'none', rownames = FALSE,
                colnames = c("ID", "RMSD", "Cluster ID"),
                  
      options = list(
        deferRender = FALSE,
        dom = 'frtiS',
        scrollY = 250,
        scrollCollapse = TRUE
        )
      )
})


####################################
####     Web GL                 ####
####################################

output$show_structs <- renderUI({

  acc <- get_acc()
  names(acc) <- acc
  
  selectInput("show_pdbids", "Visualize PDB IDs",
              choices = acc, selected = acc, multiple=TRUE)

})


output$pdbs_table1 <- renderDataTable({
  datatable(get_pdbstable(), extensions = 'Scroller', escape = FALSE,
            colnames = c("ID", "Cluster", "Length", "Name", "Species", "Ligands"),
            ##selection = "none",
            options = list(
              deferRender = TRUE,
              dom = "frtiS",
              scrollY = 200,
              scrollCollapse = TRUE,
              autoWidth = FALSE,
              columnDefs = list(list(width = '40%', targets = c(list(3))))
              ))
})


##
##-- PDBs Viewing Tab
##
output$pdbsWebGL  <- renderWebGL({
  pdbs <- fit() ## always used for visualization (can be trimmed)
  
  col.inds <- NULL  ## only set when pdbs object is trimmed
  pdbs.full <- NULL ## only set when a trimmed version occupies the pdbs object
  
  if(!is.null(input$pdbs_table1_rows_selected)) {
    inds <- as.numeric(input$pdbs_table1_rows_selected)
    
    gaps <- gap.inspect(pdbs$resno[inds, , drop = FALSE])
    col.inds <- which(gaps$col < dim(pdbs$resno[inds, , drop = FALSE])[1L])

    pdbs.full <- pdbs ## original pdbs object
    pdbs <- trim(pdbs, row.inds=inds)
  }
  else {
    inds <- 1:length(pdbs$id)  
  }

  n <- nrow(pdbs$xyz)
  grps <- cutree_rmsd()$grps[inds]
  core <- find_core()
  
  ## dims and gaps always refer to the full matrix
  if(!is.null(pdbs.full)) {
    gaps <- gap.inspect(pdbs.full$ali)
    dims <- dim(pdbs.full$ali)
  }
  else {
    gaps <- gap.inspect(pdbs$ali)
    dims <- dim(pdbs$ali)
  }
  
  
  corecol <- function() {
    col <- matrix(1, nrow=dims[1], ncol=dims[2])
    col[, core$atom] <- 2
    
    if(!is.null(pdbs.full)) 
      col <- col[, col.inds, drop=FALSE]
    
    col <- col[inds, , drop=FALSE]
    return(col)
  }
  
  gapscol <- function() {
    col <- matrix(1, nrow=dims[1], ncol=dims[2])
    if(length(gaps$t.inds)>0)
      col[, gaps$t.inds] <- 2
    
    if(!is.null(pdbs.full)) 
      col <- col[, col.inds, drop=FALSE]

    col <- col[inds, , drop=FALSE]
    return(col)
  }
  
  rmsfcol <- function() {
    if(n==1 & input$color_recalc)
      return(rep("#808080", ncol(pdbs$ali)))
    
    rf <- rmsf(pdbs$xyz)

    ## trimmed ensemble - do not recalc
    if(!input$color_recalc & !is.null(pdbs.full)) {
      rf <- rmsf(pdbs.full$xyz)[ col.inds ]
    }
    
    return(vec2color(rf))
  }

  framecol <- function() {
    if(!input$color_recalc & !is.null(pdbs.full)) {
      nstru <- nrow(pdbs.full$xyz)
      npos  <- ncol(pdbs.full$xyz)/3
      
      col <- matrix( rep(vec2color(1:nstru), times=npos), ncol=npos)
      col <- col[inds, col.inds, drop=FALSE]
    }
    else {
      nstru <- nrow(pdbs$xyz)
      npos  <- ncol(pdbs$xyz)/3

      if(nstru == 1)
        return("#0000FF")
      else
        col <- matrix( rep(vec2color(1:nstru), times=npos), ncol=npos)
    }

    return(col)
  }

  indexcol <- function() {
    if(!input$color_recalc & !is.null(pdbs.full)) {
      nstru <- nrow(pdbs.full$xyz)
      npos  <- ncol(pdbs.full$xyz)/3
      col <- matrix( rep(vec2color(1:npos), times=nstru), nrow=nstru, byrow=TRUE)
      col <- col[inds, col.inds, drop=FALSE]
    }
    else {
      return("index")
    }
  }
    
  col <- switch(input$viewColor1,
                'index' = indexcol(),
                'frame' = framecol(),
                'sse'   = 'sse',
                'cluster' = grps,
                'core'  = corecol(),
                'gaps'  = gapscol(),
                'rmsf'  = rmsfcol()
                )

  bg <- input$viewBGcolor1
  if(bg == "black") {
    col[ col==1 ] <- "grey90"
  }
  
  view(pdbs, bg.col=bg, col=col, maxframes=500, sheet="blue")
})

# observeEvent(input$viewUpdate1, {
#   updateRadioButtons(session, 'viewColor1', label='Structure color',
#                      choices=list(
#                        'By cluster ID'='cluster',
#                        'By structure ID'='struct',
#                        'Invariant core'='core',
#                        'Gap regions'='gaps'
#                        ),
#                      selected='cluster')

#   updateRadioButtons(session, 'viewBGcolor1', label='Background color',
#                      choices=list('Black'='black', 'White'='white'),
#                      selected='white')
# })



####################################
####     Download functions     ####
####################################

pdbs2rda <- reactive({
  path <- data_path()
  core <- find_core()
  pdbs <- fit()
  xyz <- pdbfit(pdbs, core)
  fn <- paste0(path, "/", "pdbs.RData")
  save(pdbs, file=fn)
  return(fn)
})

output$pdbsRData = downloadHandler(
  filename = 'pdbs.zip',
  content = function(file) {
    zip(file, files=pdbs2rda(), flags = "-9Xj")
  })


pdbs2files <- reactive({
  path <- data_path()
  core <- find_core()
  pdbs <- fit()
  xyz <- pdbfit(pdbs, core, outpath=path)
  files <- paste0(path, "/", pdbs$lab, ".pdb_flsq.pdb")
  
  ## rename files from pdb1xck_A.ent.gz to 1XCK_A
  if(configuration$pdbdir$archive) {
    files.from <- paste0(path, "/", basename(pdbs$id), "_flsq.pdb")
    file.rename(files.from, files)
  }
  return(files)
})

output$pdbsZIP = downloadHandler(
  filename = 'pdbs.zip',
  content = function(file) {
    zip(file, files=pdbs2files(), flags = "-9Xj")
})

rmsd2files <- reactive({
  path <- data_path()
  pdbs <- fit()
  rd <- round(rmsd1(), 2)
  rownames(rd) <- pdbs$lab
  colnames(rd) <- pdbs$lab
  file <- paste0(path, "/", "rmsd-mat.dat")
  write.table(rd, file=file, quote=FALSE)
  return(file)
})

output$rmsdZIP = downloadHandler(
  filename = 'rmsd-mat.zip',
  content = function(file) {
    zip(file, files=rmsd2files(), flags = "-9Xj")
})

output$rmsd_heatmap2pdf = downloadHandler(
  filename = "rmsd_heatmap.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, onefile=T, width=input$width, height=input$height)
    print(make.plot.heatmap())
    dev.off()
})

output$rmsd_dendrogram2pdf = downloadHandler(
  filename = "rmsd_dendrogram.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, onefile=T, width=input$width, height=input$height)
    print(make.plot.rmsd.dendogram())
    dev.off()
})

output$rmsd_hist2pdf = downloadHandler(
  filename = "rmsd_hist.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, onefile=T, width=input$width, height=input$height)
    print(make.plot.rmsd.hist())
    dev.off()
})

output$rmsf2pdf = downloadHandler(
  filename = "rmsf.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, onefile=T, width=input$width, height=input$height)
    print(make.plot.rmsf())
    dev.off()
})

make_pdbs_pse <- reactive({
  path <- data_path()
  core <- find_core()
  pdbs <- fit()

  ## rename files from pdb1xck_A.ent.gz to 1XCK_A
  if(configuration$pdbdir$archive) {
    files.from <- pdbs$id
    files.to <- paste0(path, "/", pdbs$lab, ".pdb")
    file.copy(files.from, files.to)
    pdbs$id <- files.to
  }
  
  col <- switch(input$viewColor1,
                "cluster" = cutree_rmsd()$grps,
                "struct" = NULL,
                "core" = core,
                "gaps" = "gaps")

  outf <- paste0(path, "/pdbs.pse")
  file <- pymol(pdbs, col=col, type="session", file=outf)
  return(outf)
})

output$pdbs2pymol = downloadHandler(
  filename = 'pdbs.pse.zip',
  content = function(file) {
    zip(file, files=make_pdbs_pse(), flags = "-9Xj")
})
