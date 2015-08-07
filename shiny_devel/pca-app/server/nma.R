### reactive stuff
init_show_trj_nma <- TRUE

rv$modes <- NULL
output$nmaIsDone <- reactive({
  return(!is.null(rv$modes))
})
outputOptions(output, 'nmaIsDone', suspendWhenHidden=FALSE)

observeEvent(input$run_nma, {
  rv$modes <- nma2()
})

observeEvent(input$filter_rmsd, {
  rv$modes <- NULL
})

############
## NMA
###########


nma2 <- reactive({
  ##pdbs <- align()
  pdbs <- pdbs4nma()
    
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  progress$set(message = 'Calculating normal modes',
               detail = 'Please wait',
               value = 0)

  rm.gaps <- TRUE
  if(is.logical(input$rm_gaps)) {
    rm.gaps <- input$rm_gaps
  }
  else {
    rm.gaps <- TRUE
  }
    
  modes <- nma(pdbs, fit=TRUE, rm.gaps=rm.gaps, progress=progress)
  
  if(init_show_trj_nma) {
    updateCheckboxInput(session, 'show_trj2', 'Show NM Trajectory', value=TRUE)
    init_show_trj_nma <<- FALSE
  }
  return(modes)
})



####################################
####   Filtering structures #######
####################################


filter_by_rmsd <- reactive({
  pdbs <- align()
  rd <- rmsd1()
  
  inds <- filter.rmsd(pdbs$xyz, rmsd.mat = rd,
                      cutoff=input$filter_rmsd, fit=FALSE,
                      method = input$hclustMethod_rmsd)$ind
  
  return(inds)
})

pdbs4nma <- reactive({
  pdbs1 <- align()
  inds <- filter_by_rmsd()
  pdbs2 <- trim.pdbs(pdbs1, row.inds=inds)
  pdbs2$lab <- pdbs1$lab[inds]

  message(pdbs2$lab)
  return(pdbs2)
})

output$rmsd_dendrogram2 <- renderPlot({
  ##pdbs <- pdbs_nma()
  pdbs <- align()
  rd <- rmsd1()
  rownames(rd) <- pdbs$lab
  colnames(rd) <- pdbs$lab

  hc <- hclust_rmsd()
  #d <- as.dist(rd)
  #hc <- hclust(d)

  inds <- filter_by_rmsd()
  message(inds)
  
  col <- rep(2, length(hc$order))
  col[inds] <- 1

  mar <- c(input$margins, 5, 3, 1)
  hclustplot(hc, k=1, col=col, main="RMSD Cluster Dendrogram")
  
})

output$filter_summary <- renderPrint({
  pdbs1 <- align()
  inds <- filter_by_rmsd()

  cat("Structures included:", length(inds), "\n")
  cat("Structures excluded:", length(pdbs1$id) - length(inds), "\n")
  
  
})


output$filter_rmsd <- renderUI({
  cut <- cutree_rmsd()
  numericInput("filter_rmsd","RMSD Cutoff", value = round(cut$h, 2), step = 0.25)
})




####################################
####   similarity measures     ####
####################################

rmsip2 <- reactive({
  if(input$rm_gaps) {
    ##pdbs <- align()
    pdbs <- pdbs4nma()
    ##modes <- nma2()
    modes <- rv$modes
    rownames(modes$rmsip) <- pdbs$lab
    colnames(modes$rmsip) <- pdbs$lab
    return(modes$rmsip)
  }
  else {
    stop("RMSIP only available with 'Omit gaps'")
  }
})


bhat <- reactive({
  if(input$rm_gaps) {
    ## Bhattacharyya coefficient
    pdbs <- align()
    ##modes <- nma2()
    modes <- rv$modes
    covs <- cov.enma(modes)
    bc <- bhattacharyya(modes, covs=covs)
    rownames(bc) <- pdbs$lab
    colnames(bc) <- pdbs$lab
    return(bc)
  }
  else {
    stop("RMSIP only available with 'Omit gaps'")
  }
})

####################################
####   Clustering functions     ####
####################################

hclust_rmsip <- reactive({
  rd <- rmsip2()
  return(hclust(as.dist(1-rd),
         method=input$hclustMethod_rmsip))
})

hclust_bhat2 <- reactive({
  rd <- bhat2()
  return(hclust(as.dist(1-rd)))
})

hclust2 <- reactive({
  if(input$group_by2 == "rmsd") {
    hc <- hclust_rmsd()
  }
  
  if(input$group_by2 == "rmsip") {
    hc <- hclust_rmsip()
  }

  if(input$group_by2 == "bhat") {
    hc <- hclust_bhat()
  }

  if(input$group_by2 == "pc_space") {
    hc <- hclust_pcspace()
  }
  
  return(hc)
})


cutree_rmsip <- reactive({
  hc <- hclust_rmsip()
  
  if(is.null(input$splitTreeK_rmsip))
    k <- NA
  else
    k <- as.numeric(input$splitTreeK_rmsip)
  
  cut <- cutreeBio3d(hc, minDistance=input$minDistance_rmsip, k=k)
  return(cut)
})


cutree_nma <- reactive({
  inds <- filter_by_rmsd()
    
  if(input$cluster_by2 == "rmsd") {
    cut <- cutree_rmsd()
    cut$grps <- cut$grps[inds]
  }

  if(input$cluster_by2 == "pc_space") {
    cut <- cutree_pca()
    cut$grps <- cut$grps[inds]
  }

  if(input$cluster_by2 == "rmsip") {
    cut <- cutree_rmsip()
  }
  
  if(input$cluster_by2 == "sequence") {
    cut <- cutree_seqide()
    cut$grps <- cut$grps[inds]
  }
 
  return(cut)
})


observeEvent(input$setk_rmsip, {
  hc <- hclust_pca()
  cut <- cutreeBio3d(hc, minDistance=input$minDistance_rmsip, k=NA)
  updateSliderInput(session, "splitTreeK_rmsip", value = cut$autok)
})


output$kslider_rmsip <- renderUI({
  cut <- cutree_pca()
  sliderInput("splitTreeK_rmsip", "Cluster/partition into K groups:",
              min = 1, max = 10, value = cut$k, step=1)
})



####################################
####     webGL functions        ####
####################################

output$struct_dropdown2 <- renderUI({
  ##pdbs <- align()
  pdbs <- pdbs4nma()
  ids <- 1:length(pdbs$id)
  names(ids) <-  pdbs$lab
  selectInput('viewStruct_nma', 'Show NMs for structure:',
              choices=ids)
})
  
output$nmaWebGL  <- renderWebGL({
  pdbs <- pdbs4nma()
  modes <- rv$modes

  if(!inherits(modes, "enma"))
    return(NULL)

  mag <- as.numeric(input$mag2)
  step <- mag/8
  
  trj <- mktrj(modes, pdbs=pdbs,
               mag=mag, step=step,
               s.inds=as.numeric(input$viewStruct_nma),
               m.inds=as.numeric(input$viewMode_nma),
               rock=FALSE,
               file = paste0(data_path(), '/mode_', input$viewMode_nma, '.pdb'))
  n <- nrow(trj)
  
  amalcol <- function(x) {
    col <- rep("grey50", length(x))
    col[1] <- "blue"
    col[length(col)] <- "red"
    return(col)
  }
  
  magcol <- function() {
    rf <- rmsf(trj)
    return(t(replicate(n, vec2color(rf, c('blue', 'red')),simplify=TRUE)))
  }

  class(trj)  <- 'xyz'
  col <- switch(input$viewColor2,
                'mag' = magcol(), # vec2color(rmsf(m)), #!! col=col, type=2
                'amalgam' = amalcol(1:n),
                'default' = colorRampPalette(c('blue', 'gray', 'red'))(n)
                )
  
    view.xyz(trj, bg.col=input$viewBGcolor2, col=col, d.cut=6)
})



####################################
####     Plotting functions     ####
####################################

output$checkboxgroup_label_ids2 <- renderUI({
  ##pdbs <- align()
  pdbs <- pdbs4nma()
  grps <- cutree_nma()$grps
  ids <- pdbs$lab
  names(ids) <- paste(ids, " (c", grps[order(grps)], ")", sep="")

  checkboxInput("toggle_all2", "Toggle all", TRUE)
  
  if(input$toggle_all2) {
    checkboxGroupInput(inputId="label_ids2", label="PDB IDs:",
                       choices=ids, selected=ids, inline=TRUE)
  }
  else {
    checkboxGroupInput(inputId="label_ids2", label="PDB IDs:",
                       choices=ids, selected=c(), inline=TRUE)
  }
})


## Fluctuation plot
make.plot.nma <- function() {
  ##pdbs_all <- align()
  pdbs <- pdbs4nma()
  ##modes <- nma2()
  modes <- rv$modes

  sse <- pdbs2sse(pdbs, ind=1, rm.gaps=TRUE, exefile=configuration$dssp$exefile)

  signif <- FALSE
  if(input$cluster) {
    col <- cutree_nma()$grps

    #if(input$signif)
    #  signif <- TRUE
  }
  else {
    col <- 1:length(pdbs$id)
  }

  if(length(input$label_ids2) > 0) {
    inds <- unlist(lapply(input$label_ids2, grep, pdbs$lab))
    show <- rep(FALSE, length(pdbs$lab))
    show[inds] <- TRUE
    col[!show] <- NA
  }

  plot(modes, pdbs, sse=sse, col=col, signif=signif,
       spread=input$spread, conservation=input$seqide)

  #if(input$toggle_rmsf2) {
  #  gaps.res <- gap.inspect(pdbs$ali)
  #  gaps.pos <- gap.inspect(pdbs$xyz)
  #  rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])
  #  
  #  par(new=TRUE)
  #  plot5 <- plot.bio3d(rf, axes=FALSE, col=2, type="l", xlab="", ylab="")
  #  axis(4, col=2)
  #}


}


output$nma_fluctplot <- renderPlot({
   print(make.plot.nma())
 })


make.plot.heatmap_rmsd2 <- function() {
  inds <- filter_by_rmsd()
  
  pdbs <- fit()
  rd <- rmsd1()
  rownames(rd) <- pdbs$lab
  colnames(rd) <- pdbs$lab

  hc <- hclust_rmsd()
  cut <- cutree_rmsd()
  grps <- cut$grps[inds]
  rd <- rd[inds, inds]
  
  mar <- as.numeric(c(input$margins3, input$margins3))
  heatmap(rd,
          hclustfun = function(x) {
            hclust(x, method=input$hclustMethod_rmsd)
          }, 
          distfun=as.dist, symm=TRUE,
          ColSideColors=as.character(grps),
          RowSideColors=as.character(grps),
          cexRow=input$cex3, cexCol=input$cex3,
          margins=mar
          )
}

output$rmsd_heatmap2 <- renderPlot({
  make.plot.heatmap_rmsd2()
})

make.plot.heatmap_rmsip2 <- function() {
  rp <- rmsip2()
  hc <- hclust_rmsip()
  cut <- cutree_rmsip()
  grps <- cut$grps
  
  mar <- as.numeric(c(input$margins3, input$margins3))
  heatmap(1-rp,
          hclustfun = function(x) {
            hclust(x, method=input$hclustMethod_rmsip)
          }, 
          distfun=as.dist, symm=TRUE,
          ColSideColors=as.character(grps),
          RowSideColors=as.character(grps),
          cexRow=input$cex3, cexCol=input$cex3,
          margins=mar)
}

output$rmsip_heatmap2 <- renderPlot({
  make.plot.heatmap_rmsip2()
})


output$bhat_heatmap2 <- renderPlot({
  bc <- bhat2()
  heatmap(1-bc, distfun=as.dist, symm=TRUE)
})


make.plot.dendrogram2 <- function() {
  ##pdbs <- align()
  pdbs <- pdbs4nma()
  pc <- pca1()

  hc <- hclust_rmsip()
  cut <- cutree_rmsip()
  grps <- cut$grps
  k <- cut$k
  ids <- pdbs$lab

  mar <- c(input$margins2, 5, 3, 1)
  main <- paste(toupper(input$group_by2), "dendrogram")
  hclustplot(hc, k=k, colors=grps, labels=ids, cex=input$cex2,
             main=main, fillbox=FALSE, mar=mar)

  if(input$show_confplot2) {
    ### colored by RMSD here
    grps <- cutree_rmsd()$grps
      
    xlim <- range(pc$z[,1])
    ylim <- range(pc$z[,2])

    par(fig=c(.65, 1, .45, 1), new = TRUE)
    plot(pc$z[,1:2], col="grey20", pch=16, cex=1.1*input$cex2, 
         ylab="", xlab="", axes=FALSE,
         xlim=xlim, ylim=ylim)
    rect(xlim[1]*1.2, ylim[1]*1.2, xlim[2]*1.2, ylim[2]*1.2, col="grey90")
    points(pc$z[,1:2], col="grey20", pch=16, cex=1.2*input$cex2)
    points(pc$z[,1:2], col=grps, pch=16, cex=0.8*input$cex2)
    box()
    mtext("PC conformer plot", side=3, line=0, adj=0, cex=0.9)
  }
}

output$dendrogram2 <- renderPlot({
  print(make.plot.dendrogram2())
})

####################################
####     Download functions     ####
####################################

output$nmaplot2pdf = downloadHandler(
  filename = "nma_fluctuations.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, width=input$width1, height=input$height1)
    make.plot.nma()
    dev.off()
})

output$nmadendrogram2pdf = downloadHandler(
  filename = "nma_dendrogram.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, width=input$width2, height=input$height2)
    make.plot.dendrogram2()
    dev.off()
})

output$nma_rmsd_heatmap2pdf = downloadHandler(
  filename = "rmsd_heatmap.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, width=input$width3, height=input$height3)
    make.plot.heatmap_rmsd2()
    dev.off()
})

output$nma_rmsip_heatmap2pdf = downloadHandler(
  filename = "rmsip_heatmap.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, width=input$width3, height=input$height3)
    make.plot.heatmap_rmsip2()
    dev.off()
})

####################################
####     Download trajectory    ####
####################################

nma2pdb  <- reactive({
  path <- data_path()
  pdbs <- align()
  ##modes <- nma2()
  modes <- rv$modes
  gaps <- gap.inspect(pdbs$ali)
  fname  <- paste0(path, '/', 'mode', input$viewMode_nma, '.pdb')
  
  mag <- as.numeric(input$mag2)
  step <- mag/8

  trj <- mktrj(modes, pdbs=pdbs,
               mag=mag, step=step,
               s.inds=as.numeric(input$viewStruct_nma),
               m.inds=as.numeric(input$viewMode_nma),
               file=fname)
    
  return(fname)
})

output$nmtraj = downloadHandler(
  filename=function() {
    paste0('mode', input$viewMode_nma, '.pdb.zip')
  },
  content=function(file) {
    ## Avoid possibility of not having write permission on server
    #src  <- normalizePath('nma-traj.pdb')
    #owd <- setwd(tempdir())
    #on.exit(setwd(owd))
    #file.copy(src, 'nma-traj.pdb')
    #file.rename(trj2pdb2(), file)
    zip(file, files=nma2pdb(), flags = "-9Xj")
  }
  )
