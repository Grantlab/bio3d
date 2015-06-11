############
## NMA
###########
init_show_trj_nma <- TRUE

nma2 <- reactive({
  pdbs <- align()
  check_aln()
    
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  progress$set(message = 'Calculating normal modes',
               detail = 'Please wait',
               value = 0)

  rm.gaps <- TRUE
  if(is.logical(input$rm.gaps)) {
    rm.gaps <- input$rm.gaps
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

rmsip2 <- reactive({
  if(input$rm.gaps) {
    pdbs <- align()
    modes <- nma2()
    rownames(modes$rmsip) <- pdbs$lab
    colnames(modes$rmsip) <- pdbs$lab
    return(modes$rmsip)
  }
  else {
    stop("RMSIP only available with 'Omit gaps'")
  }
})

hclust_rmsip2 <- reactive({
  rd <- rmsip2()
  return(hclust(as.dist(1-rd)))
})


rmsd2 <- reactive({
  pdbs <- align()
  rd <- rmsd(pdbs, fit=TRUE)
  rownames(rd) <- pdbs$lab
  colnames(rd) <- pdbs$lab
  return(rd)
})

hclust_rmsd2 <- reactive({
  rd <- rmsd2()
  return(hclust(as.dist(rd)))
})

bhat2 <- reactive({
  if(input$rm.gaps) {
    ## Bhattacharyya coefficient
    pdbs <- align()
    modes <- nma2()
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

hclust_bhat2 <- reactive({
  rd <- bhat2()
  return(hclust(as.dist(1-rd)))
})


hclust2 <- reactive({
  if(input$group_by2 == "rmsd") {
    hc <- hclust_rmsd2()
  }
  
  if(input$group_by2 == "rmsip") {
    hc <- hclust_rmsip2()
  }

  if(input$group_by2 == "bhat") {
    hc <- hclust_bhat2()
  }

  if(input$group_by2 == "pc_space") {
    pc <- pca1()
    hc <- hclust(dist(pc$z[,1:2]))
  }
  
  return(hc)
})

cutree2 <- reactive({
  hc <- hclust2()
  grps <- cutree(hc, k=input$nclusts)
  return(grps)
})



####################################
####     webGL functions        ####
####################################

output$struct_dropdown2 <- renderUI({
  pdbs <- align()
  ids <- 1:length(pdbs$id)
  names(ids) <-  pdbs$lab
  selectInput('viewStruct_nma', 'Show NMs for structure:',
              choices=ids)
})
  
output$nmaWebGL  <- renderWebGL({
  pdbs <- align()
  modes <- nma2()

  mag <- as.numeric(input$mag2)
  step <- mag/8
  
  trj <- mktrj(modes, pdbs=pdbs,
               mag=mag, step=step,
               s.inds=as.numeric(input$viewStruct_nma),
               m.inds=as.numeric(input$viewMode_nma),
               rock=FALSE)
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
  
    view.xyz(trj, bg.col=input$viewBGcolor2, col=col)
})



####################################
####     Plotting functions     ####
####################################

output$checkboxgroup_label_ids2 <- renderUI({
  pdbs <- align()
  grps <- cutree2()
  ids <- pdbs$lab[order(grps)]
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


make.plot.nma <- function() {
  pdbs <- align()
  modes <- nma2()

  signif <- FALSE
  if(input$cluster) {
    col <- cutree2()

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

  plot(modes, pdbs, col=col, signif=signif,
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
  rd <- rmsd2()
  hc <- hclust_rmsd2()
  grps <- cutree(hc, k=input$nclusts)
  mar <- as.numeric(c(input$margins3, input$margins3))
  heatmap(rd, distfun=as.dist, symm=TRUE,
          ColSideColors=as.character(grps),
          RowSideColors=as.character(grps),
          cexRow=input$cex3, cexCol=input$cex3,
          margins=mar)
}

output$rmsd_heatmap2 <- renderPlot({
  make.plot.heatmap_rmsd2()
})

make.plot.heatmap_rmsip2 <- function() {
  rp <- rmsip2()
  hc <- hclust_rmsip2()
  grps <- cutree(hc, k=input$nclusts)
  mar <- as.numeric(c(input$margins3, input$margins3))
  heatmap(1-rp, distfun=as.dist, symm=TRUE,
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
  pdbs <- align()
  pc <- pca1()

  hc <- hclust2()
  grps <- cutree2()
  ids <- pdbs$lab

  mar <- c(input$margins2, 5, 3, 1)
  main <- paste(toupper(input$group_by2), "dendrogram")
  hclustplot(hc, k=input$nclusts, colors=grps, labels=ids, cex=input$cex2,
             main=main, fillbox=FALSE, mar=mar)

  if(input$show_confplot2) {
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
  modes <- nma2()
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
