############
## NMA
###########
nma2 <- reactive({
  pdbs <- align()
    
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  
  progress$set(message = 'Calculating normal modes',
               detail = 'Please wait',
               value = 0)
  ##modes <- nma(pdbs, fit=TRUE, rm.gaps=input$rm.gaps, progress=progress)

  rm.gaps <- TRUE
  if(is.logical(input$rm.gaps))
    rm.gaps <- input$rm.gaps
    
  modes <- nma(pdbs, fit=TRUE, rm.gaps=rm.gaps, progress=progress)
  return(modes)
})

rmsip2 <- reactive({
  if(input$rm.gaps) {
    pdbs <- align()
    modes <- nma2()
    rownames(modes$rmsip) <- basename.pdb(rownames(modes$rmsip))
    colnames(modes$rmsip) <- basename.pdb(colnames(modes$rmsip))
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
  rownames(rd) <- basename.pdb(pdbs$id)
  colnames(rd) <- basename.pdb(pdbs$id)
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
    rownames(bc) <- basename.pdb(pdbs$id)
    colnames(bc) <- basename.pdb(pdbs$id)
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
  if(input$group_by == "rmsd") {
    hc <- hclust_rmsd2()
  }
  
  if(input$group_by == "rmsip") {
    hc <- hclust_rmsip2()
  }

  if(input$group_by == "bhat") {
    hc <- hclust_bhat2()
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

output$struct_dropdown <- renderUI({
  pdbs <- align()
  ids <- 1:length(pdbs$id)
  names(ids) <-  basename.pdb(pdbs$id)
  selectInput('viewStruct', 'Choose Structure:',
              choices=ids)
})
  
output$nmaWebGL  <- renderWebGL({
  pdbs <- align()
  modes <- nma2()
  trj <- mktrj(modes, pdbs=pdbs,
               s.inds=as.numeric(input$viewStruct),
               m.inds=as.numeric(input$viewMode),
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
  
  typ <- switch(input$viewColor2,
                'mag' = 2,
                'amalgam' = 1,
                'default' = 1
                )
  
    view.xyz(trj, bg.col=input$viewBGcolor2, col=col, add=TRUE, type=typ)
})



####################################
####     Plotting functions     ####
####################################

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

  plot(modes, pdbs, col=col, signif=signif,
       spread=input$spread, conservation=input$seqide)

}


output$nma_plot <- renderPlot({
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
  ids <- basename.pdb(pdbs$id)

  mar <- c(input$margins2, 5, 3, 1)
  main <- paste(toupper(input$group_by), "dendrogram")
  hclustplot(hc, k=input$nclusts, colors=grps, labels=ids, cex=input$cex2,
             main=main, fillbox=FALSE, mar=mar)

  par(fig=c(.65, 1, .45, 1), new = TRUE)
  plot(pc$z[,1:2], col="grey50", pch=16, cex=1.1*input$cex2, 
       ylab="", xlab="", axes=FALSE)
  points(pc$z[,1:2], col=grps, pch=16, cex=0.7*input$cex2)
  box()
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

trj2pdb2  <- reactive({
    dir <- data_path()
    pdbs <- align()
    modes <- nma2()
    gaps <- gap.inspect(pdbs$ali)
    fname  <- paste0(dir, '/', 'pc', input$viewPC, '.pdb')
    trj <- mktrj(modes, pdbs=pdbs,
                 s.inds=as.numeric(input$viewStruct),
                 m.inds=as.numeric(input$viewMode),
                 file=fname)
                 #resno=pdbs$resno[1, gaps$f.inds],
                 #resid=pdbs$resid[1, gaps$f.inds],
                 #chain=pdbs$chain[1, gaps$f.inds])
    
    return(fname)
})

output$nmtraj = downloadHandler(
    filename=function() {
        paste0('pc', input$viewMode, '.pdb')
    },
    content=function(file) {
        # Avoid possibility of not having write permission on server
        src  <- normalizePath('nma-traj.pdb')
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        file.copy(src, 'nma-traj.pdb')
        file.rename(trj2pdb2(), file)
    }

)
