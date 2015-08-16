init_show_trj_nma <- TRUE

########################
## NMA TAB observers  ##
########################

observeEvent(input$run_nma, {
  rv$pdbs4nma <- pdbs4nma()
  rv$modes <- nma2()
})

observeEvent(input$filter_rmsd, {
  if(nma_allowed()$value)
    rv$nma_allowed2 <- TRUE
  else
    rv$nma_allowed2 <- FALSE
})

observeEvent(input$filter_clusterBy, {

  if(!is.null(input$filter_rmsd)) {
    
    cut <- 1
    lab <- "RMSD Cutoff"
    step <- 0.25

    if(input$filter_clusterBy == "pca") {
      #cut <- cutree_pca2()
      cut <- 1
      lab <- "PC Cutoff"
      step <- 0.5
    }
    
    if(input$filter_clusterBy == "seqide") {
      #cut <- cutree_seqide()
      cut <- 0.1
      lab <- "Identity Cutoff"
      step <- 0.1
    }
    
    updateNumericInput(session, "filter_rmsd", lab, value = cut, step = step)
  }
})

## rv for output
output$nma_allowed2 <- reactive({
  nma_allowed()$value
})
outputOptions(output, 'nma_allowed2', suspendWhenHidden=FALSE)

output$nmaIsDone <- reactive({
  return(!is.null(rv$modes))
})
outputOptions(output, 'nmaIsDone', suspendWhenHidden=FALSE)


############
## NMA
###########
nma2 <- reactive({

  if(!nma_allowed()$value)
    return(NULL)
  
  pdbs <- rv$pdbs4nma
    
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

nma_allowed <- reactive({
  if(!rv$aligned | !rv$fitted) {
    return(list(value=FALSE, msg="<strong style='color: red;'>Unfinshed required calculations</strong>. Please go back to the ALIGN and/or FIT tabs. "))
  }
  
  pdbs <- align()
  inds <- filter_pdbs()
  pdbs <- trim.pdbs(pdbs, row.inds=inds)
  gaps <- gap.inspect(pdbs$ali)

  lens <- ncol(pdbs$ali) - gaps$row
  size <- sum(lens)

  if(!length(inds)>1){
    return(list(value=FALSE, msg="<strong style='color: red;'>> 1 structure needed for esemble NMA</strong>. Decrease cutoff value for filtering. "))
  }
  
  if(size > 5000)
    return(list(value=FALSE, msg="<strong style='color: red;'>Ensemble too large for NMA</strong>. Increase cutoff value for filtering. "))
  else
    return(list(value=TRUE, msg=""))
})

## returns indices for filtering
filter_pdbs <- reactive({
  pdbs <- align()
  rd <- rmsd1()

  cut <- as.numeric(input$filter_rmsd)
  if(length(cut)>0) {
    if(is.null(cut) | is.na(cut)) {
      cut <- 1
    }
  }
  else {
    cut <- 1
  }

  if(input$filter_clusterBy == "rmsd") {
    message("filter by rmsd")
    inds <- filter.rmsd(pdbs$xyz, rmsd.mat = rd,
                        cutoff=cut, fit=FALSE,
                        method = input$hclustMethod_rmsd)$ind
  }

  if(input$filter_clusterBy == "pca") {
    message("filter by pca")
    rd <- pc_dist()
    inds <- filter.mat(as.matrix(rd), dist.fun=as.dist,
                       cutoff=cut, 
                       method = input$hclustMethod_pca)$ind
  }

  if(input$filter_clusterBy == "seqide") {
    message("filter by seqide")
    ide <- seqide()
    rownames(ide) <- pdbs$lab
    colnames(ide) <- pdbs$lab
    
    inds <- filter.mat(1-ide, dist.fun=as.dist, 
                       cutoff=cut, 
                       method = input$hclustMethod)$ind
  }

  message("inds by function filter_pdbs():")
  message(paste(inds, sep=" "))
  message(paste(pdbs$lab[inds], sep=" "))
  
  return(inds)
})

## returns a new pdbs object for NMA
pdbs4nma <- reactive({
  message("pdbs4nma()")
  
  pdbs1 <- align()
  inds <- filter_pdbs()
  pdbs2 <- trim.pdbs(pdbs1, row.inds=inds)
  pdbs2$lab <- pdbs1$lab[inds]

  message(paste(pdbs2$lab, sep=" "))
  return(pdbs2)
})

## dendrogram for the top row (filtering structures)
output$rmsd_dendrogram2 <- renderPlot({
  pdbs <- fit()
  rd <- rmsd1()
  rownames(rd) <- pdbs$lab
  colnames(rd) <- pdbs$lab

  hc <- hclust_rmsd()
  inds <- filter_pdbs()
  
  col <- rep(2, length(hc$order))
  col[inds] <- 1

  mar <- c(input$margins, 5, 3, 1)
  hclustplot(hc, k=1, col=col, main="Cluster Dendrogram")
  
})


output$filter_summary <- renderUI({
  invisible(capture.output(pdbs1 <- align()))
  invisible(capture.output(inds <- filter_pdbs()))

  str <- paste("<strong>Structures included</strong>:", length(inds), "/", length(pdbs1$id))

  if(!nma_allowed()$value)
    str <- c(str, nma_allowed()$msg)

  HTML(paste(str, collapse="<br>"))
})


output$filter_rmsdInput <- renderUI({
  #cut <- cutree_rmsd()
  lab <- "RMSD Cutoff"
  numericInput("filter_rmsd", lab, value = 1, step = 0.25)
})





####################################
####   similarity measures     ####
####################################


## RMSIP
rmsip2 <- reactive({
  if(input$rm_gaps) {
    pdbs <- rv$pdbs4nma
    modes <- rv$modes
    rownames(modes$rmsip) <- pdbs$lab
    colnames(modes$rmsip) <- pdbs$lab
    return(modes$rmsip)
  }
  else {
    stop("RMSIP only available with 'Omit gaps'")
  }
})


## Bhattacharyya coefficient
bhat <- reactive({
  if(input$rm_gaps) {
    pdbs <- rv$pdbs4nma
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
  message("hclust_rmsip")
  
  rd <- rmsip2()
  return(hclust(as.dist(1-rd),
         method=input$hclustMethod_rmsip))
})

hclust_bhat2 <- reactive({
  rd <- bhat2()
  return(hclust(as.dist(1-rd)))
})

hclust2 <- reactive({
  message("hclust2")
  
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
  message("cutree_rmsip")
  hc <- hclust_rmsip()
  
  if(is.null(input$splitTreeK_rmsip))
    k <- NA
  else
    k <- as.numeric(input$splitTreeK_rmsip)
  
  cut <- cutreeBio3d(hc, minDistance=input$minDistance_rmsip, k=k)
  return(cut)
})


cutree_nma <- reactive({
  message("cutree_nma")
  inds <- filter_pdbs()
    
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

  print(cut$grps)
  
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
  pdbs <- rv$pdbs4nma
  ids <- 1:length(pdbs$id)
  names(ids) <-  pdbs$lab
  selectInput('viewStruct_nma', 'Show NMs for structure:',
              choices=ids)
})


## this is to avoid calling RGL functions twice
## needed after hiding control panel before the the NMA calc is done
rv$viewStruct_nma <- 1
rv$viewMode_nma <- 1
rv$mag2 <- 5
rv$viewColor2 <- "amalgam"
rv$viewBGColor2 <- "white"

observeEvent(input$viewStruct_nma, {
  rv$viewStruct_nma <- as.numeric(input$viewStruct_nma)
})

observeEvent(input$viewMode_nma, {
  rv$viewMode_nma <- as.numeric(input$viewMode_nma)
})

observeEvent(input$viewColor2, {
  rv$viewColor2 <- input$viewColor2
})

observeEvent(input$viewBGColor2, {
  rv$viewBGColor2 <- input$viewBGColor2
})

observeEvent(input$mag2, {
  rv$mag2 <- as.numeric(input$mag2)
})
  
output$nmaWebGL  <- renderWebGL({
  message("nmaWebGL called")
  pdbs <- rv$pdbs4nma
  modes <- rv$modes

  if(!inherits(modes, "enma") | !inherits(pdbs, "pdbs"))
    return(NULL)
  
  mag <- as.numeric(rv$mag2)
  step <- mag/8
  
  trj <- mktrj(modes, pdbs=pdbs,
               mag=mag, step=step,
               s.inds=rv$viewStruct_nma,
               m.inds=rv$viewMode_nma,
               rock=FALSE,
               file = paste0(data_path(), '/mode_', rv$viewMode_nma, '.pdb'))
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
  col <- switch(rv$viewColor2,
                'mag' = magcol(), # vec2color(rmsf(m)), #!! col=col, type=2
                'amalgam' = amalcol(1:n),
                'default' = colorRampPalette(c('blue', 'gray', 'red'))(n)
                )

  view.xyz(trj, bg.col=rv$viewBGColor2, col=col, d.cut=6)
})



####################################
####     Plotting functions     ####
####################################

## checkbox labels for filtering structures in fluctuation plot
output$checkboxgroup_label_ids2 <- renderUI({
  pdbs <- rv$pdbs4nma
  grps <- cutree_nma()$grps

  message(grps)
  
  ids <- pdbs$lab
  message(ids)
  
  ##names(ids) <- paste(ids, " (c", grps[order(grps)], ")", sep="")
  names(ids) <- paste(ids, " (c", grps, ")", sep="")
  ids <- ids[ order(grps) ]
  print(ids)

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
  pdbs <- rv$pdbs4nma
  modes <- rv$modes
  
  gaps.res <- gap.inspect(pdbs$ali)
  ##gaps.pos <- gap.inspect(pdbs$xyz)
  
  resno <- pdbs$resno[1, gaps.res$f.inds]
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

  par(mar=c(4, 5, 4, 2))
  plot(modes, pdbs, sse=sse, col=col, signif=signif,
       spread=input$spread, conservation=input$seqide, resno=resno,
       main = "NMA derived fluctuations")

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
  inds <- filter_pdbs()
  
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
  pdbs <- rv$pdbs4nma
  pc <- pca1()

  hc <- hclust_rmsip()
  cut <- cutree_rmsip()
  grps <- cut$grps
  k <- cut$k
  ids <- pdbs$lab

  mar <- c(input$margins2, 5, 3, 1)
  main <- paste(toupper(input$group_by2), "RMSIP Cluster Dendrogram")
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
  pdbs <- rv$pdbs4nma
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

