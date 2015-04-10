#####################
##-- Superimpose and RMSD stuff
####################

find_core <- reactive({
  pdbs <- align()
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Finding invariant core',
               detail = 'This may take some time...')
  progress$set(value = 2)
  
  core <<- core.find(pdbs)
  progress$set(value = 5)
  return(core)
})

fit <- reactive({
  if(input$fit_type == "full") {
    pdbs$xyz <<- pdbfit(pdbs)
  }
  else {
    core <- find_core()
    pdbs$xyz <<- pdbfit(pdbs, core)
  }
  
  return(pdbs)
})

rmsd1 <- reactive({
  pdbs <- fit()
  rd <- rmsd(pdbs)
})

####################################
####     Plotting functions     ####
####################################
output$rmsd_heatmap <- renderPlot({
  pdbs <- fit()
  rd <- rmsd1()

  rownames(rd) <- basename.pdb(pdbs$id)
  colnames(rd) <- basename.pdb(pdbs$id)

  hc <- hclust(as.dist(rd))
  grps <- cutree(hc, k=input$clusters)
  
  heatmap(rd, distfun=as.dist, symm=TRUE,
          ColSideColors=as.character(grps),
          RowSideColors=as.character(grps)
          )
})

output$rmsd_hist <- renderPlot({
  pdbs <- fit()
  rd <- rmsd1()
  hist(rd[upper.tri(rd)], breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")
})

output$rmsd_dendrogram <- renderPlot({
  pdbs <- fit()
  rd <- rmsd1()
  pdbs$id <- basename.pdb(pdbs$id)
  hc <- hclust(as.dist(rd))
  hclustplot(hc, k=input$clusters, labels=pdbs$id, cex=input$cex,
             ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)

})

output$rmsf_plot <- renderPlot({
  pdbs <- fit()
  rd <- rmsd1()
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  resno <- pdbs$resno[1, gaps.res$f.inds]
  ##sse <- .pdbs2sse(pdbs, ind=1, rm.gaps=TRUE)
  
  rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])
  plot.bio3d(rf, resno=resno, 
             ylab="RMSF (Å)", xlab="Residue No.")
})

####################################
####     Data tables     ####
####################################

output$print_core <- renderDataTable({
  core <- find_core()
  gaps <- gap.inspect(pdbs$ali)

  ##ind <- grep(input$reference_id, pdbs$id)
  ind <- 1
  resid <- pdbs$resid[ind, core$atom]
  resno <- pdbs$resno[ind, core$atom]
  bds <- bounds(resno)

  resid.inds <- which(resno %in% bds[,"start"])
  start <- paste0(resid[resid.inds], bds[, "start"])

  resid.inds <- which(resno %in% bds[,"end"])
  end <- paste0(resid[resid.inds], bds[, "end"])
  
  out <- data.frame(start, end, length=bds[, "length"])
  return(out)
  
}, options = list(searching=FALSE, lengthChange=FALSE, paging=FALSE))

output$reference_selector <- renderUI({
  ids <- basename.pdb(pdbs$id)
  names(ids) <- ids
  selectInput("reference_id", "Reference PDB id:", ids)
})

output$rmsd_table <- renderDataTable({
  if(!is.null(input$reference_id))
    ind <- grep(input$reference_id, pdbs$id)
  else
    ind <- 1
  
  gaps <- gap.inspect(pdbs$xyz)
  rd <- rmsd(pdbs$xyz[ind, gaps$f.inds],
             pdbs$xyz[, gaps$f.inds],
             fit = TRUE)
  names(rd) <- basename.pdb(pdbs$id)
  
  return(data.frame(ids=names(rd), rmsd=rd))
}, options = list(searching=FALSE, lengthChange=FALSE, paging=TRUE))



####################################
####     Download functions     ####
####################################
pdbs2files <- reactive({
  core <- find_core()
  pdbs$xyz <<- pdbfit(pdbs, core)
  pdbfit(pdbs, core, outpath="fitlsq")
  return(paste0("fitlsq/", basename(pdbs$id), "_flsq.pdb"))
})

output$pdbsZIP = downloadHandler(
  filename = 'pdbs.zip',
  content = function(file) {
    zip(file, files=pdbs2files())
  })
