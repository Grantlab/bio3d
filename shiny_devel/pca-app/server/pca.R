####################
##-- PCA
####################

pca1 <- reactive({
  pdbs <- fit()
  pc <<- pca(pdbs)
  return(pc)
})

clust <- reactive({
  pdbs <- fit()
  rd <- rmsd(pdbs)
  hc <- hclust(as.dist(rd))
  return(hc)
})

clustgrps <- reactive({
  hc <- clust()
  grps <- cutree(hc, k=input$nclust)
  return(grps)
})

output$checkboxgroup_label_ids <- renderUI({
  ids <- basename.pdb(pdbs$id)
  names(ids) <- ids
  checkboxGroupInput("label_ids", "PDB IDs:",
                     ids, selected=ids, inline=TRUE)
})


## normal conformer plot
output$pca_plot1 <- renderPlot({
  invisible(capture.output( pc <- pca1() ))
  
  col <- 1
  if(input$nclust>1)
    col <- clustgrps()
  
  if(input$screeplot)
    par(mfrow=c(1,2))
  
  xax <- as.numeric(input$pcx)
  yax <- as.numeric(input$pcy)
  
  p <- paste0("PC", c(xax, yax),
              " (", round((pc$L[c(xax, yax)]/sum(pc$L)) * 
                          100, 2), "%)")
  
  plot(pc$z[, xax], pc$z[, yax],
       col="grey50", pch=16,
       xlab=p[1], ylab=p[2],
       cex=1.6)

  points(pc$z[, xax], pc$z[, yax],
         col=col, pch=16,
         cex=1.2)
  
  abline(h = 0, col = "gray", lty = 2)
  abline(v = 0, col = "gray", lty = 2)

  if(input$labelplot) {
    if(length(input$label_ids)>0) {
      inds <- unlist(lapply(input$label_ids, grep, pdbs$id))
      text(pc$z[inds, xax], pc$z[inds, yax],
           labels=basename.pdb(pdbs$id[inds]),
           pos=1, offset=input$offset)
    }
  }
  
  if(input$screeplot) {
    plot.pca.scree(pc$L)
  }
  
})

## fancy plot using rChart
output$pca_plot2 <- renderChart2({
  invisible(capture.output( pc <- pca1() ))
  
  col <- 1
  if(input$nclust>1)
    col <- clustgrps()
  
  if(input$screeplot)
    par(mfrow=c(1,2))
  
  xax <- as.numeric(input$pcx)
  yax <- as.numeric(input$pcy)
  
  p <- paste0("PC", c(xax, yax),
              " (", round((pc$L[c(xax, yax)]/sum(pc$L)) * 
                          100, 2), "%)")
  
  ## generate a dataframe: pc$z + pdbids + group
  x <- as.data.frame(cbind(pc$z, substr(basename(pdbs$id),1,6), col))
  colnames(x)[c(xax,yax,dim(pc$z)[2]+1, dim(pc$z)[2]+2)] <- c(p,"id","group")
  
  ## use dPlot() to generate the interactive plot
  p1 <- dPlot(x = p[1], y = p[2], groups = "id", data = x, type = "bubble")
  
  p1$colorAxis(type = "addColorAxis", colorSeries = "group",
               palette = c("red", "blue", "black") )
  p1$xAxis(type = "addMeasureAxis")
  
  return(p1) 
})


output$pdbs_table <- renderDataTable({
  anno <- db_get_anno(basename.pdb(pdbs$id))
  
  url <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=", substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  anno <- cbind(anno, url)
  anno <- cbind(anno, id=1:nrow(anno))

  return(anno[, c("id", "url", "compound", "source", "ligandId", "chainLength")])
},  escape=FALSE)

