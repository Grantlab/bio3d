####################
##-- PCA
####################

pca1 <- reactive({
  pdbs <- fit()
  pc <- pca(pdbs)
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
  pdbs <- fit()
  grps <- clustgrps()
  ids <- basename.pdb(pdbs$id)[order(grps)]
  names(ids) <- paste(ids, " (c", grps[order(grps)], ")", sep="")
  

  checkboxInput("toggle_all", "Toggle all", TRUE)

  if(input$toggle_all) {
    checkboxGroupInput("label_ids", "PDB IDs:",
                       ids, selected=ids, inline=TRUE)
  }
  else {
    checkboxGroupInput("label_ids", "PDB IDs:",
                       ids, selected=c(), inline=TRUE)
  }
})


## normal conformer plot
  
output$pca_plot1_conf <- renderPlot({
  invisible(capture.output( pdbs <- fit() ))
  invisible(capture.output( pc <- pca1() ))
  col <- 1
  if(input$nclust>1)
    col <- clustgrps()
  
  op <- par(pty="s")
    
  xax <- as.numeric(input$pcx)
  yax <- as.numeric(input$pcy)
  
  p <- paste0("PC", c(xax, yax),
              " (", round((pc$L[c(xax, yax)]/sum(pc$L)) * 
                          100, 2), "%)")
  
  plot(pc$z[, xax], pc$z[, yax],
       col="grey50", pch=16,
       xlab=p[1], ylab=p[2],
       cex=input$cex_points*1.5)
  points(pc$z[, xax], pc$z[, yax],
         col=col, pch=16,
         cex=input$cex_points*1)
  
  abline(h = 0, col = "gray", lty = 2)
  abline(v = 0, col = "gray", lty = 2)
  if(input$labelplot) {
    if(length(input$label_ids)>0) {
      inds <- unlist(lapply(input$label_ids, grep, pdbs$id))

      if(input$distribute_labels) {
        pointLabel(pc$z[inds, xax], pc$z[inds, yax],
             labels=basename.pdb(pdbs$id[inds]),
             pos=1, offset=input$offset, cex=input$cex_labels)
      }
      else {
        text(pc$z[inds, xax], pc$z[inds, yax],
             labels=basename.pdb(pdbs$id[inds]),
             pos=1, offset=input$offset, cex=input$cex_labels)
      }
    }
  }
  invisible(par(op))
  
})
output$pca_plot1_scree <- renderPlot({
  invisible(capture.output( pdbs <- fit() ))
  invisible(capture.output( pc <- pca1() ))

  op <- par(pty="s")
  plot.pca.scree(pc$L)
  invisible(par(op))
})

## fancy plot using rChart
output$pca_plot2_conf <- renderChart2({
  invisible(capture.output( pdbs <- fit() ))
  invisible(capture.output( pc <- pca1() ))
  
  col <- 1
  if(input$nclust>1)
    col <- clustgrps()

  xax <- as.numeric(input$pcx)
  yax <- as.numeric(input$pcy)
  
  p <- paste0("PC", c(xax, yax),
              " (", round((pc$L[c(xax, yax)]/sum(pc$L)) * 
                          100, 2), "%)")
  
  ## generate a dataframe: pc$z + pdbids + group
  x <- as.data.frame(cbind(pc$z, substr(basename(pdbs$id),1,6), col))
  colnames(x)[c(xax,yax,dim(pc$z)[2]+1, dim(pc$z)[2]+2)] <- c(p,"id","group")
  
  ## use dPlot() to generate the interactive plot
  p1 <- dPlot(
    x = p[1], y = p[2], 
    groups = "id", 
    data = x, 
    type = "bubble",
    height=420,
    width=440,
    bounds = list(x=45, y=20, height=360, width=365)
  )
  
  p1$colorAxis(type = "addColorAxis", colorSeries = "group",
               palette = c("red", "blue", "black") )
  p1$xAxis(type = "addMeasureAxis")
  return(p1) 
})

output$pca_plot2_scree <- renderChart2({
  invisible(capture.output( pdbs <- fit() ))
  invisible(capture.output( pc <- pca1() ))

  PC <- c(1:length(pc$L))
  percent <- (pc$L/sum(pc$L))*100
  cumv<-cumsum(percent)
  data <- data.frame(rank=PC, var=percent, cumvar=cumv)
  #xy <- xy.coords(x, y)
  r1 <- nPlot(
    var ~ rank,
    data = data[1:20,],
    type = "lineChart"#,
    #bounds = list(x=80, y=75, height=425, width=400)
  )
  r1$xAxis(axisLabel = "Eigenvalue rank")
  r1$yAxis(axisLabel = "Proportion of variance (%)", width=50 )
  r1$chart(color = c("black","black") )
  r1$chart(forceY = c(0,100))
  r1$xAxis( tickValues=c(1:5, seq(5,20,5)) )
  r1$yAxis( tickValues=c(round(data[1:5,2],2), seq(80,100,10)) )
  r1$set(height=420,width=450)
  r1$chart(tooltipContent = "#! function(item, x, y, e){
      return 'Cumulative variance: ' + e.point.cumvar.toFixed(2) + '%'
    } !#")
  r1$chart(showLegend = FALSE)
  return(r1)
})


output$pdbs_table <- renderDataTable({
  pdbs <- fit()
  grps <- clustgrps()
  anno <- get_annotation(basename.pdb(pdbs$id))
  
  url <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=", substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  anno <- cbind(anno, url)
  anno <- cbind(anno, id=1:nrow(anno))
  anno$cluster <- grps

  return(anno[, c("id", "url", "cluster", "compound", "source", "ligandId", "chainLength")])
},  escape=FALSE)



####################################
####     Download functions     ####
####################################
traj2files <- reactive({
  dir <- data_path()
  
  pdbs <- fit()
  gaps <- gap.inspect(pdbs$ali)
  pc <- pca1()
  files <- rep(NA, 5)
  for(i in 1:5) {
    f <- paste0(dir, "/", "pc", i, ".pdb")
    trj <- mktrj(pc, pc=i, file=f,
                 resno=pdbs$resno[1, gaps$f.inds],
                 resid=pdbs$resid[1, gaps$f.inds],
                 chain=pdbs$chain[1, gaps$f.inds])
    files[i] <- f
  }
  return(files)
})

output$pctrajZIP = downloadHandler(
  filename = 'pc-traj.zip',
  content = function(file) {
    zip(file, files=traj2files(), flags = "-9Xj")
})
