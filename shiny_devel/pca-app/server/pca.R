####################
##-- PCA
####################
init_show_trj <- TRUE

pca1 <- reactive({
  pdbs <- fit()
  pc <- pca(pdbs)

  if(init_show_trj) {
    updateCheckboxInput(session, 'show_trj', 'Show PC Trajectory', value=TRUE)
    init_show_trj <<- FALSE
  }

  return(pc)
})

clust <- reactive({
  pdbs <- fit()
  pc <- pca1()

  if(input$cluster_by == "rmsd") {
    rd <- rmsd(pdbs)
    hc <- hclust(as.dist(rd))
  }

  if(input$cluster_by == "pc_space") {
    hc <- hclust(dist(pc$z[,1:as.numeric(input$clust_npcs)]))
  }
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
  ids <- pdbs$lab[order(grps)]
  names(ids) <- paste(ids, " (c", grps[order(grps)], ")", sep="")

  checkboxInput("toggle_all", "Toggle all", TRUE)

  if(input$toggle_all) {
    checkboxGroupInput(inputId="label_ids", label="Label PDB IDs:",
                       choices=ids, selected=ids, inline=TRUE)

    ##lapply(1:input$nclust, function(x) {
    ##do.call(checkboxGroupInput,list("label_ids1", paste0("Cluster: ", x), ids[grps[order(grps)]==x], selected=ids, inline=TRUE))
    ##})
  }
  else {
        checkboxGroupInput(inputId="label_ids", label="PDB IDs:",
                           choices=ids, selected=c(), inline=TRUE)
  }
})


## normal conformer plot

output$pca_plot1_conf <- renderPlot({
  invisible(capture.output( pdbs <- fit() ))
  invisible(capture.output( pc <- pca1() ))
  col <- 1
  if(as.numeric(input$nclust) > 1)
    col <- clustgrps()

  op <- par(pty="s")

  xax <- as.numeric(input$pcx)
  yax <- as.numeric(input$pcy)

  p <- paste0("PC", c(xax, yax),
              " (", round((pc$L[c(xax, yax)]/sum(pc$L)) *
                          100, 2), "%)")
  if(input$nclust < 9) {
    cluster.colors <- col2hex(palette())[col]
  }
  else {
    cluster.colors <- col2hex(rainbow(input$nclust))[col]
  }

  xlim <- range(pc$z[, xax]) * as.numeric(input$inner_margin)
  ylim <- range(pc$z[, yax]) * as.numeric(input$inner_margin)

  plot(pc$z[, xax], pc$z[, yax], xlab=p[1], ylab=p[2],
       bg=cluster.colors, pch=21,
       cex=as.numeric(input$cex_points), col='grey50',
       xlim=xlim, ylim=ylim)

  abline(h = 0, col = "gray", lty = 2)
  abline(v = 0, col = "gray", lty = 2)

  if(input$labelplot) {
    if(length(input$label_ids)>0) {
      inds <- unlist(lapply(input$label_ids, grep, pdbs$lab))

      if(input$distribute_labels) {
        pointLabel(pc$z[inds, xax], pc$z[inds, yax],
                   labels=pdbs$lab[inds],
                   pos=1, offset=as.numeric(input$offset),
                   cex=as.numeric(input$cex_labels))
      }
      else {
        text(pc$z[inds, xax], pc$z[inds, yax],
             labels=pdbs$lab[inds],
             pos=1, offset=input$offset, cex=as.numeric(input$cex_labels))
      }
    }
  }

  if(!is.null(input$pdbs_table_rows_selected)) {
    inds <- as.numeric(input$pdbs_table_rows_selected)
    if(length(inds)>0)
      points(pc$z[inds, xax], pc$z[inds, yax], col=2, cex=2.5*as.numeric(input$cex_points))
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

output$scatterplot3d_webgl <- renderWebGL({
  invisible(capture.output( pdbs <- fit() ))
  invisible(capture.output( pc <- pca1() ))
  col <- 1
  if(input$nclust>1)
    col <- clustgrps()

  op <- par(pty="s")

  xax <- as.numeric(input$pcx)
  yax <- as.numeric(input$pcy)
  zax <- as.numeric(input$pcz)

  p <- paste0("PC", c(xax, yax, zax),
              " (", round((pc$L[c(xax, yax, zax)]/sum(pc$L)) *
                          100, 2), "%)")
  if(input$nclust < 9) {
    cluster.colors <- col2hex(palette())[col]
  }
  else {
    cluster.colors <- col2hex(rainbow(input$nclust))[col]
  }

  points3d(pc$z[, xax], pc$z[, yax], pc$z[, zax], col=cluster.colors, size=5)
  axes3d()
  title3d('','', p[1], p[2], p[3])
})

output$scatterplot3d_rthreejs <- renderScatterplotThree({
  invisible(capture.output( pdbs <- fit() ))
  invisible(capture.output( pc <- pca1() ))

  col <- 1
  if(input$nclust>1)
    col <- clustgrps()

  op <- par(pty="s")

  xax <- as.numeric(input$pcx)
  yax <- as.numeric(input$pcy)
  zax <- as.numeric(input$pcz)

  p <- paste0("PC", c(xax, yax, zax),
              " (", round((pc$L[c(xax, yax, zax)]/sum(pc$L)) *
                          100, 2), "%)")
  if(input$nclust < 9) {
    cluster.colors <- col2hex(palette())[col]
  }
  else {
    cluster.colors <- col2hex(rainbow(input$nclust))[col]
  }

  labs <- pdbs$lab

  names(labs) <- NULL
  df <- data.frame(x=pc$z[,xax], y=pc$z[,yax], z=pc$z[,zax])
  names(df) <- p
  scatterplot3js(x=df,
                 color=cluster.colors,
                 size=1,
                 labels=labs,
                 label.margin="80px 10px 10px 10px",
                 renderer=input$renderer,
                 grid=input$grid)
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
  x <- as.data.frame(cbind(pc$z, pdbs$lab, col))
  colnames(x)[c(xax,yax,dim(pc$z)[2]+1, dim(pc$z)[2]+2)] <- c(p,"id","cluster")

  ## use rainbow palette for more than 8 clusters
  if(as.numeric(input$nclust) < 9) {
      cluster.colors <- col2hex(palette())
  } else {
      cluster.colors <- col2hex(rainbow(input$nclust))
  }

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

  p1$colorAxis(type = "addColorAxis", colorSeries = "cluster",
               palette = if(input$nclust==1) rep('black',3) else cluster.colors[1:input$nclust] )
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

make.plot.loadings <- function(){
  pdbs <- fit()
  pc <- pca1()
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  resno <- pdbs$resno[1, gaps.res$f.inds]
  sse <- pdbs2sse(pdbs, ind=1, rm.gaps=TRUE)
  rf <- rmsf(pdbs$xyz[, gaps.pos$f.inds])

  pcs <- as.numeric(input$loadings_pc)

  if(length(pcs)>1)
    par(mfrow=c(length(pcs), 1))
  if(length(pcs)==0)
    pcs <- 1

  if(input$toggle_rmsf1)
    par(mar=c(5, 4, 4, 5))
  else
    par(mar=c(5, 4, 4, 2))

  for(i in pcs) {
    plot4 <- plot.bio3d(pc$au[, i], resno=resno, sse=sse,
                        ylab=paste0("PC-", i, " (Å)"), xlab="Residue No.")

    if(input$toggle_rmsf1) {
      par(new=TRUE)
      plot5 <- plot.bio3d(rf, resno=resno, sse=sse, axes=FALSE, col=2, type="l", xlab="", ylab="")
      axis(4, col=2)
      mtext("RMSF (Å)", side=4, line=3)
    }
  }

  return(plot4)
}

output$loadings_plot <- renderPlot({
  print(make.plot.loadings())
})

output$pcaWebGL  <- renderWebGL({
    #ptm <- proc.time()
    #pdbs <- fit()
    pc <- pca1()
    
    mag <- as.numeric(input$mag1)
    step <- mag/8
    
    trj <- mktrj(pc, pc=as.numeric(input$viewPC), rock=FALSE,
                 mag=mag, step=step)
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
    col <- switch(input$viewColor,
                  'mag' = magcol(), # vec2color(rmsf(m)), #!! col=col, type=2
                  'amalgam' = amalcol(1:n),
                  'default' = colorRampPalette(c('blue', 'gray', 'red'))(n)
                 )

    view.xyz(trj, bg.col=input$viewBGcolor, col=col, d.cut=6)
    #proc.time()
    #cat(proc.time() - ptm)
})

observeEvent(input$viewUpdate, {
    updateSelectInput(session, 'viewPC', label='Choose Principal Component:',
        choices=c(1:10))
    updateRadioButtons(session, 'viewColor', label='Structure color',
        choices=list(
          'Amalgam' = 'amalgam',
          'Magnitude'='mag',
          'By Frame (blue->gray->red)'='default'
          ),
        selected='amalgam')
    updateRadioButtons(session, 'viewBGcolor', label='Background color',
        choices=list('Black'='black', 'White'='white'),
        selected='white')
})

get_pdbstable <- reactive({
  pdbs <- fit()
  grps <- clustgrps()
  anno <- get_annotation(pdbs$lab)

  pdbId <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=", substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  anno <- cbind(anno, pdbId)
  anno <- cbind(anno, id=1:nrow(anno))
  if(input$nclust < 9) {
      cluster.colors <- col2hex(palette())
  } else {
      cluster.colors <- col2hex(rainbow(input$nclust))
  }
  anno$cluster <- paste0("<span style=\"color:",
    cluster.colors[grps], "; font-size:large\">&#x25CF;</span>&nbsp;",
    grps)
  
  rownames(anno) <- NULL
  show.cols <- c("acc", "cluster", "compound", "source", "ligandId")
  col.inds <- sapply(show.cols, grep, colnames(anno))
  
  return(anno[, col.inds])
})

output$pdbs_table <- renderDataTable({
  datatable(get_pdbstable(), extensions = 'Scroller', escape = FALSE,
            colnames = c("ID", "Cluster", "Name", "Species", "Ligands"),
            selection = "none",
            options = list(
              deferRender = TRUE,
              dom = "frtiS",
              scrollY = 200,
              scrollCollapse = TRUE,
              autoWidth = FALSE,
              columnDefs = list(list(width = '40%', targets = c(list(3))))
              ))
})



####################################
####     Download functions     ####
####################################
#traj2files <- reactive({
#  dir <- data_path()
#
#  pdbs <- fit()
#  gaps <- gap.inspect(pdbs$ali)
#  pc <- pca1()
#  files <- rep(NA, 5)
#  for(i in 1:5) {
#    f <- paste0(dir, "/", "pc", i, ".pdb")
#    trj <- mktrj(pc, pc=i, file=f,
#                 resno=pdbs$resno[1, gaps$f.inds],
#                 resid=pdbs$resid[1, gaps$f.inds],
#                 chain=pdbs$chain[1, gaps$f.inds])
#    files[i] <- f
#  }
#  return(files)
#})
#
#output$pctrajZIP = downloadHandler(
#  filename = 'pc-traj.zip',
#  content = function(file) {
#    zip(file, files=traj2files(), flags = "-9Xj")
#})

pca2pdb  <- reactive({
    path <- data_path()
    pdbs <- fit()
    gaps <- gap.inspect(pdbs$ali)
    pc <- pca1()
    fname  <- paste0(path, '/', 'pc', input$viewPC, '.pdb')
    mktrj.pca(pca=pc, pc=as.numeric(input$viewPC), file=fname,
              resno=pdbs$resno[1, gaps$f.inds],
              resid=pdbs$resid[1, gaps$f.inds],
              chain=pdbs$chain[1, gaps$f.inds])
    return(fname)
})

output$pctraj = downloadHandler(
  filename=function() {
    paste0('pc', input$viewPC, '.pdb..zip')
  },
  content=function(file) {
    ## Avoid possibility of not having write permission on server
    #src  <- normalizePath('pc-traj.pdb')
    #owd <- setwd(tempdir())
    #on.exit(setwd(owd))
    #file.copy(src, 'pc-traj.pdb')
    #file.rename(trj2pdb(), file)
    zip(file, files=pca2pdb(), flags = "-9Xj")
    }

  )

make_pca_pse <- reactive({
  path <- data_path()
  pc <- pca1()
  outf <- paste0(path, "/pc", as.numeric(input$viewPC), ".pse")
  file <- pymol.modes(pc, mode=as.numeric(input$viewPC), type="session",
                file=outf)
  return(outf)
})

output$pca2pymol = downloadHandler(
  filename=function() {
    paste0('pc', input$viewPC, '.pse.zip')
  },
  content = function(file) {
    zip(file, files=make_pca_pse(), flags = "-9Xj")
})


output$pcloadings2pdf = downloadHandler(
  filename = "pc_loadings.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, width=input$width2, height=input$height2)
    make.plot.loadings()
    dev.off()
})
