#####################
##-- Superimpose and RMSD stuff
####################
init_show_pdbs <- TRUE
find_core <- reactive({
  pdbs <- align()

  progress <- shiny::Progress$new()
  on.exit(progress$close())

  progress$set(message = 'Finding core',
               detail = 'Please wait ...',
               value = 0)

  core <- core.find(pdbs, progress=progress)
  progress$close()
  return(core)
})

fit <- reactive({
  pdbs <- align()
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
  return(pdbs)
})

rmsd1 <- reactive({
  pdbs <- fit()
  rd <- rmsd(pdbs)

  return(rd)
})

hclust1 <- reactive({
  rd <- rmsd1()
  hc <- hclust(as.dist(rd))
})

cutree1 <- reactive({
  hc <- hclust1()
  grps <- cutree(hc, k=input$clusters)
})


representatives <- reactive({
  pdbs <- fit()
  rd <- rmsd1()
  hc <- hclust1()
  grps <- cutree1()

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

  return(pdbs$id[rep.inds])
})

####################################
####     Plotting functions     ####
####################################
make.plot.heatmap <- function() {
  pdbs <- fit()
  rd <- rmsd1()
  rownames(rd) <- basename.pdb(pdbs$id)
  colnames(rd) <- basename.pdb(pdbs$id)
  hc <- hclust(as.dist(rd))
  grps <- cutree(hc, k=input$clusters)
  mar <- as.numeric(c(input$margins, input$margins))
  plot1 <- heatmap(rd, distfun=as.dist, symm=TRUE,
          ColSideColors=as.character(grps),
          RowSideColors=as.character(grps),
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
  pdbs$id <- basename.pdb(pdbs$id)
  hc <- hclust(as.dist(rd))
  mar <- c(input$margins, 5, 3, 1)
  plot3 <- hclustplot(hc, k=input$clusters, labels=pdbs$id, cex=input$cex,
             ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE,
             mar = mar)
  return(plot3)
}

output$rmsd_dendrogram <- renderPlot({
  print(make.plot.rmsd.dendogram())
})

make.plot.rmsf <- function(){
  pdbs <- fit()
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  resno <- pdbs$resno[1, gaps.res$f.inds]
  sse <- pdbs2sse(pdbs, ind=1, rm.gaps=TRUE)

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

output$representatives <- renderPrint({
  pdbs <- fit()
  reps <- basename.pdb(representatives())

  cat(reps, sep="\n")

})

output$print_core <- renderDataTable({
  pdbs <- fit()
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
  pdbs <- fit()
  ids <- basename.pdb(pdbs$id)
  names(ids) <- ids
  selectInput("reference_id", "Reference PDB id:", ids)
})

output$rmsd_table <- renderDataTable({
  pdbs <- fit()
  rd <- rmsd1()
  hc <- hclust1()
  grps <- cutree1()

  if(!is.null(input$reference_id))
    ind <- grep(input$reference_id, pdbs$id)
  else
    ind <- 1

  gaps <- gap.inspect(pdbs$xyz)
  rd <- rmsd(pdbs$xyz[ind, gaps$f.inds],
             pdbs$xyz[, gaps$f.inds],
             fit = TRUE)
  names(rd) <- basename.pdb(pdbs$id)

  return(data.frame(ids=names(rd), rmsd=round(rd,1), "cluster"=grps))
}, options = list(searching=FALSE, lengthChange=FALSE, paging=TRUE))


output$pdbsWebGL  <- renderWebGL({
  pdbs <- fit()
  xyz <- pdbs$xyz
  n <- nrow(xyz)
  grps <- cutree1()
  core <- find_core()
  dims <- dim(pdbs$ali)
  gaps <- gap.inspect(pdbs$ali)

  corecol <- function() {
    col <- matrix(1, nrow=dims[1], ncol=dims[2])
    col[, core$atom] <- 2
    return(col)
  }

  gapscol <- function() {
    col <- matrix(1, nrow=dims[1], ncol=dims[2])
    if(length(gaps$t.inds)>0)
      col[, gaps$t.inds] <- 2
    return(col)
  }

  col <- switch(input$viewColor1,
                'struct' = vec2color(1:n),
                'cluster' = grps,
                'core' = corecol(),
                'gaps' = gapscol()
                )

  typ <- switch(input$viewColor1,
                'struct' = 1,
                'cluster' = 1,
                'core' = 2,
                'gaps' = 2)

  view.xyz(xyz, bg.col=input$viewBGcolor1, col=col, type=typ, add=TRUE)
})

observeEvent(input$viewUpdate1, {
  updateRadioButtons(session, 'viewColor1', label='Structure color',
                     choices=list(
                       'By cluster ID'='cluster',
                       'By structure ID'='struct',
                       'Invariant core'='core',
                       'Gap regions'='gaps'
                       ),
                     selected='cluster')

  updateRadioButtons(session, 'viewBGcolor1', label='Background color',
                     choices=list('Black'='black', 'White'='white'),
                     selected='white')
})



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
  files <- paste0(path, "/", basename.pdb(pdbs$id), ".pdb_flsq.pdb")
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
  rownames(rd) <- basename.pdb(pdbs$id)
  colnames(rd) <- basename.pdb(pdbs$id)
  file <- paste0(path, "/", "rmsd-mat.dat")
  write.table(rd, file=file, quote=FALSE)
  return(file)
})

output$rmsdZIP = downloadHandler(
  filename = 'rmsd-mat.zip',
  content = function(file) {
    zip(file, files=rmsd2files())
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
  
  col <- switch(input$viewColor1,
                "cluster" = cutree1(),
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
