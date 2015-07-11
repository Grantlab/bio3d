
############################
## Overlap analysis
###########################

observeEvent(input$run_overlap, {
  rv$overlap <- calc_overlap()
})

calc_overlap <- reactive({
  message("overlapping")
  hits <- rv$blast

  if(length(input$blast_table_rows_selected)>0) {
    ids <- hits[as.numeric(input$blast_table_rows_selected)]
  }
  else {
    return(NULL)
  }

  message(ids)

  pdblist <- vector("list", length(ids)+1)
  pdblist[[1]] <- trim(rv$final_pdb, "calpha")

  for(i in 1:length(ids)) {
    pdblist[[i+1]] <- read_pdb(ids[i])
    pdblist[[i+1]] <- trim(pdblist[[i+1]], "calpha", chain=substr(ids[i], 6, 6))
  }

  ## pdb for nma
  pdb <- pdblist[[1]]

  pdbs <- pdbaln(pdblist, exefile=configuration$muscle$exefile,
                 outfile=tempfile(pattern="aln", fileext=".fasta"))
  pdbs$id[1] <- input$pdbid
  pdbs$id[2:length(pdbs$id)] <- ids
  pdbs$xyz <- pdbfit(pdbs)

  gaps.pos <- gap.inspect(pdbs$xyz)
  gaps.res <- gap.inspect(pdbs$ali)

  ## indices for nma calc
  tomatch <- paste(pdbs$resno[1, gaps.res$f.inds],
                   pdbs$chain[1, gaps.res$f.inds], sep="|")
  resids <-  paste(pdb$atom$resno, pdb$atom$chain, sep="|")
  nmainds <- as.select(resids %in% tomatch)


  ## re-calculate normal modes
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Calculating normal modes',
               detail = 'Please wait')
  progress$set(value = 3)
  
  if(input$forcefield %in% c("calpha", "sdenm", "reach"))
    modes <- nma(pdb, ff=rv$forcefield, mass=FALSE, temp=300,
                 outmodes = nmainds)
  
  if(input$forcefield %in% c("anm", "pfanm"))
    modes <- nma(pdb, ff=rv$forcefield, cutoff=rv$cutoff,
                 mass=FALSE, temp=NULL, outmodes = nmainds)

  progress$set(value = 5)
  progress$close()
  
  
  ## calculate overlap
  progress <- shiny::Progress$new(session, min=1, max=length(pdbs$id))
  on.exit(progress$close())
  progress$set(message = 'Calculating overlap',
               detail = 'Please wait')

  overlaps <- vector("list", length(pdbs$id)-1)
  for(i in 2:length(pdbs$id)) {
    inds <- c(1, i)
    dv <- difference.vector(pdbs$xyz[inds, gaps.pos$f.inds, drop=FALSE])
    ov <- overlap(modes=modes, dv=dv)
    ov$id <- pdbs$id[i]
    overlaps[[i-1]] <- ov
    progress$set(value = i)
  }
  progress$close()

  rv$overlap <- overlaps
  out <- list(values=overlaps, pdbid=get_pdbid6())
  
  return(out)
})


output$overlap_plot <- renderPlot({
  ##input$run_overlap

  ov <- rv$overlap
  pdbid6 <- get_pdbid6()

  if(!is.null(ov)) {
    if(ov$pdbid==pdbid6)
      overlap_plot2(ov$values)
  }

})

overlap_plot2 <- function(o) {
  ids <- sapply(o, function(x) x$id)

  ## overlap values
  cex3 <- as.numeric(input$cex3)
  lty3 <- as.numeric(input$lty3)
  lwd3 <- as.numeric(input$lwd3)
  typ3 <- input$typ3

  if(!length(typ3)>0)
    typ3 <- "h"

  ## cummulative values
  cex4 <- as.numeric(input$cex4)
  lty4 <- as.numeric(input$lty4)
  lwd4 <- as.numeric(input$lwd4)
  typ4 <- input$typ4
  
  if(!length(typ4)>0)
    typ4 <- "b"
  
  plot.new()
  plot.window(xlim = c(1,20), ylim = c(0,1))

  ## show only overlap values for one structure
  if(input$show_overlap & !length(ids) > 1) {
    for(j in 1:length(typ3))
      points(o[[1]]$overlap, col=1, lty=lty3, lwd=lwd3, cex=cex3, type=typ3[j])
  }

  ## multiple cumulative values allowd
  if(input$show_overlap_cum) {
    for(i in 1:length(o)) {
      ov <- o[[i]]
      for(j in 1:length(typ4))
        points(ov$overlap.cum, col=i+1, lty=lty4, lwd=lwd4, cex=cex4, type=input$typ4[j])
    }
  }

  if(input$show_legend3 & length(ids) > 1)
    legend("bottomright", ids, col=1:length(ids)+1, lty=1)

  if(input$show_legend3 & length(ids) == 1 &
     input$show_overlap & input$show_overlap_cum) {
    legend("bottomright", c("overlap", "cumulative overlap"),
           col=c(1,2), lty=1)
  }


  box()
  axis(1)
  axis(2)
  mtext("Mode index", side=1, line=3)

  if(length(ids) > 1 | !input$show_overlap) {
    mtext("Cumulative Overlap", side=2, line=3)
  }
  else {
    mtext("Overlap", side=2, line=3)
  }

  main <- paste0("Overlap analysis for PDB ID ",  rv$pdbid, " (",
                 paste(rv$chains_selected, collapse=""), ")")
  sub <- paste("Difference vector(s) calculated from PDB IDs",
               paste(ids, collapse=", "))
  
  title(main)
  mtext(sub, line=.5)
}


output$overlapplot2pdf = downloadHandler(
  filename = "overlap.pdf",
  content = function(FILE=NULL) {
    ov <- rv$overlap$values
    pdf(file=FILE, width=input$width3, height=input$height3)
    overlap_plot2(ov)
    dev.off()
})

