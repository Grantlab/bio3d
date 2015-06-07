########################
##-- NMA Calculation   #
########################s

calcModes <- reactive({
  pdb <- final_pdb()
  ff <- rv$forcefield
  cutoff <- rv$cutoff
  
  if(sum(pdb$calpha)>600)
    stop("maximum 600 C-alpha atoms for NMA")
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Calculating normal modes',
               detail = 'Please wait')
  progress$set(value = 2)
  
  if(ff %in% c("calpha", "sdenm", "reach"))
    modes <- nma(pdb, ff=ff, mass=TRUE, temp=300)
  
  if(ff %in% c("anm", "pfanm"))
    modes <- nma(pdb, ff=ff, cutoff=cutoff,
                  mass=FALSE, temp=NULL)
  progress$set(value = 5)
  
  return(modes)
})

##########################
##- Fluctuation plot     #
##########################
output$fluct_plot <- renderPlot({
  make_fluct_plot()
})

make_fluct_plot <- function() {
  pdb <- final_pdb()
  modes <- calcModes()

  mode.inds <- input$mode_inds
  
  if(any(mode.inds == "all")) {
    x <- modes$fluctuations
    modes_str <- ""
  }
  else {
    x <- fluct.nma(modes, as.numeric(mode.inds))
    bs <- bounds(as.numeric(mode.inds))
    modes_str <- paste(apply(bs, 1, function(x) if(x[1]==x[2]) x[1] else paste(x[1], x[2], sep="-")), collapse=", ")
  }
  

  if(input$fluxs2 & input$fluxs3)
    par(mar=c(5, 4, 4, 5))
  else
    par(mar=c(5, 4, 4, 2))
  
  if(input$fluxs3) { 
    main <- paste0("NMA derived fluctuations for PDB id ", rv$pdbid, " (",
                  paste(rv$chainids, collapse=""), ")")
    sub <- paste("Mode(s)", mode.inds)
    plot.bio3d(x, sse=pdb, main=main, resno=pdb,
               xlab="Residue No.", ylab="Fluctions (Å^2)", 
               col=input$col1, type=input$typ1,
               pch=input$pch1, lty=input$lty1, cex=input$cex1, lwd=input$lwd1)

    if(modes_str!="") {
      if(length(mode.inds)>1)
        mtext(paste("Modes ", modes_str), line=.5)
      else
        mtext(paste("Mode ", modes_str), line=.5)
    }
  }
  if(input$fluxs2) {
    if(input$fluxs3) {
      par(new=TRUE)
      plot.bio3d(pdb$atom$b[ pdb$calpha ], resno=pdb, axes=F, xlab="",ylab="",
                 pch=input$pch2, col=input$col2, typ=input$typ2,
                 lty=input$lty2, lwd=input$lwd2)
      axis(4, col=input$col2)
      mtext("B-factors (Å)", side=4, line=3)
    }
    else {
      main <- paste("B-factors for PDB id", input$pdbid)
      plot.bio3d(pdb$atom$b[ pdb$calpha ], sse=pdb, main=main, resno=pdb,
                 xlab="Residue No.", ylab="B-factors (Å)", 
                 pch=input$pch2, col=input$col2, typ=input$typ2, lty=input$lty2, cex=input$cex2,
                 lwd=input$lwd2)
    }
  }
}


output$fluctplot2pdf = downloadHandler(
  filename = "nma_fluctuations.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, width=input$width1, height=input$height1)
    make_fluct_plot()
    dev.off()
})




####################
## Trajectory
###################

nma2pdb <- reactive({
  path <- data_path()
  pdb <- final_pdb()
  pdb <- trim(pdb, "calpha")
  modes <- calcModes()
  mag <- as.numeric(input$mag)
  step <- mag/8
  
  i <- as.numeric(input$mode_choice)
  fname  <- paste0(path, '/', 'mode', input$mode_choice, '.pdb')
  x <- mktrj(modes, mode=i, file=fname,
             mag=mag, step=step,
             b=modes$fluctuations,
             resno=pdb$atom$resno,
             resid=pdb$atom$resid,
             chain=pdb$atom$chain)
  
  return(fname)
})


output$trj2zip = downloadHandler(
  filename=function() {
    paste0('mode', input$mode_choice, '.pdb.zip')
  },
  content = function(file) {
    zip(file, files=nma2pdb(), flags = "-9Xj")
  })

####################
## Vector field
###################

make_nma_pse <- reactive({
  path <- data_path()
  pdb <- final_pdb()
  pdb <- trim(pdb, "calpha")
  modes <- calcModes()
  mag <- as.numeric(input$mag)
  
  outf <- paste0(path, "/mode", as.numeric(input$mode_choice), ".pse")
  file <- pymol.modes(modes, mode=as.numeric(input$mode_choice), type="session",
                file=outf, scale=mag)
  return(outf)
})

output$nma2pymol = downloadHandler(
  filename=function() {
    paste0('mode', input$mode_choice, '.pse.zip')
  },
  content = function(file) {
    zip(file, files=make_nma_pse(), flags = "-9Xj")
})


## summary
output$modeSummary <- renderPrint({
  invisible(capture.output( modes <- calcModes()))
  print(modes)
})
