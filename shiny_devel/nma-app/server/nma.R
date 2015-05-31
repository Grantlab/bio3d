########################
##-- NMA Calculation   #
########################s
output$resetable_nma_input <- renderUI({
  ## used as a trigger for reset
  reset <- input$reset_nma_input
  div(
    selectInput("forcefield", "Choose a forcefield:",
                choices = c("calpha", "sdenm", "reach", "anm", "pfanm")),
    
    sliderInput("cutoff", "Cutoff value:",
                min = 7, max = 50, value = 15)
    )
})

calcModes <- reactive({
  pdb <- get_pdb()
  
  if(sum(pdb$calpha)>600)
    stop("maximum 600 C-alpha atoms for NMA")
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Calculating normal modes',
               detail = 'This may take while...')
  progress$set(value = 2)
  
  if(input$forcefield %in% c("calpha", "sdenm", "reach"))
    modes <- nma(pdb, ff=input$forcefield, mass=TRUE, temp=300)
  
  if(input$forcefield %in% c("anm", "pfanm"))
    modes <- nma(pdb, ff=input$forcefield, cutoff=input$cutoff,
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
  pdb <- get_pdb()
  modes <- calcModes()
  x <- modes$fluctuations

  if(input$fluxs2 & input$fluxs3)
    par(mar=c(5, 4, 4, 5))
  else
    par(mar=c(5, 4, 4, 2))
  
  if(input$fluxs3) { 
    main <- paste("NMA derived fluctuations for PDB id", input$pdbid)
    plot.bio3d(x, sse=pdb, main=main, resno=pdb,
               xlab="Residue No.", ylab="Fluctions (Å^2)", 
               col=input$col1, type=input$typ1,
               pch=input$pch1, lty=input$lty1, cex=input$cex1, lwd=input$lwd1)
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
mktrj2 <- reactive({
  pdb <- get_pdb()
  modes <- calcModes()

  progress <- shiny::Progress$new(session, min=1, max=input$trj_nmodes)
  on.exit(progress$close())
  
  progress$set(message = 'Calculation in progress',
               detail = 'This may take a moment...')
  
  
  files <- c()
  for(i in 7:(6+input$trj_nmodes)) {
    f <- paste0("mode_", i, ".pdb")
    x <- mktrj(modes, mode=i, file=f,
               b=modes$fluctuations,
               resno=pdb.ca$atom$resno,
               resid=pdb.ca$atom$resid,
               chain=pdb.ca$atom$chain)
    
    files <- c(files, f)
    progress$set(value = i)
  }
  print(files)
  zip("tmp_trj.zip", files, flags="-FS")
  return("tmp_trj.zip")
})

output$trj2zip = downloadHandler(
  filename = "trj.zip",
  content = function(file) {
    file.copy(mktrj2(), file)
  })


