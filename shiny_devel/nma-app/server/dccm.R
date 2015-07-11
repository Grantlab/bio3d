
######################
##- DCCM calculation #
######################

observeEvent(input$run_dccm, {
  rv$dccm <- calcDCCM()
})

calcDCCM <- reactive({
  message("calculating dccm")

  progress <- shiny::Progress$new()
  on.exit(progress$close())

  progress$set(message = 'Calculating DCCM',
               detail = 'Please wait ...',
               value = 0)
  
  modes <- rv$modes
  cij <- dccm(modes, ncore=1, progress=progress)
  progress$close()
  
  out <- list(cij=cij, pdbid=get_pdbid6())
  rv$dccm <- cij
  return(out)
})

##################
##- DCCM plot    #
##################
output$dccm_plot <- renderPlot({
  input$run_dccm
  pdb <- rv$final_pdb
 
  pdbid6 <- get_pdbid6()
  dccm <- rv$dccm
  cij <- dccm$cij

  if(!is.null(rv$dccm)) {
    if(!is.null(dim(cij)) & dccm$pdbid==pdbid6)
      make_dccm_plot(pdb, cij)
  }
})

make_dccm_plot <- function(pdb, cij) {
  message("make_dccm_plot")
  
  if(input$contourplot==FALSE) {
    contour <- FALSE
    cols <- bwr.colors(20)
    at <- seq(-1, 1, 0.1)
  }
  else {
    contour <- TRUE
    cols <- NULL
    at <- c(-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1)
  }
  
  if(input$sse)
    sse <- pdb
  else
    sse <- NULL
  
  main <- paste0("DCCM, PDB id ", rv$pdbid, " (",
                paste(rv$chains_selected, collapse=""), ")")
  plot(cij, sse=sse, main=main,
       contour=contour,
       col.regions=cols, colorkey=input$colorkey, at=at)
  
}



output$dccmplot2pdf = downloadHandler(
  filename = "dccm.pdf",
  content = function(FILE=NULL) {
    pdb <- rv$final_pdb
    pdf(file=FILE, width=input$width2, height=input$height2)
    make_dccm_plot(pdb, rv$dccm$cij)
    dev.off()
})




dccm_pymol <- reactive({
  path <- data_path()
  pdb <- final_pdb()
  cij <- rv$dccm$cij
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Running PyMOL script',
               detail = 'Please wait')
  
  progress$set(value = 1)
  progress$set(value = 3)
  f <- paste0(path, "/dccm.pse")
  pymol.dccm(cij, pdb, type="session", file=f)
  return(f)
})

output$dccm2pymol = downloadHandler(
  filename = "dccm_pymol.zip",
  content = function(file) {
    zip(file, dccm_pymol(), flags="-9Xj")
  })
