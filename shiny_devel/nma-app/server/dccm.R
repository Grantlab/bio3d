
######################
##- DCCM calculation #
######################
calcDCCM <- reactive({
  message("calculating dccm")
  modes <- calcModes()
  cij <- dccm(modes, ncore=1)
  return(cij)
})

output$modeSummary <- renderPrint({
  input$nmaaction
  
  invisible(capture.output( modes <- calcModes()))
  print(modes)
})

##################
##- DCCM plot    #
##################
output$dccm_plot <- renderPlot({
  make_dccm_plot()
})

make_dccm_plot <- function() {
  pdb <- final_pdb()
  
  if(input$calc_dccm) {
    
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
    
    progress <- shiny::Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating DCCM',
                 detail = 'Please wait')
    
    for (i in 1:5) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    
    main <- paste("DCCM, PDB id", input$pdbid)
    cij <- calcDCCM()
    
    for (i in 6:10) {
      progress$set(value = i)
      Sys.sleep(0.05)
    }
    
    
    plot(cij, sse=sse, main=main,
         contour=contour,
         col.regions=cols, colorkey=input$colorkey, at=at)
  }
}



output$dccmplot2pdf = downloadHandler(
  filename = "dccm.pdf",
  content = function(FILE=NULL) {
    pdf(file=FILE, width=input$width2, height=input$height2)
    make_dccm_plot()
    dev.off()
})




dccm_pymol <- reactive({
  path <- data_path()
  pdb <- final_pdb()
  cij <- calcDCCM()
  
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

output$dccm2zip = downloadHandler(
  filename = "dccm_pymol.zip",
  content = function(file) {
    zip(file, dccm_pymol(), flags="-9Xj")
  })
