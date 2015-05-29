
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
  dccm_plot2()
})

dccm_plot2 <- reactive({
  dccm_plot1()
})


dccm_plot1 <- function() {
  pdb <- get_pdb()
  
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
    
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
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


output$plot2pdf = downloadHandler(
  filename = 'dccm.pdf',
  content = function(file) {
    pdf(file, w=10, h=6)
    dccm_plot1()
    dev.off()
  })

dccm_pymol <- reactive({
  pdb <- get_pdb()
  
  progress <- shiny::Progress$new(session, min=1, max=5)
  on.exit(progress$close())
  
  progress$set(message = 'Calculation in progress',
               detail = 'This may take a moment...')
  
  progress$set(value = 1)
  cij <- calcDCCM()
  progress$set(value = 3)
  view.dccm(cij, pdb, launch=FALSE)
  files <- c("corr.py", "corr.inpcrd.pdb")
  zip("tmp_dccm.zip", files, flags="-FS")
  progress$set(value = 5)
  return("tmp_dccm.zip")
})

output$dccm2zip = downloadHandler(
  filename = "dccm_pymol.zip",
  content = function(file) {
    file.copy(dccm_pymol(), file)
  })
