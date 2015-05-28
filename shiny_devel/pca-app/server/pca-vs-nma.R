
refit <- reactive({
  pdbs <- align()
  pdbs$xyz <- pdbfit(pdbs)
  return(pdbs)
})

pca3 <- reactive({
  pdbs <- refit()
  pc <- pca(pdbs)
  return(pc)
})

nma3 <- reactive({
  if(input$rm.gaps == TRUE) {
    modes <- nma2()
  }
  else {
    pdbs <- refit()

    progress <- shiny::Progress$new()
    on.exit(progress$close())

    progress$set(message = 'Calculating normal modes',
                 detail = 'Please wait',
                 value = 0)
    modes <- nma(pdbs, fit=TRUE, rm.gaps=TRUE, progress=progress)
  }
  return(modes)
})


pcanma_rmsip <- reactive({
  pc <- pca3()
  modes <- nma3()

  if(is.null(input$viewStruct))
    stop()

  r <- rmsip(pc$U, modes$U.subspace[,, as.numeric(input$viewStruct)])
  return(r)
})

output$rmsip_plot3 <- renderPlot({
  r <- pcanma_rmsip()
  plot(r, xlab="PCA", ylab="NMA")
})

output$rmsip_print3 <- renderPrint({
  r <- pcanma_rmsip()
  r$overlap <- round(r$overlap, 2)
  #class(r) <- NULL
  print(r$overlap)
})

output$struct_dropdown3 <- renderUI({
  pdbs <- refit()
  ids <- 1:length(pdbs$id)
  names(ids) <-  pdbs$lab
  selectInput('viewStruct', 'Compare NMs of structure:',
              choices=ids)
})
 
