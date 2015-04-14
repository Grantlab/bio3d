  ############
  ## NMA
  ###########
  nma1 <- reactive({
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Calculating normal modes',
                 detail = 'This may take some time...')
    progress$set(value = 2)

    pdbs <- align()
    modes <- nma(pdbs, fit=TRUE)
    progress$set(value = 5)
    return(modes)
  })

  output$nma_plot <- renderPlot({
    pdbs <- align()
    modes <- nma1()
    plot(modes, pdbs)
  })
  
  output$rmsip_plot <- renderPlot({
    pdbs <- align()
    modes <- nma1()
    rownames(modes$rmsip) <- basename.pdb(rownames(modes$rmsip))
    colnames(modes$rmsip) <- basename.pdb(colnames(modes$rmsip))
    heatmap(1-modes$rmsip, distfun=as.dist, symm=TRUE)
  })
