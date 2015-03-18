library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  finalPDB <- function() {
    ## PDB for NMA
    pdb <<- fetchPDB()

    ## unmodified PDB
    pdb2 <<- pdb
    
    chains <<- chainPDB()
    #print(chains)

    ## trim PDB
    if(is.vector(input$chains))
      pdb <<- trimPDB()

    #print(pdb)
    #print(pdb2)
    
    return(pdb)
  }
    

  fetchPDB <- reactive({
    if(nchar(input$pdbid)==4) {
      file <- get.pdb(input$pdbid)
      pdb <- read.pdb(file, verbose=FALSE)
      return(pdb)
    }
    else {
      stop("Provide a 4 character PDB code")
    }
  })
    
  trimPDB <- reactive({
    ##message("Trimming PDB")
    if(is.vector(input$chains)) {
      pdb <- trim(pdb, chain=input$chains)
      return(pdb)
    }
  })
  
  chainPDB <- reactive({
    invisible(capture.output( pdb2 <- fetchPDB() ))
    chains <- unique(trim(pdb2, "protein", elety="CA")$atom$chain)
    names(chains) <- chains
    return(chains)
  })
  
  output$pdbSummary <- renderPrint({
    input$pdbaction
     
    invisible(capture.output( pdb <- finalPDB() ))
    print(pdb)
  })

  ## all available chains in PDB
  output$chains1 <- renderPrint({
    input$pdbaction
    
    invisible(capture.output(  chains <- chainPDB() ))
     cat( chains, sep=", ")
    })


  ## checkbox 
  output$chains2 <- renderUI({
    input$pdbaction
    
    chains <- chainPDB()
    checkboxGroupInput("chains", label = "Limit to chain IDs:", 
                       choices = chains,
                       selected = chains[1:length(chains)])
  })
  
  
  calcModes <- reactive({
    pdb <- finalPDB()
    
    temp <- as.numeric(input$temp)
    if(temp==0)
      temp <- NULL
    
    if(input$forcefield %in% c("calpha", "sdenm"))
      modes <<- nma(pdb, ff=input$forcefield, mass=input$mass, temp=temp)

    if(input$forcefield %in% c("anm", "pfanm"))
      modes <<- nma(pdb, ff=input$forcefield, cutoff=input$cutoff,
                    mass=input$mass, temp=temp)
    
    return(modes)
    
  })

  calcDCCM <- reactive({
    modes <- calcModes()
    cij <- dccm(modes)
    return(cij)
  })
  
  output$modeSummary <- renderPrint({
    input$nmaaction
    
    invisible(capture.output( modes <- calcModes()))
    print(modes)
  })

  output$fluctPlot <- renderPlot({
    input$nmaaction
    
    pdb <- finalPDB()
    modes <- calcModes()

    x <- modes$fluctuations
      
    ## draw fluctuations
    if(input$fluxs1 & input$fluxs2 )
      par(mfrow=c(2,1))

    if(input$fluxs1) {
      main <- paste("NMA derived fluctuations for PDB id", input$pdbid)
      plot.bio3d(x, sse=pdb, main=main)
    }

    if(input$fluxs2) {
      main <- paste("B-factors for PDB id", input$pdbid)
      plot.bio3d(pdb$atom$b[ pdb$calpha ], sse=pdb, main=main)
    }
  })


  
  output$dccmPlot <- renderPlot({
    if(input$cijs) {
    
    #if(input$style=="Red-blue") {
    #  contour <- FALSE
    #  cols <- bwr.colors(20)
    #  at <- seq(-1, 1, 0.1)
    #    }
    #else {
    #  contour <- TRUE
    #  cols <- NULL
    #  at <- c(-1, -0.75, -0.5,  -0.25, 0.25, 0.5, 0.75, 1)
    #}
    
    #if(input$sse)
    #  sse <- pdb
    #else
    #  sse <- NULL
    
    progress <- shiny::Progress$new(session, min=1, max=10)
    on.exit(progress$close())
    
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:5) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    
    main <- paste("DCCM, PDB id", input$pdbid)
    cij <- calcDCCM()
    
    for (i in 6:10) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    
    
    plot(cij, sse=pdb, main=main)
         #contour=contour,
         #col.regions=cols, colorkey=input$colorkey, at=at)
  }
    
  })
  
})
