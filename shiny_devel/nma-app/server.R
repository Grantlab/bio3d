#options(rgl.useNULL=TRUE) ## Required for a headless server
library(shiny)
library(shinyRGL)
library(rgl)

## Define server logic for NMA shiny demo
## TODO: Simplify!

shinyServer(function(input, output, session) {


  ##- PDB input UI that is responsive to 'reset button' below
  output$resetable_pdb_input <- renderUI({
    ## 'input$reset_pdb_input' is just used as a trigger for reset
    reset <- input$reset_pdb_input
    textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21") #)
  })



  finalPDB <- function() {
    pdb <<- fetchPDB()

    ## unmodified PDB, CONFUSING
    pdb2 <<- pdb
    chains <<- chainPDB()

    if(is.vector(input$chains))
      pdb <<- trimPDB()

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
                       choices = chains)
                       ##selected = chains[1:length(chains)])
  })
  


  ########################
  ##-- NMA Calculation   #
  ########################s
  output$resetable_nma_input <- renderUI({
    ## used as a trigger for reset
    reset <- input$reset_nma_input
    div(
    selectInput("forcefield", "Choose a forcefield:",
      choices = c("calpha", "anm", "pfanm", "sdenm")),


    sliderInput("cutoff", "Cutoff value:",
      min = 6, max = 50, value = 10),

    sliderInput("temp", "Temperature scaling:",
      min = 0, max = 350, value = 300),

    checkboxInput("mass", "Mass-weighting", value=TRUE) )
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
  
  ######################
  ##- DCCM calculation #
  ######################
  calcDCCM <- reactive({
    message("calculating dccm")
    modes <- calcModes()
    cij <- dccm(modes)
    return(cij)
  })
  
  output$modeSummary <- renderPrint({
    input$nmaaction
    
    invisible(capture.output( modes <- calcModes()))
    print(modes)
  })

  ##########################
  ##- Fluctuation plot     #
  ##########################
  output$fluct_plot <- renderPlot({
    fluct_plot2()
  })
  
  fluct_plot2<- reactive({
    fluct_plot1()
  })
  
  fluct_plot1 <- function() {
    pdb <- finalPDB()
    modes <- calcModes()
    x <- modes$fluctuations
    
    if(input$fluxs1 & input$fluxs3) { 
      main <- paste("NMA derived fluctuations for PDB id", input$pdbid)
      plot.bio3d(x, sse=pdb, main=main, resno=pdb,
                 col=input$col1, type=input$typ1,
                 pch=input$pch1, lty=input$lty1, cex=input$cex1, lwd=input$lwd1)
    }
    if(input$fluxs2) {
      if(input$fluxs3) {
        par(new=TRUE)
        plot.bio3d(pdb$atom$b[ pdb$calpha ], resno=pdb, axes=F, xlab="",ylab="",
             pch=input$pch2, col=input$col2, typ=input$typ2, lty=input$lty2, lwd=input$lwd2)
        axis(4, col=input$col2)
      }
      else {
        main <- paste("B-factors for PDB id", input$pdbid)
        plot.bio3d(pdb$atom$b[ pdb$calpha ], sse=pdb, main=main, resno=pdb,
             pch=input$pch2, col=input$col2, typ=input$typ2, lty=input$lty2, cex=input$cex2,
                   lwd=input$lwd2)
      }
    }
  }
  
  output$plot1pdf = downloadHandler(
    filename = 'fluct.pdf',
    content = function(file) {
      pdf(file, w=10, h=6)
      fluct_plot1()
      dev.off()
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
    message("running plot function")
    
    if(input$cijs) {
    
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
                 detail = 'This may take a moment...')
    
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
  
    
  

  
  
})
