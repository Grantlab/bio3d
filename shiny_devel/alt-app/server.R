library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

   fetchPDB <- reactive({
      if(nchar(input$pdbid)==4) {
         file <- get.pdb(input$pdbid)
         pdb <- read.pdb(file, verbose=FALSE)
      } else {
         stop("Provide a 4 character PDB code")
      }
      return(pdb)
   })

   output$pdbSummary <- renderPrint({
      invisible(capture.output( 
         #pdb <- checkPDB()
         pdb <- fetchPDB()
         ))
      print(pdb)
   })

   chainPDB <- reactive({
   ch <- unique(trim(fetchPDB(),"protein", elety="CA")$atom$chain)
      names(ch) <- ch
      return( ch )
   })

   output$chains <- renderPrint({
      cat( chainPDB(), sep=", ")
    })
   
   ## Would like to use on ui.r but don't yet know how!
   output$chainIDs <- reactive({ as.list(chainPDB())} )

    calcModes <- reactive({
        #pdb <- checkPDB()
        pdb <- fetchPDB()
        modes <- nma(pdb)
        return(modes)
   })

    output$modeSummary <- renderPrint({
        invisible(capture.output( modes <- calcModes()))
        print(modes)
    })

   output$fluctPlot <- renderPlot({
      #pdb <- checkPDB()
      pdb <- fetchPDB()
      modes <- calcModes()

      x <- modes$fluctuations
      main <- paste("NMA derived fluctuations for PDB id", input$pdbid)

      # draw fluctuations
      #par(mfrow=c(2,1))
      plot.bio3d(x, sse=pdb, main=main)

      #main <- paste("B-factors for PDB id", input$pdbid)
      #plot.bio3d(pdb$atom$b[ pdb$calpha ], sse=pdb, main=main)
  })
})
