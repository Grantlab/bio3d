library(shiny)
library(bio3d)

# Define server logic required to calculate normal modes
# and plot fluctuations 
shinyServer(function(input, output) {

  output$distPlot <- renderPlot({
    pdbid <- input$pdbid
    
    if(nchar(pdbid)==4) {
      # read pdb
      #pdb <- read.pdb( system.file("examples/1hel.pdb", package="bio3d") )
      pdb <- read.pdb(pdbid)

      if(!is.pdb(pdb))
        stop(paste("pdb id", pdbid, "not found"))

      if(sum(pdb$calpha)>300)
        stop(">300 residues not allowed")
      if(sum(pdb$calpha)<10)
        stop("<10 residues not allowed")

      temp <- as.numeric(input$temp)
      if(temp==0)
        temp <- NULL

      if(input$mass=="Yes")
        mass <- TRUE
      else
        mass <- FALSE
       
      # Calculate normal modes
      if(input$forcefield %in% c("calpha", "sdenm"))
        modes <- nma(pdb, ff=input$forcefield, mass=mass, temp=temp)
      if(input$forcefield %in% c("anm", "pfanm"))
        modes <- nma(pdb, ff=input$forcefield, cutoff=input$cutoff, mass=mass, temp=temp)


      if(input$plot=="Fluctuations") {
        x <- modes$fluctuations
        main <- paste("NMA on PDB id", pdbid)

        # draw fluctuations
        plot.bio3d(x, sse=pdb, main=main)
      }
      else {
        main <- paste("DCCM, PDB id", pdbid)        
        cij <- dccm(modes)
        plot(cij, sse=pdb, main=main)
      }
    }
    else {
      stop("Enter a valid PDB of length 4")
    }

    
  })
})

