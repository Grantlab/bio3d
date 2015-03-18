library(shiny)
library(bio3d)

# Define server logic required to calculate normal modes
# and plot fluctuations
shinyServer(function(input, output, session) {


    fetchPDB <- reactive({
        if(nchar(input$pdbid)==4) {
            print(input)
            file <- get.pdb(input$pdbid)
            pdb <- read.pdb(file)
        }
        else {
            stop("Provide a 4 character PDB code")
        }
        return(pdb)
    })


    checkPDB <- reactive({
        pdb <- fetchPDB()

        if(input$chain!="") {
            chain <- input$chain
            chain <- gsub(" ", "", chain)
            chain <- gsub(",", "", chain)
            chain <- unlist(strsplit(chain, ""))

            if(chain[1] != "") {
                sele <- atom.select(pdb, chain=chain)

                if(!length(sele$atom)>0)
                    stop("insufficient atoms in selection")

                pdb <- trim.pdb(pdb, sele)
            }
        }

        max_atoms <- 500
        min_atoms <- 2

        if(sum(pdb$calpha) > max_atoms)
            stop(paste(">", max_atoms, "residues not allowed"))
        if(sum(pdb$calpha) < min_atoms)
            stop(paste("<", min_atoms, "residues not allowed"))

        return(pdb)
    })


    calcModes <- reactive({
        pdb <- checkPDB()

        temp <- as.numeric(input$temp)
        if(temp==0)
            temp <- NULL

        if(input$mass=="Yes")
            mass <- TRUE
        else
            mass <- FALSE

        if(input$forcefield %in% c("calpha", "sdenm"))
            modes <- nma(pdb, ff=input$forcefield, mass=mass, temp=temp)

        if(input$forcefield %in% c("anm", "pfanm"))
            modes <- nma(pdb, ff=input$forcefield, cutoff=input$cutoff,
                         mass=mass, temp=temp)

        return(modes)

   })

    calcDCCM <- reactive({
        modes <- calcModes()
        cij <- dccm(modes)
        return(cij)
    })

    output$fluctPlot <- renderPlot({
        pdb <- checkPDB()
        modes <- calcModes()

        x <- modes$fluctuations
        main <- paste("NMA on PDB id", input$pdbid)

        # draw fluctuations
        par(mfrow=c(2,1))
        plot.bio3d(x, sse=pdb, main=main)

        main <- paste("B-factors for PDB id", input$pdbid)
        plot.bio3d(pdb$atom$b[ pdb$calpha ], sse=pdb, main=main)
  })




    output$dccmPlot <- renderPlot({
        pdb <- checkPDB()

        if(input$style=="Red-blue") {
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
            Sys.sleep(0.1)
        }


        main <- paste("DCCM, PDB id", input$pdbid)
        cij <- calcDCCM()

        for (i in 6:10) {
            progress$set(value = i)
            Sys.sleep(0.1)
        }


        plot(cij, sse=sse, main=main, contour=contour,
             col.regions=cols, colorkey=input$colorkey, at=at)

  })

    output$pdbSummary <- renderPrint({
        invisible(capture.output( pdb <- checkPDB()))
        print(pdb)
    })

    output$modeSummary <- renderPrint({
        invisible(capture.output( modes <- calcModes()))
        print(modes)
    })

})

