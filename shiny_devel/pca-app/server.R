library(bio3d)
library(lattice)
library(shiny)

## Define server logic for PCA shiny demo

shinyServer(function(input, output, session) {

  ##- PDB input UI that is responsive to 'reset button' below
  output$resetable_pdb_input <- renderUI({
    ## 'input$reset_pdb_input' is just used as a trigger for reset
    reset <- input$reset_pdb_input
    textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21") #)
  })
  
  finalPDB <- function() {
    pdb <<- fetchPDB()
    pdb.orig <<- pdb
    chains <<- chainPDB()

    if(is.vector(input$chains)) {
      pdb <<- trimPDB()
    }

    pdb.ca <<- trim(pdb, "calpha")
    return(pdb)
  }
    

  fetchPDB <- reactive({
    print(input$pdbid)
    
    if(nchar(input$pdbid)==4) {
      progress <- shiny::Progress$new(session, min=1, max=5)
      on.exit(progress$close())

      progress$set(message = 'Fetching PDB',
                   detail = 'This should be quick...')
      progress$set(value = 2)
      file <- get.pdb(input$pdbid, path="raw_files")
      progress$set(value = 3)

      progress$set(message = 'Parsing PDB',
                   detail = 'Please wait...')
      pdb <- read.pdb(file, verbose=FALSE)
      for(i in 4:5)
        progress$set(value = i)
      
      return(pdb)
    }
    else {
      stop("Provide a 4 character PDB code")
    }
  })
    
  trimPDB <- reactive({
    if(is.vector(input$chains)) {
      pdb <- trim(pdb, chain=input$chains)
      return(pdb)
    }
  })

  chainPDB <- reactive({
    invisible(capture.output( pdb <- fetchPDB() ))
    chains <- unique(trim(pdb, "protein", elety="CA")$atom$chain)
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
  ##-- BLAST / HMMER     #
  ########################
  blast <- reactive({
    pdb <- finalPDB()

    if(toupper(input$pdbid)=="4Q21")
      load(file="4Q21_hmmer.RData")
    else
      hmm <- hmmer(pdbseq(pdb))

    hmm <<- hmm
    return(hmm)
  })

  output$blast_plot <- renderPlot({
    plot(blast(), cutoff=input$cutoff)
  })

  hits1 <- function() {
    hits <- plot(blast(), cutoff=input$cutoff)$acc
    
    if(!is.null((input$limit_hits))) {
      n <- as.numeric(input$limit_hits)
      if(n > 0) {
        if(n>length(hits))
          n <- length(hits)
        hits <- hits[1:n]
      }
    }
    return(hits)
  }
    
  output$hits2 <- renderPrint({
    invisible(capture.output( hits <- hits1()) )
    print(hits)
  })

  output$hits_slider <- renderUI({
    hmm <- blast()
    invisible(capture.output(nhits <- plot(hmm)$acc))
    sliderInput("limit_hits", "Limit hits:",
                min = 1, max = length(nhits), value = 5)
    })
  
  
  output$pdbids_checkboxgroup <- renderUI({
    hits <- hits1()
    names(hits) <- hits
    checkboxGroupInput("pdb_ids", "PDB IDs:",
                       hits, selected=hits, inline=TRUE)
  })


  #####################
  ##-- Align
  ####################
  
  fetch_pdbs <- reactive({
    ids <- input$pdb_ids
    unq <- unique(substr(ids, 1,4))
    
    progress <- shiny::Progress$new(session, min=1, max=length(unq))
    
    
    progress$set(message = 'Fetching PDBs',
                 detail = 'This may take some time...')
    
    ##raw.files <- get.pdb(ids, gzip=TRUE)
    raw.files <- vector("character", length(unq))
    for(i in 1:length(unq)) {
      raw.files[i] <- get.pdb(unq[i], path="raw_files", gzip=TRUE)
      progress$set(value = i)
    }
    progress$close()
    
    progress <- shiny::Progress$new(session, min=1, max=length(ids))
    
    
    progress$set(message = 'Splitting PDBs',
                 detail = 'This may take some time...')

    ##files <- pdbsplit(raw.files, ids)
    
    ## this is possibly error prone
    files <- vector("character", length(ids))
    for(i in 1:length(unq)) {
      inds <- grep(unq[i], ids)
      
      files[inds] <- pdbsplit(raw.files[i], ids[inds])
      progress$set(value = i)
    }
    
    progress$close()
    return(files)
  })

  align <- reactive({
    files <- fetch_pdbs()
        
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Aligning PDBs',
                 detail = 'This may take some time...')
    progress$set(value = 2)

    pdbs <- pdbaln(files)
    progress$set(value = 5)
    
    rownames(pdbs$ali) <- basename.pdb(rownames(pdbs$ali))
    pdbs <<- pdbs
    return(pdbs)
  })
  
   output$alignment <- renderPrint({
     input$action_fetchpdbs
     
     invisible(capture.output( aln <- align() ))
     print(aln)
   })

  ####################
  ##-- PCA
  ####################
  find_core <- reactive({
    pdbs <- align()

    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Finding invariant core',
                 detail = 'This may take some time...')
    progress$set(value = 2)

    core <<- core.find(pdbs)
    progress$set(value = 5)
    return(core)
  })

  fit <- reactive({
    core <- find_core()
    pdbs$xyz <<- pdbfit(pdbs, core)
    return(pdbs)
  })
  
  pca1 <- reactive({
    pdbs <- fit()
    pc <<- pca(pdbs)
    return(pc)
  })

  clust <- reactive({
    pdbs <- fit()
    rd <- rmsd(pdbs)
    hc <- hclust(as.dist(rd))
    return(hc)
  })

  clustgrps <- reactive({
    hc <- clust()
    grps <- cutree(hc, k=input$nclust)
    return(grps)
  })

  output$checkboxgroup_label_ids <- renderUI({
    
    ids <- input$pdb_ids
    names(ids) <- ids
    checkboxGroupInput("label_ids", "PDB IDs:",
                       ids, selected=ids, inline=TRUE)
    
  })

  output$pca_plot <- renderPlot({
    invisible(capture.output( pc <- pca1() ))

    col <- 1
    if(input$nclust>1)
      col <- clustgrps()
              
    if(input$screeplot)
      par(mfrow=c(1,2))
    
    xax <- as.numeric(input$pcx)
    yax <- as.numeric(input$pcy)

    p <- paste0("PC", c(xax, yax),
                " (", round((pc$L[c(xax, yax)]/sum(pc$L)) * 
                            100, 2), "%)")
    
    plot(pc$z[, xax], pc$z[, yax],
         col=col, pch=16,
         xlab=p[1], ylab=p[2])
    
    abline(h = 0, col = "gray", lty = 2)
    abline(v = 0, col = "gray", lty = 2)

    if(input$labelplot) {
      if(length(input$label_ids)>0) {
        inds <- unlist(lapply(input$label_ids, grep, pdbs$id))
        text(pc$z[inds, xax], pc$z[inds, yax],
             labels=basename.pdb(pdbs$id[inds]),
             pos=1, offset=input$offset)
      }
    }
    
    if(input$screeplot) {
      plot.pca.scree(pc$L)
    }
    
  })

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
  
})
