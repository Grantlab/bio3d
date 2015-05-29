
############################
## Overlap analysis
###########################

raw_pdb2 <- reactive({
  ## downloads and reads the raw PDB
  id <- substr(input$pdbid2, 1, 4)
  
  if(nchar(id)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Fetching PDB',
                 detail = 'This should be quick...')
    progress$set(value = 2)
    
    if(configuration$pdbdir$archive) {
      ids <- paste0(tolower(substr(id, 1, 4)), "_", substr(id, 6, 6))
      raw.files <- paste0(configuration$pdbdir$splitfiles, "/", substr(ids, 2, 3), "/pdb", ids, ".ent.gz")
      
      if(!file.exists(raw.files))
        stop("PDB not found")
      
      pdb <- read.pdb(raw.files)
    }
    else {
      file <- get.pdb(id)
    }
    progress$set(value = 3)
    
    progress$set(message = 'Parsing PDB',
                 detail = 'Please wait...')
    pdb <- read.pdb(file, verbose=FALSE)
    progress$set(value = 5)

    return(pdb)
  }
  else {
    stop("Provide a 4 character PDB code")
  }
})

## returns the final PDB object
get_pdb2 <- reactive({
  pdb <- raw_pdb2()
  
  if(is.vector(input$chains2)) {
    pdb <- trim(pdb, chain=input$chains2)
  }
  
  return(pdb)
})



## return chain IDs
chain_pdb2 <- reactive({
  invisible(capture.output( pdb <- raw_pdb2() ))
  chains <- unique(trim(pdb, "protein", elety="CA")$atom$chain)
  names(chains) <- chains
  return(chains)
})


## all available chains in PDB
output$chains3 <- renderPrint({
  invisible(capture.output(  chains <- chain_pdb2() ))
  cat( chains, sep=", ")
})

## checkbox 
output$chains4 <- renderUI({
  chains <- chain_pdb2()
  checkboxGroupInput("chains2", label = "Limit to chain IDs:", 
                     choices = chains, inline = TRUE)
})


calc_overlap <- reactive({
  pdb1 <- get_pdb()
  pdb2 <- get_pdb2()
  pdb1 <- trim(pdb1, "calpha")
  pdb2 <- trim(pdb2, "calpha")
  
  id1 <- toupper(paste(input$pdbid, input$chains, sep="_"))
  id2 <- toupper(paste(input$pdbid, input$chains, sep="_"))

  t1 <- tempfile()
  t2 <- tempfile()

  write.pdb(pdb1, file=t1)
  write.pdb(pdb2, file=t2)

  pdbs <- pdbaln(c(t1, t2))
  pdbs$xyz <- pdbfit(pdbs)

  gaps.pos <- gap.inspect(pdbs$xyz)
  gaps.res <- gap.inspect(pdbs$ali)
  inds <- as.select(gaps.res$f.inds)

  if(input$forcefield %in% c("calpha", "sdenm", "reach"))
    modes <- nma(pdb1, ff=input$forcefield, mass=FALSE, temp=300,
                 outmodes = inds)
  
  if(input$forcefield %in% c("anm", "pfanm"))
    modes <- nma(pdb1, ff=input$forcefield, cutoff=input$cutoff,
                 mass=FALSE, temp=NULL, outmodes = inds)
  
  xyz <- pdbs$xyz[, gaps.pos$f.inds, drop=FALSE]
  dv <- difference.vector(xyz)
  ov <- overlap(modes=modes, dv=dv)
  return(ov)
})


output$overlap_plot <- renderPlot({
  overlap_plot2()
})

overlap_plot2 <- function() {
  o <- calc_overlap()
  plot(o$overlap, type='h', ylim=c(0,1), xlab="Overlap", ylab="Mode index")
  points(o$overlap)
  lines(o$overlap.cum, type='b', col='red')
}

