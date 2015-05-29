
##- PDB input UI that is responsive to 'reset button' below
output$resetable_pdb_input <- renderUI({
  ## 'input$reset_pdb_input' is just used as a trigger for reset
  reset <- input$reset_pdb_input
  textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21") #)
})

## downloads and reads the raw PDB
raw_pdb <- reactive({
  id <- substr(input$pdbid, 1, 4)
  
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
get_pdb <- reactive({
  pdb <- raw_pdb()
  
  if(is.vector(input$chains)) {
    pdb <- trim(pdb, chain=input$chains)
  }
  
  return(pdb)
})

## returns the final PDB object
get_pdbCA <- reactive({
  pdb <- get_pdb()
  return(trim(pdb, "calpha"))
})

## return chain IDs
chain_pdb <- reactive({
  invisible(capture.output( pdb <- raw_pdb() ))
  chains <- unique(trim(pdb, "protein", elety="CA")$atom$chain)
  names(chains) <- chains
  return(chains)
})

output$pdbSummary <- renderPrint({
  input$pdbaction
  invisible(capture.output( pdb <- get_pdb() ))
  print(pdb)
})

## all available chains in PDB
output$chains1 <- renderPrint({
  input$pdbaction
  
  invisible(capture.output(  chains <- chain_pdb() ))
  cat( chains, sep=", ")
})

## checkbox 
output$chains2 <- renderUI({
  input$pdbaction
  
  chains <- chain_pdb()
  checkboxGroupInput("chains", label = "Limit to chain IDs:", 
                     choices = chains, inline = TRUE )
  ##selected = chains[1:length(chains)])
})

output$pdbWebGL  <- renderWebGL({
  pdb <- get_pdb()
  view.pdb(pdb, as="overview", col="sse")
})
