
## set user data to store stuff
data_path <- reactive({
  dir <- paste0(format(Sys.time(), "%Y-%m-%d"), "_", randstr())
  path <- paste0(configuration$user_data, "/", dir)
  dir.create(path)
  return(path)
})

##- PDB input UI that is responsive to 'reset button' below
output$resetable_pdb_input <- renderUI({
  ## 'input$reset_pdb_input' is just used as a trigger for reset
  reset <- input$reset_pdb_input
  textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21") #)
})

## downloads and reads a PDB
read_pdb <- function(pdbid) {
  pdbid <- substr(pdbid, 1, 4)
  
  if(nchar(pdbid)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Fetching PDB',
                 detail = 'Please wait')
    progress$set(value = 2)
    
    if(configuration$pdbdir$archive) {
      pdbid <- tolower(pdbid)
      file <- paste0(configuration$pdbdir$rawfiles, "/", substr(pdbid, 2, 3),
                          "/pdb", pdbid, ".ent.gz")

      if(!file.exists(file))
        stop("PDB not found")
    }
    else {
      file <- get.pdb(pdbid, path=configuration$pdbdir$rawfiles)
    }

    progress$set(value = 3)
    progress$set(message = 'Parsing PDB',
                 detail = 'Please wait')
    pdb <- read.pdb(file, verbose=FALSE)
    progress$set(value = 5)
    
    return(pdb)
  }
  else {
    stop("Provide a 4 character PDB code")
  }
}

## returns the final PDB object
raw_pdb <- reactive({
  pdb <- read_pdb(input$pdbid)
  return(pdb)
})

## returns the final all-atom PDB object
final_pdb <- reactive({
  pdb <- raw_pdb()

  if(!length(input$chains)>0)
    stop()
  
  if(is.vector(input$chains)) {
    pdb <- trim(pdb, chain=input$chains)
  }
  
  return(pdb)
})

## return chain IDs
chain_pdb <- reactive({
  invisible(capture.output( pdb <- raw_pdb() ))
  chains <- unique(trim(pdb, "protein", elety="CA")$atom$chain)
  names(chains) <- chains
  return(chains)
})

## checkbox 
output$chain_checks <- renderUI({
  input$pdbaction
  
  chains <- chain_pdb()
  print(chains)
  checkboxGroupInput("chains", label = "Limit to chain IDs:", 
                     choices = chains, selected = chains, 
                     inline = TRUE )
})

output$pdbWebGL  <- renderWebGL({
  pdb <- final_pdb()
  view.pdb(pdb, as="overview", col="sse")
})

output$pdbSummary <- renderPrint({
  input$pdbaction
  invisible(capture.output( pdb <- final_pdb() ))
  print(pdb)
})
