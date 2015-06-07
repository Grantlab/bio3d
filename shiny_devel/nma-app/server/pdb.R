###################################
##-- User path for storing stuff  #
################################### 
data_path <- reactive({
  dir <- paste0(format(Sys.time(), "%Y-%m-%d"), "_", randstr())
  path <- paste0(configuration$user_data, "/", dir)
  dir.create(path)
  return(path)
})


###########################
##-- PDB AND BLAST INPUT  #
###########################

### set default vaules
rv <- reactiveValues()
rv$pdbid <- "4Q21"
rv$chainids <- "A"
rv$forcefield <- "calpha"
rv$cutoff <- 7

observeEvent(input$pdbid, {
  rv$pdbid <- input$pdbid
})

observeEvent(input$pdb_chains, {
  rv$chainids <- input$pdb_chains
})

observeEvent(input$cutoff, {
  rv$cutoff <- input$cutoff
})

observeEvent(input$forcefield, {
  rv$forcefield <- input$forcefield
})

observeEvent(input$reset_pdbid, {
  updateTextInput(session, "pdbid", value = "4Q21")
})

observeEvent(input$reset_nma_input, {
  updateSelectInput(session, "forcefield", value = "calpha")
  updateSelectInput(session, "cutoff", value = 7)
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

      message(file)

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

## returns PDB code (4-characters)
get_pdbid <- reactive({
  if (is.null(rv$pdbid))
    return()
  else {
    pdbid <- rv$pdbid
    
    if(!nchar(pdbid)==4) {
      stop("Provide a PDB code of 4 characters")
    }
    
    return(pdbid)
  }
})

get_pdbid6 <- reactive({
  return(paste(rv$pdbid, paste(rv$pdb_chains, collapse=""), sep="_"))
})

## returns a vector of selected chain IDs
get_chainids <- reactive({
  if (is.null(rv$chainids))
    return()
  else
    return(unlist(strsplit(rv$chainids, "")))
})

## returns the final PDB object
raw_pdb <- reactive({
  pdbid <- get_pdbid()
  pdb <- read_pdb(pdbid)
  return(pdb)
})

## returns the final all-atom PDB object
final_pdb <- reactive({
  pdb <- raw_pdb()
  chainids <- get_chainids()
  
  if(is.vector(chainids)) {
    pdb <- trim(pdb, chain=chainids)
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
output$chain_input <- renderUI({
  input$pdbaction
  
  chains <- chain_pdb()

  #checkboxGroupInput("pdb_chains", label = "Limit to chain IDs:", 
  #                   choices = chains, selected = chains, 
  #                   inline = TRUE )

  selectInput("pdb_chains", "Limit to chain IDs:",
              choices = chains, selected = chains[1], multiple=TRUE)
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
