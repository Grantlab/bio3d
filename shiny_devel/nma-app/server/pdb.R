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
##-- Reactive values
###########################

### set default vaules
rv <- reactiveValues()
rv$pdbid <- "4Q21"
rv$chains_all <- "A"
rv$chains_selected <- "A"

rv$forcefield <- "calpha"
rv$cutoff <- 7

rv$dccm <- NULL
rv$blast <- NULL
rv$overlap <- NULL

rv$modes <- readRDS("4q21_modes.RDS")
rv$raw_pdb <- readRDS("4q21_pdb.RDS")
rv$final_pdb <- readRDS("4q21_pdb.RDS")


observeEvent(input$pdbid, {
  if(nchar(input$pdbid)>3) {
    rv$pdbid <- substr(input$pdbid, 1, 4)

    
    updateButton(session, "run_dccm", disabled = FALSE)
    updateButton(session, "run_blast", disabled = FALSE)
    updateButton(session, "run_overlap", disabled = TRUE,
                 icon = icon("gears")) ##, label = " Calculate overlap")
    
    rv$dccm <- NULL
    rv$blast <- NULL
    rv$overlap <- NULL

    if(toupper(rv$pdbid) == "4Q21") {
      message("ID is 4q21 - loading data")
      rv$raw_pdb <- readRDS("4q21_pdb.RDS")
      rv$final_pdb <- readRDS("4q21_pdb.RDS")
      rv$modes <- readRDS("4q21_modes.RDS")

      rv$chains_all <- "A"
      rv$chains_selected <- "A"
    }
    else {
      rv$raw_pdb <- raw_pdb()

      chains <- chain_pdb()
      rv$chains_all <- chains
      rv$chains_selected <- chains[1]
      
      rv$final_pdb <- final_pdb()

      if(!pdb_isok1()) {
        updateButton(session, "run_dccm", disabled = TRUE)
        updateButton(session, "run_blast", disabled = TRUE)
        updateButton(session, "run_overlap", disabled = TRUE)
      }
      else {
        rv$modes <- calcModes()
      }
    }
  }
})

observeEvent(input$chainid, {
  if(length(input$chainid) > 0) {
    rv$chains_selected <- input$chainid
    message(rv$chains_selected)

    if(toupper(rv$pdbid) == "4Q21") {
      rv$modes <- readRDS("4q21_modes.RDS")
    }
    else {
      rv$final_pdb <- final_pdb()
      
      if(!pdb_isok1()) {
        updateButton(session, "run_dccm", disabled = TRUE)
        updateButton(session, "run_blast", disabled = TRUE)
        updateButton(session, "run_overlap", disabled = TRUE)
      }
      else {
        rv$modes <- calcModes()
      }
    }
  }
})

observeEvent(input$forcefield, {
  rv$forcefield <- input$forcefield
  if(toupper(rv$pdbid) == "4Q21" & input$forcefield == "calpha") {
    rv$modes <- readRDS("4q21_modes.RDS")
  }
  else {
    rv$modes <- calcModes()
  }
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
  updateSelectInput(session, "forcefield", selected = "calpha")
  updateSliderInput(session, "cutoff", value = 15)
})

observeEvent(input$run_dccm, {
  updateButton(session, "run_dccm", disabled = TRUE)
})

observeEvent(input$run_blast, {
  updateButton(session, "run_blast", disabled = TRUE)
  updateButton(session, "run_overlap", disabled = FALSE)
})

observeEvent(input$run_overlap, {
  updateButton(session, "run_overlap", disabled = FALSE,
               icon = icon("refresh")) ##, label = " Re-calculate overlap")
})

output$dccm_isdone <- reactive({
  return(!is.null(rv$dccm))
})
outputOptions(output, 'dccm_isdone', suspendWhenHidden=FALSE)

output$blast_isdone <- reactive({
  return(!is.null(rv$blast))
})
outputOptions(output, 'blast_isdone', suspendWhenHidden=FALSE)

output$overlap_isdone <- reactive({
  return(!is.null(rv$overlap))
})
outputOptions(output, 'overlap_isdone', suspendWhenHidden=FALSE)

## check if PDB is ok
pdb_isok1 <- reactive({
  pdb <- rv$final_pdb

  message(sum(pdb$calpha))

  if(sum(pdb$calpha) > 600 | sum(pdb$calpha) < 10)
    return(FALSE)
  else
    return(TRUE)
})

output$pdb_isok <- reactive({
  pdb_isok1()
})
outputOptions(output, 'pdb_isok', suspendWhenHidden=FALSE)


###########################
##-- PDB AND BLAST INPUT  #
###########################



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
      pdbid <- format_pdbids(pdbid, casefmt=tolower)

      file <- paste0(configuration$pdbdir$rawfiles, "/", substr(pdbid, 2, 3),
                     "/pdb", pdbid, ".ent.gz")
      

    }
    else {
      file <- get.pdb(pdbid, path=configuration$pdbdir$rawfiles)
    }

    if(!file.exists(file))
      return(NULL)

    progress$set(value = 3)
    progress$set(message = 'Parsing PDB',
                 detail = 'Please wait')
    pdb <- read.pdb(file, verbose=FALSE)
    progress$set(value = 5)
    
    return(pdb)
  }
  else {
    return(NULL)
  }
}

## returns PDB code (4-characters)
get_pdbid <- reactive({
  if(is.null(rv$pdbid))
    return(NULL)
  else {
    pdbid <- rv$pdbid
    return(pdbid)
  }
})

get_pdbid6 <- reactive({
  return(paste(rv$pdbid, paste(rv$chains_selected, collapse=""), sep="_"))
})

## returns the final PDB object
raw_pdb <- reactive({
  pdbid <- rv$pdbid
  pdb <- read_pdb(pdbid)
  return(pdb)
})

## returns the final all-atom PDB object
final_pdb <- reactive({
  message("final_pdb")
  pdb <- rv$raw_pdb
  
  if(is.pdb(pdb)){
    if(is.vector(rv$chains_selected)) {
      if(length(rv$chains_selected)>0)
        pdb <- trim(pdb, chain=rv$chains_selected)
    }
  }

  ##rv$final_pdb <- pdb
  return(pdb)
})

## return all chain IDs
chain_pdb <- reactive({
  if(is.pdb(rv$raw_pdb)){
    chains <- unique(trim(rv$raw_pdb, "protein", elety="CA")$atom$chain)
    names(chains) <- chains
    rv$chains_all <- chains
    return(chains)
  }
  else {
    return(NULL)
  }
})

## checkbox 
output$chain_input <- renderUI({
  selectInput("chainid", "Limit to chain IDs:",
              choices = rv$chains_all,
              selected = rv$chains_selected,
              multiple=TRUE)
})


##init_show_pdbs <- TRUE
output$pdbWebGL  <- renderWebGL({
  pdb <- rv$final_pdb

  as <- input$view_inpdb_as
  bg <- "white"
  
  if("calpha" %in% as) {
    view.pdb(pdb, as="calpha", col=input$view_inpdb_col, bg.col=bg, sheet="blue")
  }

  if("overview" %in% as) {
    view.pdb(pdb, as="calpha", col=input$view_inpdb_col, lwd=5, bg.col=bg, sheet="blue")
    view.pdb(pdb, as="all", col="atom", add=TRUE)
  }
})


## prints a longer log of the input PDB
output$pdb_log <- renderPrint({
  pdbsum(rv$final_pdb, pdbid=rv$pdbid, chainid=rv$chains_selected)
})

#output$pdbSummary <- renderPrint({
#  input$pdbaction
#  invisible(capture.output( pdb <- final_pdb() ))
#  print(pdb)
#})
