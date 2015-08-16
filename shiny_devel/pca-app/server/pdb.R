
## returns PDB code (4-characters)
get_pdbid <- reactive({
  return(rv$pdbid)
})

## returns the selected chain ID
get_chainid <- reactive({
  if (is.null(rv$chainid))
    return()
  else
    return(rv$chainid)
})

## returns all chain IDs in PDB
get_chainids <- reactive({
  pdbid <- get_pdbid()

  if(is.null(pdbid))
    return()

  anno <- input_pdb_annotation()
  return(anno$chainId)
})

## parse and return pdb codes for multipdb input type
get_multipdbids <- reactive({
  acc <- unique(trim(unlist(strsplit(rv$pdb_codes, ","))))
  inds <- nchar(acc) > 3 & nchar(acc) < 7
  acc <- acc[inds]
  acc <- format_pdbids(acc)
  
  anno <- get_annotation(acc, use_chain=FALSE)
  inds <- unlist(sapply(acc, grep, anno$acc))
  anno <- anno[inds, ]
  acc <- anno$acc
  
  return(acc)
})

## returns the PDB object (trimmed to chain ID)
get_pdb <- reactive({
  message("get_pdb called")
  pdbid <- get_pdbid()
  chainid <- get_chainid()

  anno <- input_pdb_annotation()

  if(!nrow(anno) > 0)
    stop("no annotation data found")

  if(is.vector(rv$chainId)) {
    ind <- which(anno$chainId == rv$chainid)
    anno <- anno[ind,]
  }
  else {
    anno <- anno[1,]
  }

  id <- anno$acc
  if(configuration$pdbdir$archive) {
    ids <- paste0(tolower(substr(id, 1, 4)), "_", substr(id, 6, 6))
    raw.files <- paste0(configuration$pdbdir$splitfiles, "/", substr(ids, 2, 3), "/pdb", ids, ".ent.gz")

    if(!file.exists(raw.files))
      stop("PDB not found")

    pdb <- read.pdb(raw.files)
  }
  else {
    raw.files <- get.pdb(unique(substr(id, 1,4)), path=configuration$pdbdir$rawfiles, gzip=TRUE)
    pdb <- read.pdb(raw.files)
    pdb <- trim.pdb(pdb, chain=substr(id, 6,6))
  }

  return(pdb)
})

## returns annotation data for input PDB
input_pdb_annotation <- reactive({
  message("input_pdb_annotation called")

  pdbid <- get_pdbid()

  if(nchar(pdbid)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())

    progress$set(message = 'Fetching PDB data',
                 detail = 'Please wait')
    progress$set(value = 3)

    anno <- get_annotation(pdbid, use_chain=FALSE)
    progress$set(value = 5)

    return(anno)
  }
  else {
    stop("Provide a 4 character PDB code")
  }
})

## prints a short summary of the input PDB
output$input_pdb_summary <- renderPrint({
  message("input_pdb_summary called")
  pdbid <- get_pdbid()
  chainid <- get_chainid()
  if(is.null(pdbid)) {
    return()
  }

  input$input_type
  anno <- input_pdb_annotation()

  if(is.vector(chainid)) {
    ind <- which(anno$chainId==chainid[1])
    anno <- anno[ind,]
  }
  else {
    anno <- anno[1,]
  }

  cat(anno$compound[1], "\n",
      "(", anno$source[1], ")")

})



## prints a longer log of the input PDB
output$pdb_log <- renderPrint({
  pdbid <- get_pdbid()
  chainid <- get_chainid()
  invisible(capture.output( pdb <- get_pdb() ))

  pdbsum(pdb, pdbid=pdbid, chainid=chainid)
    
})

##
##-- Structure Viewer tab
##
output$pdbWebGL  <- renderWebGL({
  pdb <- get_pdb()

  as <- input$view_inpdb_as
  bg <- "white"
  
  if("calpha" %in% as) {
    view.pdb(pdb, as="calpha", col=input$view_inpdb_col, bg.col=bg, sheet="blue")
  }
  
  if("overview" %in% as) {
    view.pdb(pdb, as="calpha", col=input$view_inpdb_col,
             bg.col=bg, sheet="blue")
    view.pdb(pdb, as="ligand", col="atom", add=TRUE)
  }
  
  if("allatoms" %in% as) {
    view.pdb(pdb, as="calpha", col=input$view_inpdb_col,
             lwd=5, bg.col=bg, sheet="blue")
    view.pdb(pdb, as="ligand", col="atom", add=TRUE)
    view.pdb(pdb, as="protein", col="atom", add=TRUE)
  }
  
})

output$pdb_chains <- renderUI({
  chains <- get_chainids()

  selectInput("chainId", "Limit to chain id:",
              choices = chains, selected = chains[1], multiple=FALSE)
})
