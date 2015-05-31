
############################
## Overlap analysis
###########################

get_pdb2 <- function(acc) {
  ## downloads and reads the raw PDB
  id <- substr(acc, 1, 4)
  
  if(nchar(id)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Fetching PDB',
                 detail = 'This should be quick...')
    progress$set(value = 2)
    
    if(configuration$pdbdir$archive) {
      id2 <- paste0(tolower(substr(id, 1, 4)), "_", substr(id, 6, 6))
      raw.files <- paste0(configuration$pdbdir$rawfiles, "/", substr(id2, 2, 3), "/pdb", ids, ".ent.gz")
      
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
}

## returns the final PDB object
#get_pdb2 <- reactive({
#  pdb <- raw_pdb2()
#  
#  if(is.vector(input$chains2)) {
#    pdb <- trim(pdb, chain=input$chains2)
#  }
#  
#  return(pdb)
#})



## return chain IDs
#chain_pdb2 <- reactive({
#  invisible(capture.output( pdb <- raw_pdb2() ))
#  chains <- unique(trim(pdb, "protein", elety="CA")$atom$chain)
#  names(chains) <- chains
#  return(chains)
#})


## all available chains in PDB
#output$chains3 <- renderPrint({
#  invisible(capture.output(  chains <- chain_pdb2() ))
#  cat( chains, sep=", ")
#})

## checkbox 
#output$chains4 <- renderUI({
#  chains <- chain_pdb2()
#  checkboxGroupInput("chains2", label = "Limit to chain IDs:", 
#                     choices = chains, inline = TRUE)
#})


calc_overlap <- reactive({
  input$goButton
  
  ## get list of blast hits
  hits <- filter_hits()

  if(length(input$blast_table_rows_selected)>0) {
    ids <- hits[input$blast_table_rows_selected]
  }
  else {
    stop("no pDBs selected")
  }

  pdblist <- vector("list", length(ids)+1)
  pdblist[[1]] <- get_pdb()

  for(i in 1:length(ids)) {
    pdblist[[i+1]] <- get_pdb2(ids[i])
  }

  for(i in 1:length(pdblist)) {
    pdblist[[i]] <- trim(pdblist[[i]], "calpha")
  }

  ## pdb for nma
  pdb <- pdblist[[1]]

  pdbs <- pdbaln(pdblist)
  pdbs$id[1] <- input$pdbid
  pdbs$id[2:length(pdbs$id)] <- ids
  pdbs$xyz <- pdbfit(pdbs)

  gaps.pos <- gap.inspect(pdbs$xyz)
  gaps.res <- gap.inspect(pdbs$ali)

  ## indices for nma calc
  tomatch <- paste(pdbs$resno[1, gaps.res$f.inds],
                   pdbs$chain[1, gaps.res$f.inds], sep="|")
  resids <-  paste(pdb$atom$resno, pdb$atom$chain, sep="|")
  ##print(resids %in% tomatch)
  nmainds <- as.select(resids %in% tomatch)
  
  if(input$forcefield %in% c("calpha", "sdenm", "reach"))
    modes <- nma(pdb, ff=input$forcefield, mass=FALSE, temp=300,
                 outmodes = nmainds)
  
  if(input$forcefield %in% c("anm", "pfanm"))
    modes <- nma(pdb, ff=input$forcefield, cutoff=input$cutoff,
                 mass=FALSE, temp=NULL, outmodes = nmainds)

  print(pdbs)
  print(modes)
  
  overlaps <- vector("list", length(pdbs$id)-1)
  for(i in 2:length(pdbs$id)) {
    inds <- c(1, i)
    dv <- difference.vector(pdbs$xyz[inds, gaps.pos$f.inds, drop=FALSE])
    ov <- overlap(modes=modes, dv=dv)
    ov$id <- pdbs$id[i]
    overlaps[[i-1]] <- ov
  }
  return(overlaps)
})


output$overlap_plot <- renderPlot({
  overlap_plot2()
})

overlap_plot2 <- function() {
  o <- calc_overlap()
  ids <- sapply(o, function(x) x$id)
  
  #plot.new()
  #plot.window(xlim = c(0,20), ylim = c(0,1))

  plot(o[[1]]$overlap, type='h', ylim=c(0,1),
       ylab="Overlap", xlab="Mode index", main="Overlap analysis")
  points(o[[1]]$overlap)
  
  for(i in 1:length(o)) {
    ov <- o[[i]]
    lines(ov$overlap.cum, type='b', col=i+1)
  }

  legend("bottomright", ids, col=1:length(ids)+1, lty=1)
  
  ##plot(o$overlap, type='h', ylim=c(0,1), ylab="Overlap", xlab="Mode index")
  ##points(o$overlap)
  ##lines(o$overlap.cum, type='b', col='red')
}

