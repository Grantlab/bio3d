##- PDB input UI that is responsive to 'reset button' below
output$resetable_pdb_input <- renderUI({
  ## 'input$reset_pdb_input' is just used as a trigger for reset
  reset <- input$reset_pdb_input
  textInput("pdbid", label="Enter RCSB PDB code/ID:", value = "4Q21") #)
})

final_sequence <- function() {
  ## option 1 - PDB code provided
  if(input$input_type == "pdb") {
    anno <<- pdb_anno()
  
    if(is.vector(input$chainId)) {
      ind <- which(anno$chainId==input$chainId[1])
      seq <- unlist(strsplit(anno$sequence[ind], ""))
    }
    else {
      seq <- unlist(strsplit(anno$sequence[1], ""))
    }
    return(as.fasta(seq))
  }

  ## option 2 - sequence provided
  if(input$input_type == "sequence") {
    go <- input$action_input

    inp <- unlist(strsplit(input$sequence, "\n"))
    inds <- grep("^>", inp, invert=TRUE)

    if(!length(inds)>0)
      stop("Error reading input sequence")

    inp <- toupper(paste(inp[inds], collapse=""))

    print(inp)
    
    seq <- as.fasta(unlist(strsplit(inp, "")))
    print(seq)
    return(seq)
  }
    
}


pdb_anno <- reactive({
  print(input$pdbid)
  
  if(nchar(input$pdbid)==4) {
    progress <- shiny::Progress$new(session, min=1, max=5)
    on.exit(progress$close())
    
    progress$set(message = 'Fetching PDB data',
                 detail = 'This should be quick...')
    progress$set(value = 3)

    anno <- pdb.annotate(input$pdbid)
    for(i in 4:5)
      progress$set(value = i)
    
    return(anno)
  }
  else {
    stop("Provide a 4 character PDB code")
  }
})


## all available chains in PDB
output$chains1 <- renderPrint({
  invisible(capture.output(anno <- pdb_anno()))
  cat( anno$chainId, sep=", ")
})

## checkbox 
output$chains2 <- renderUI({
  invisible(capture.output(anno <- pdb_anno()))
  chains <- anno$chainId
  radioButtons("chainId", label="Choose chain ID:",
               choices=chains, inline=TRUE)
  
})
