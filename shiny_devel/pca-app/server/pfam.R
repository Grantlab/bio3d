
## returns PFAM annotation for input PDB
get_pfam_annotation <- reactive({
  message("get_pfam_anotation called")
  
  if(input$input_type == "multipdb") {
    ids <- get_multipdbids()
    anno <- get_pfam(ids)
    return(anno[ anno$acc %in% ids, ])
  }
  else {
    pdbid <- get_pdbid()
    return(get_pfam(pdbid))
  }

})


##-- PFAM annotation of single or multiple PDBs
output$pfam_table <- DT::renderDataTable({
  message("input_pdb_pfam called")
  pdbid <- get_pdbid()
  #chainid <- get_chainid()
  if(is.null(pdbid)) {
    return()
  }
  get_pfam_table()
})

output$pfam_table_multi <- DT::renderDataTable({
  message("pdb_pfam_multi called")
  if(is.null(input$pdb_codes)) {
    return()
  }
  get_pfam_table()

})

get_pfam_table <- reactive({
  pfam <- get_pfam_annotation()

  if(!nrow(pfam) > 0) {
    return(DT::datatable(data.frame("Pfam data not found"),
                         class = 'compact row-border',
                         selection = "none",
                         rownames = FALSE, colnames = FALSE,
                         options = list( dom = "t", autoWidth = TRUE)
                         )
           )
  }
  
  pfam$ID <- paste(pfam$structureId, pfam$chainId, sep = "_")
  pfam$PFAM <- paste0(pfam$pfamName, " (", pfam$pfamAcc, ")")
  pfam <- pfam[, c("ID", "PFAM", "pfamDesc", "eValue")]
  
  pfam['PFAM'] <- lapply(pfam['PFAM'], function(x) {
      pfid <- regmatches(x, gregexpr("(?<=\\().*?(?=\\))", x, perl=T))[[1]]
      link <- gsub(pfid,
                   tags$a(href=paste0("http://pfam.xfam.org/family/",pfid), target="_blank", pfid),
                   x)
  })
  colnames(pfam)=c("ID", "PFAM", "Annotation","eValue")
  DT::datatable(pfam, escape = FALSE,
                class = 'compact row-border',
                selection = "none", rownames = FALSE,
                options = list( dom = "t", autoWidth = TRUE,
                    columnDefs = list( 
                      list( orderable = 'false', targets = c(0,1,2) )
        ),
    initComplete = JS(
    "function(settings, json) {",
    '$(this.api().table().header()).find("th").removeClass("sorting");',
    '$(this.api().table().header()).find("th").prop("onclick",null).off("click");',
    "}")
    )
  )
})
##To Do:
##       Add this to multiple PDB IDs well/div also.

