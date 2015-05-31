
get_blasttable <- reactive({
  blast <- run_blast()
  acc <- filter_hits()
  
  anno <- get_annotation(acc)
  anno$struct_nr <- 1:nrow(anno)

  anno$acc <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=",
                     substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  anno$compound <- substr(anno$compound, 1, 20)
  
  rownames(anno) <- NULL
  show.cols <- c("acc", "compound")
  col.inds <- sapply(show.cols, grep, colnames(anno))
  
  print(col.inds)
  print(head(anno[, col.inds]))
  return(anno[, col.inds])
})

output$blast_table <- renderDataTable({
  datatable(get_blasttable(), extensions = 'Scroller', escape = FALSE,
            colnames = c("ID", "Name"),
            options = list(
              deferRender = TRUE,
              dom = "frtiS",
              scrollY = 200,
              scrollCollapse = TRUE
              ))
})



