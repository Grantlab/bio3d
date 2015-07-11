
##get_blasttable <- reactive({
output$blast_table <- renderDataTable({
  message("blasttable1")
  acc <- rv$blast

  if(is.null(acc))
    return()
  
  anno <- get_annotation(acc)
  anno$struct_nr <- 1:nrow(anno)

  anno$acc <- paste0("<a href=\"", "http://pdb.org/pdb/explore/explore.do?structureId=",
                     substr(anno$acc, 1, 4), "\" target=\"_blank\">", anno$acc, "</a>")
  anno$compound <- substr(anno$compound, 1, 20)

  
  ## http://xray.bmc.uu.se/hicup/ATP/
  anno$ligandId <- unlist(lapply(anno$ligandId,
         function(x) {
           if(is.na(x))
             return("")

           spl <- unlist(strsplit(x, ","))
           if(length(spl) > 0) {
             return(paste(
                         sapply(spl, function(x) { as.character(tags$a(href = paste0("http://www.rcsb.org/pdb/ligand/ligandsummary.do?hetId=", x), target = '_blank', x))}),
                          collapse=", "))
           }
           else {
             return("")
           }

         }))
  
  rownames(anno) <- NULL
  show.cols <- c("acc", "compound", "ligandId")
  col.inds <- sapply(show.cols, grep, colnames(anno))
  
  anno <- anno[, col.inds]
  
  datatable(anno, extensions = 'Scroller', escape = FALSE,
            colnames = c("ID", "Name", "Ligand ID"),
            options = list(
              deferRender = TRUE,
              dom = "frtiS",
              scrollY = 200,
              scrollCollapse = TRUE
              ))
})



