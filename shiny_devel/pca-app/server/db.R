## 

db_con <- function() {
  con <- dbConnect(RMySQL::MySQL(), dbname = "test", username="root",
                   password="", host="localhost")
  return(con)
}

db_discon <- function(con) 
  dbDisconnect(con)

db_get_anno <- function(acc) {
  

  missing_inds <- c()
  
  con <- db_con()
  for(i in 1:length(acc)) {
    query <- paste0("SELECT acc FROM pdb_annotation WHERE acc='", toupper(acc[i]), "'")

    res <- dbGetQuery(con, query)
    
    if(!nrow(res)>0)
      missing_inds <- c(missing_inds, i)
  }
  
  if(length(missing_inds)>0)
    res <- db_add_anno(acc[missing_inds], con=con)
  
  where <- paste0("WHERE acc IN ('", paste(acc, collapse="', '"), "')")
  query <- paste("SELECT acc, compound, source, ligandId, chainLength FROM pdb_annotation", where)

  res <- dbGetQuery(con, query)
   
  db_discon(con)
  return(res)

}

db_add_anno <- function(acc, con=NULL) {

  message(acc)
  anno <- pdb.annotate(acc)
  anno$acc <- rownames(anno)
  message(anno$acc)

  close <- FALSE
  if(is.null(con)) {
    con <- db_con()
    close <- TRUE
  }
  
  res <- dbWriteTable(con, "pdb_annotation",
                      
                      value=anno[, c("acc", "structureId", "chainId", "source", "compound",
                        "experimentalTechnique", "resolution", "ligandId", "ligandName",
                        "chainLength", "db_id")], 
                      
                      overwrite=FALSE,
                      row.names=FALSE,
                      append=TRUE)

  if(close)
    db_discon(con)
  
  return(res)
}
