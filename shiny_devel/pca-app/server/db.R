
db_connect <- function() {
  db <- configuration$db
  if(!db$use)
    return(FALSE)
  
  con <- tryCatch(
    {
      dbConnect(RMySQL::MySQL(),
                dbname = db$dbname, username = db$username,
                password = db$password, host = db$host)
    },
    error = function(e) {
      con <- FALSE
    }
    )
  
  return(con)
}

db_disconnect <- function(con) 
  dbDisconnect(con)

get_annotation <- function(acc) {
  con <- db_connect()
  
  if(!inherits(con, "MySQLConnection")) {
    warning("could not connect to database :(")
    anno <- pdb.annotate(acc, anno.terms=c(
                                "structureId", "chainId",
                                "source", "compound",
                                "experimentalTechnique", "resolution",
                                "ligandId", "ligandName",
                                "chainLength"
                                )
                         )
    anno <- cbind(acc, anno)
    return(anno)
  }
  
  missing_inds <- c()
  for(i in 1:length(acc)) {
    query <- paste0("SELECT acc FROM pdb_annotation WHERE acc='", toupper(acc[i]), "'")

    res <- dbGetQuery(con, query)
    
    if(!nrow(res)>0)
      missing_inds <- c(missing_inds, i)
  }
  
  if(length(missing_inds)>0)
    res <- db_add_annotation(acc[missing_inds], con=con)
  
  where <- paste0("WHERE acc IN ('", paste(acc, collapse="', '"), "')")
  query <- paste("SELECT acc, compound, source, ligandId, chainLength FROM pdb_annotation", where)

  res <- dbGetQuery(con, query)
   
  db_disconnect(con)
  return(res)

}

db_add_annotation <- function(acc, con=NULL) {
  anno <- pdb.annotate(acc)
  anno$acc <- rownames(anno)

  close <- FALSE
  if(is.null(con)) {
    con <- db_connect()
    close <- TRUE
  }
  
  res <- dbWriteTable(con, "pdb_annotation",
                      value=anno[, c("acc", "structureId",
                        "chainId", "source", "compound",
                        "experimentalTechnique", "resolution",
                        "ligandId", "ligandName",
                        "chainLength", "db_id")],
                      overwrite=FALSE,
                      row.names=FALSE,
                      append=TRUE)
  
  if(close)
    db_disconnect(con)
  
  return(res)
}
