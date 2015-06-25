
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

get_annotation <- function(acc, use_chain=TRUE) {
  acc <- format_pdbids(acc)
  unq <- unique(substr(acc, 1, 4))
  
  if(!use_chain) {
    acc <- unq
  }

  con <- db_connect()

  if(!inherits(con, "MySQLConnection")) {
    warning("could not connect to database :(")
    anno <- pdb.annotate(acc, anno.terms=c(
                                "structureId", "chainId",
                                "source", "compound",
                                "experimentalTechnique", "resolution",
                                "ligandId", "ligandName",
                                "chainLength", "sequence"
                                )
                         )
    anno$acc <- rownames(anno)
    return(anno)
  }

  missing_inds <- c()
  for(i in 1:length(unq)) {
    query <- paste0("SELECT acc FROM pdb_annotation WHERE structureId='", unq[i], "'")
    res <- dbGetQuery(con, query)

    if(!nrow(res)>0)
      missing_inds <- c(missing_inds, i)
  }

  if(length(missing_inds) > 0)
    res <- db_add_annotation(unq[missing_inds], con=con)

  if(use_chain) {
    where <- paste0("WHERE acc IN ('", paste(acc, collapse="', '"), "')")
  }
  else {
    where <- paste0("WHERE structureId IN ('", paste(unq, collapse="', '"), "')")
  }
  query <- paste("SELECT acc, structureId, chainId, compound, source, ligandId, chainLength, sequence FROM pdb_annotation", where, "ORDER BY acc")

  res <- dbGetQuery(con, query)

  if(use_chain) {
    rownames(res) <- res$acc
    res <- res[acc, ]
  }

  db_disconnect(con)
  return(res)
}

db_add_annotation <- function(acc, con=NULL) {
  unq <- unique(substr(acc, 1, 4))
  anno <- pdb.annotate(unq)
  anno$acc <- row.names(anno)
  
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
                        "chainLength", "db_id", "sequence")],
                      overwrite=FALSE,
                      row.names=FALSE,
                      append=TRUE)

  if(close)
    db_disconnect(con)

  return(res)
}
