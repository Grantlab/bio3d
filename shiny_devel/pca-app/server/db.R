
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


get_annotation_from_pdb <- function(acc) {
  anno <- pdb.annotate(acc, anno.terms=c(
                              "structureId", "chainId",
                              "source", "compound",
                              "experimentalTechnique", "resolution",
                              "ligandId", "ligandName",
                              "chainLength", "db_id", "sequence"
                              )
                       )
  anno$acc <- rownames(anno)
  return(anno)
}

get_annotation <- function(acc, use_chain=TRUE) {
  acc <- format_pdbids(acc)
  unq <- unique(substr(acc, 1, 4))

  if(all(nchar(acc)==4))
    use_chain <- FALSE
  
  if(!use_chain) {
    acc <- unq
  }

  con <- db_connect()

  if(!inherits(con, "MySQLConnection")) {
    warning("could not connect to database :(")
    return(get_annotation_from_pdb(acc))
  }

  if(use_chain) {
    where <- paste0("WHERE acc IN ('", paste(acc, collapse="', '"), "')")
  }
  else {
    where <- paste0("WHERE structureId IN ('", paste(unq, collapse="', '"), "')")
  }
  
  query <- paste("SELECT acc, structureId, chainId, compound, source, ligandId, ligandName, chainLength, experimentalTechnique, resolution, sequence FROM pdb_annotation", where, "ORDER BY acc")

  res <- dbGetQuery(con, query)
 
  if(nrow(res) > 0) {
    if(use_chain) {
      rownames(res) <- res$acc
      res <- res[acc, ]
    }

    missing_inds <- !(unq %in% unique(res$structureId))
  }
  else {
    missing_inds <- rep(TRUE, length(unq))
  }

  if(any(missing_inds)) {

    print("missing in db:")
    print(unq[missing_inds])
    
    ## add missing to SQL db
    res <- db_add_annotation(unq[missing_inds], con=con)

    ## fetch again from SQL db
    res <- dbGetQuery(con, query)
  }

  db_disconnect(con)
  if(use_chain) {
    inds <- res$acc %in% acc
    return(res[inds,, drop=FALSE])
  }
  else {
    inds <- res$structureId %in% unq
    return(res[inds,, drop=FALSE])
  }

}

db_add_annotation <- function(acc, con=NULL) {
  unq <- unique(substr(acc, 1, 4))
  anno <- get_annotation_from_pdb(unq)
  ##print(head(anno))
  
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



get_pfam <- function(acc) {
  acc <- format_pdbids(acc)
  unq <- unique(substr(acc, 1, 4))
  
  con <- db_connect()
  
  if(!inherits(con, "MySQLConnection")) {
    warning("could not connect to database :(")
    pfam <- pdb.pfam(unq, compact = FALSE)
    pfam$acc <- paste(pfam$structureId, pfam$chainId, sep = "_")
    return(pfam)
  }

  missing_inds <- c()
  for(i in 1:length(unq)) {
    query <- paste0("SELECT acc FROM pfam WHERE structureId='", unq[i], "'")
    res <- dbGetQuery(con, query)

    if(!nrow(res)>0)
      missing_inds <- c(missing_inds, i)
  }

  if(length(missing_inds) > 0)
    res <- db_add_pfam(unq[missing_inds], con=con)
  
  #if(use_chain) {
  #  where <- paste0("WHERE acc IN ('", paste(acc, collapse="', '"), "')")
  #}
  #else {
  #  where <- paste0("WHERE structureId IN ('", paste(unq, collapse="', '"), "')")
  #}

  where <- paste0("WHERE structureId IN ('", paste(unq, collapse="', '"), "')")
  query <- paste("SELECT acc, structureId, chainId, pfamAcc, pfamName, pfamDesc, eValue FROM pfam", where, "ORDER BY acc")

  res <- dbGetQuery(con, query)

  #if(use_chain) {
  #  ##rownames(res) <- res$acc
  #  res <- res[acc, ]
  #}

  db_disconnect(con)
  return(res)
}


db_add_pfam <- function(acc, con=NULL) {
  unq <- unique(substr(acc, 1, 4))
  pfam <- pdb.pfam(unq, compact=FALSE)
  pfam$acc <- paste(pfam$structureId, pfam$chainId, sep = "_")
  
  close <- FALSE
  if(is.null(con)) {
    con <- db_connect()
    close <- TRUE
  }
  
  res <- dbWriteTable(con, "pfam",
                      value=pfam[, c(
                        "acc", "structureId", "chainId", 
                        "pdbResNumStart", "pdbResNumEnd",
                        "pfamAcc", "pfamName", "pfamDesc", "eValue")],
                      overwrite=FALSE,
                      row.names=FALSE,
                      append=TRUE)

  if(close)
    db_disconnect(con)

  return(res)
}
