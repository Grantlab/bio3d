pdb.pfam <- function(ids, best.only=TRUE, compact=TRUE) {

  ##-- Annotate PDBs with PFAM hits from RCSB PDB
  ##    Input 'pdbid' can be a vector of ids with chainID:
  ##  
  ##     pdb.pfam(c("1bg2", "2kin", "4q21", "3kin_A"))
  ## 
  ##   There are often multiple hits to a structure. The option
  ##    'best.only' will keep the lowest eValue one only.
  ##

  oops <- requireNamespace("XML", quietly = TRUE)
  if(!oops) { stop("Please install the XML package from CRAN") }

  ##- 'ids' may optionally contain chain ID (e.g '3kin_A')
  pdbid <- ids    ## keep a copy of input 'ids' for later reporting
  ids   <- toupper(ids) 
  chains <- rep(NA, length(ids))

  standard.ids <- nchar(ids) == 4
  if(any(!standard.ids)) {
    ids[!standard.ids] <- sub("\\.pdb$","", basename(ids[!standard.ids]))
    chains[!standard.ids] <- substr(ids[!standard.ids],6,6)
    ids[!standard.ids] <- substr(ids[!standard.ids],1,4) 
  }

	url <- paste0("http://www.rcsb.org/pdb/rest/hmmer?structureId=", 
                paste0(ids, collapse=","))
	xl <-  try(XML::xmlToList(url), silent=TRUE)

	##- Exit with NULL if PFAM annotation page not found (cant load)
	if(!is.list(xl)) { return(NULL) }

	## Else process and order data.frame
	df <- data.frame(matrix(unlist(xl), ncol=8, byrow=TRUE), stringsAsFactors=FALSE)
	colnames(df) <- names(xl[1]$pfamHit)
  df <- df[order(df$eValue),]

  ## Add ID in from of 'pdbID_chainID' and PFAMname_PFAMacc
  df$ID = paste(df$structureId, df$chainId,sep="_")
  df$PFAM = paste0(df$pfamName," (",df$pfamAcc,")")

  ##- Optionally keep top hit for each structure only 
  if(best.only) { df <- df[!duplicated(df$ID),] }

  ##- If chain Id was found in input keep only entries for that chain 
  ind <- !is.na(chains)
  if(any(ind)) {
    ## Exclude entries that match PDBID but not PDBID_chainID
    filter.id1 <- ids[ind]
    filter.id2 <- paste(ids[ind], chains[ind],sep="_")
    filter.ind1 <- which(df$structureId %in% filter.id1)
    filter.ind2 <- which(df$ID %in% filter.id2)
    filter.ind  <- filter.ind1[!filter.ind1 %in% filter.ind2]
    if(length(filter.ind) > 0) { df <- df[-filter.ind,] }
  }

  ##- Check if any 'pdbids' were not found and add a blank row to df
  missing <- !toupper(ids) %in% df$structureId
  if( any(missing) ){
    warning(paste("PFAM annotation not found for pdbid(s):", paste(pdbid[missing], collapse=", ")))
    blank <- matrix(rep(NA, ncol(df)*sum(missing)), nrow=sum(missing))
    colnames(blank) <- colnames(df)
    blank[,"ID"] <- toupper(pdbid[missing])
    df <- rbind(df,blank)
  }

  ## Optionally rtn a subset of data cols only
  if(compact) {
    return(df[,c("ID","PFAM","pfamDesc","eValue")])
  } else { return(df) }

}

