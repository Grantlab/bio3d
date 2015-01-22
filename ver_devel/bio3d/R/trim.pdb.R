"trim.pdb" <-function(pdb, ..., inds=NULL, sse=TRUE) {
  if(!is.pdb(pdb))
    stop("input 'pdb' must be a PDB list object as returned from 'read.pdb'")

  cl <- match.call()
  extra.args <- list(...)
  
  if(length(extra.args)>0) {
     if(!is.null(inds))
        warning("Multiple atom selections are provided. Only the one indicated by 'inds' will be used")
     else if(is.select(extra.args[[1]])) 
        # to be back-compatible with the habit calling trim.pdb(pdb, inds)
        inds = extra.args[[1]]
     else
        inds = atom.select(pdb, ...)
  }

  if(is.null(inds))
    stop("no selection indices provided")

  if(!is.list(inds))
    stop("selection indices must be provided i.e. from 'atom.select'")

  if(is.null(inds$atom) || is.null(inds$xyz))
    stop("selection indices must be provided i.e. from 'atom.select'")

  ## Trim main components
  atom <- pdb$atom[inds$atom,]
  ca.inds <- atom.select(pdb, "calpha", verbose=FALSE)
  calpha <- inds$atom %in% ca.inds$atom
#  calpha = (atom[,"elety"]=="CA") & (atom[,"resid"]!="CA") & (atom[,"type"]=="ATOM")
  
  
  if(is.null(nrow(pdb$xyz))) {
    xyz <-  pdb$xyz[inds$xyz, drop=FALSE]
  } else {
    xyz <-  pdb$xyz[,inds$xyz, drop=FALSE]
  }
  
  helix <- NULL; sheet <- NULL;
  
  if(sse) {
    ss <- pdb2sse(pdb)

    ##- Trimed positions
    calpha2 <- ca.inds$atom %in% inds$atom
    ss <- ss[calpha2]
    
    ##- New sse
    new.sse <- bounds.sse(ss)

    helix <- new.sse$helix
    ##- add back other components
    add <- pdb$helix[!names(pdb$helix) %in% names(new.sse$helix)]
    ##- match sse number in case some sse are completely removed
    add <- lapply(add, function(x) x[new.sse$helix$id])
    helix <- c(helix, add)

    sheet <- new.sse$sheet
    ##- add back other components
    add <- pdb$sheet[!names(pdb$sheet) %in% names(new.sse$sheet)]
    ##- match sse number in case some sse are completely removed
    add <- lapply(add, function(x) x[new.sse$sheet$id])
    sheet <- c(sheet, add)

    ##- remove 'id'; Maybe we don't need it?
    helix$id <- NULL
    sheet$id <- NULL
  }

  calpha = (atom[,"elety"]=="CA") & (atom[,"resid"]!="CA") & (atom[,"type"]=="ATOM")
  output<-list(atom=atom,
               helix=helix, 
               sheet=sheet, 
               seqres=pdb$seqres, ## return unmodified
               xyz=xyz,
               calpha = calpha,
               call = cl)
  
  class(output) <- c("pdb", "sse")
  class(output$xyz) <- c("numeric","xyz")
  return(output)
}
