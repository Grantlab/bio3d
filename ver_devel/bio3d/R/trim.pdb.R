"trim.pdb" <-function(pdb, inds=NULL, sse=TRUE) {
  if(inherits(pdb, "pdb"))
    stop("input 'pdb' must be a PDB list object as returned from 'read.pdb'")

  if(is.null(inds))
    stop("no selection indices provided")

  if(!is.list(inds))
    stop("selection indices must be provided i.e. from 'atom.select'")

  if(is.null(inds$atom) || is.null(inds$xyz))
    stop("selection indices must be provided i.e. from 'atom.select'")

  ## Trim main components
  atom <- pdb$atom[inds$atom,]
  calpha <- as.logical(atom[,"elety"]=="CA")

  if(is.null(nrow(pdb$xyz))) {
    xyz <-  pdb$xyz[inds$xyz]
  } else {
    ## Multi-model records
    xyz <-  pdb$xyz[,inds$xyz]
  }

  helix <- NULL; sheet <- NULL;

  if(sse) {
    ##-- Build reference SSE sequence matrix 'ss'
    ref <- paste(pdb$atom[pdb$calpha, "resno"], 
                 pdb$atom[pdb$calpha, "chain"], sep="_")
    ss <- matrix(NA, ncol=length(ref), nrow=3) ## sse vector
    colnames(ss) = pdb$atom[pdb$calpha, "resno"]

    ##- Trimed positions
    ref.trim <- paste(atom[calpha, "resno"], 
                      atom[calpha, "chain"], sep="_")

    if(length(pdb$helix$start)>0) {
      ##- HELIX chain, type and resno
      chain.h <- rep(pdb$helix$chain, (pdb$helix$end-pdb$helix$start+1))
      type.h <- rep(pdb$helix$type, (pdb$helix$end-pdb$helix$start+1))
      resno.h <- unbound(pdb$helix$start, pdb$helix$end)

      ##- Use reference vector for filling 'ss'
      ref.h <- paste(resno.h, chain.h, sep="_")
      ss[1,ref %in% ref.h]="H"
      ss[2,ref %in% ref.h]=chain.h
      ss[3,ref %in% ref.h]=type.h
    }


    if(length(pdb$sheet$start)>0) {
      ##- SHEET chain, type and resno
      chain.e <- rep(pdb$sheet$chain, (pdb$sheet$end-pdb$sheet$start+1))
      type.e <- rep(pdb$sheet$sense, (pdb$sheet$end-pdb$sheet$start+1))
      resno.e <- unbound(pdb$sheet$start, pdb$sheet$end)

      ##- Use reference vector for filling 'ss'
      ref.e <- paste(resno.e, chain.e, sep="_")
      ss[1,ref %in% ref.e]="E"
      ss[2,ref %in% ref.e]=chain.e
      ss[3,ref %in% ref.e]=type.e
    }

    ##- Lookup trimed positions 'ref.trim' in 'ss'
    s <- ss[,ref %in% ref.trim, drop=FALSE]

    if( any(s[1,] =="H",na.rm=TRUE) ) {
      h.bounds <- bounds( as.numeric(names(s[1,][s[1,]=="H"])) )
      helix$start = h.bounds[,"start"]
      helix$end = h.bounds[,"end"]
      helix$chain = ss[2,as.character(h.bounds[,"start"])]
      helix$type = ss[3,as.character(h.bounds[,"start"])]
    } #else {
      #hexix$start=NULL; helix$end=NULL
      #helix$chain=NULL; helix$type=NULL
    #}

    if( any(s[1,] =="E",na.rm=TRUE) ) {
      e.bounds <- bounds( as.numeric(names(s[1,][s[1,]=="E"])) )
      sheet$start = e.bounds[,"start"]
      sheet$end = e.bounds[,"end"]
      sheet$chain = ss[2,as.character(e.bounds[,"start"])]
      sheet$sense = ss[3,as.character(e.bounds[,"start"])]
    } #else {
      #sheet$start=NULL; sheet$end=NULL
      #sheet$chain=NULL; sheet$sense=NULL
    #}
  }

  output<-list(atom=atom,
               helix=helix, 
               sheet=sheet, 
               seqres=pdb$seqres, ## return unmodified
               xyz=xyz,
               calpha = calpha)
  
  class(output) <- c("pdb", "sse")
  class(output$xyz) <- c("numeric","xyz")
  return(output)
}
