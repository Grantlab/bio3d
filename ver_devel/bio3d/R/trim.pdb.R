"trim.pdb" <-function(pdb, inds=NULL, sse=TRUE) {
  if(!is.pdb(pdb))
    stop("input 'pdb' must be a PDB list object as returned from 'read.pdb'")

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
    ##-- Build reference SSE sequence matrix 'ss'
    ref <- paste(pdb$atom[ca.inds$atom, "resno"], 
                 pdb$atom[ca.inds$atom, "chain"], sep="_")
    ss <- matrix(NA, ncol=length(ref), nrow=4) ## sse reference matrix
    ss[4,] <- pdb$atom[ca.inds$atom, "resno"]    ## resno
    colnames(ss) <- ref

    ##- Trimed positions
    ref.trim <- paste(atom[calpha, "resno"], 
                      atom[calpha, "chain"], sep="_")

    if(length(pdb$helix$start)>0) {
      ##- HELIX chain, type and resno
      chain.h <- rep(pdb$helix$chain, (pdb$helix$end-pdb$helix$start+1))
      type.h <- rep(pdb$helix$type, (pdb$helix$end-pdb$helix$start+1))
      resno.h <- unbound(pdb$helix$start, pdb$helix$end)

      ## catch for missing chain identifier
      chain.h[chain.h == " "] = NA

      ##- Use reference vector for filling 'ss'
      ref.h <- paste(resno.h, chain.h, sep="_")
      i <- ref.h %in% ref ## Fix for HETATOM in SSE (1P3Q)
      ss[1, ref.h[i]]="H"
      ss[2, ref.h[i]]=chain.h[i]
      ss[3, ref.h[i]]=type.h[i]
    }


    if(length(pdb$sheet$start)>0) {
      ##- SHEET chain, type and resno
      chain.e <- rep(pdb$sheet$chain, (pdb$sheet$end-pdb$sheet$start+1))
      if(is.null(pdb$sheet$sense))
         type.e <- rep("", sum(pdb$sheet$end-pdb$sheet$start+1))
      else 
         type.e <- rep(pdb$sheet$sense, (pdb$sheet$end-pdb$sheet$start+1))
      resno.e <- unbound(pdb$sheet$start, pdb$sheet$end)

      ## catch for missing chain identifier
      chain.e[chain.e == " "] = NA

      ##- Use reference vector for filling 'ss'
      ref.e <- paste(resno.e, chain.e, sep="_")
      i <- ref.e %in% ref
      ss[1, ref.e[i]]="E"
      ss[2, ref.e[i]]=chain.e[i]
      if(!is.null(type.e))
        ss[3, ref.e[i]]=type.e[i]
    }

    ##- Lookup trimed positions 'ref.trim' in 'ss'
    s <- ss[, ref.trim, drop=FALSE]

    if( any(s[1,] =="H",na.rm=TRUE) ) {
      sh <- (s[, s[1,] %in% "H"])
      h.bounds <- bounds( as.numeric(sh[4,]), pre.sort=FALSE)
      h.inds <- rle2(rep(1:nrow(h.bounds), h.bounds[,"length"]))$inds

      helix$start = h.bounds[,"start"]
      helix$end = h.bounds[,"end"]

      helix$chain = sh[2, h.inds]
      helix$type = sh[3, h.inds]
    }

    if( any(s[1,] =="E",na.rm=TRUE) ) {
      se <- (s[, s[1,] %in% "E"])
      e.bounds <- bounds( as.numeric(se[4,]), pre.sort=FALSE)
      e.inds <- rle2(rep(1:nrow(e.bounds), e.bounds[,"length"]))$inds

      sheet$start = e.bounds[,"start"]
      sheet$end = e.bounds[,"end"]

      sheet$chain = se[2, e.inds]
      sheet$sense = se[3, e.inds]
    } 
  }

  ## Filtering the 'connectivity' component
  if(!is.null(pdb$con)) {
    r <-     pdb$con$eleno.1 %in% atom$eleno
    r <- r & pdb$con$eleno.2 %in% atom$eleno
    if(any(r))
      pdb$con <- pdb$con[r, , drop = FALSE]
    else
      pdb["con"] <- list(NULL)
  }

  calpha = (atom[,"elety"]=="CA") & (atom[,"resid"]!="CA") & (atom[,"type"]=="ATOM")
  output<-list(atom=atom,
               helix=helix, 
               sheet=sheet, 
               seqres=pdb$seqres, ## return unmodified
               xyz=xyz,
               con=pdb$con,
               calpha = calpha)
  
  class(output) <- c("pdb", "sse")
  class(output$xyz) <- c("numeric","xyz")

  return(output)
}
