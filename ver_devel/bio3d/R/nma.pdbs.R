"nma.pdbs" <- function(pdbs, fit=FALSE, full=FALSE, 
                       trim.inds=NULL, outpath = "pdbs_nma", ...) {

  if(class(pdbs)!="3dalign")
    stop("input 'pdbs' should be a list object as obtained from 'read.fasta.pdb'")

  #if(!is.null(xyz)) {
  #  if(class(xyz)=="matrix") {
  #    if(!identical(dim(pdbs$xyz), dim(xyz.core)))
  #      stop("dimension mismatch: input 'xyz' should have same dimensions as 'pdbs$xyz'")
  #  }
  #  else
  #    stop("input 'xyz' should be of type 'matrix'")
  #}
  
  if(!is.null(outpath))
    dir.create(outpath, FALSE)
  
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)
  
  f.inds <- NULL
  if(!is.null(trim.inds)) {
    strip <- TRUE
    ##f.inds$res <- intersect(gaps.res$f.inds, trim.inds)
    f.inds$res <- trim.inds
    f.inds$pos <- atom2xyz(f.inds$res)
  }
  else {
    strip <- FALSE
    ## indices only used for fit.xyz
    f.inds$res <- gaps.res$f.inds
    f.inds$pos <- gaps.pos$f.inds
  }

  keep <- 20
  if ((length(f.inds$res)*3) < (keep+6))
    keep <- length(f.inds$res)*3 - 6

  if(nrow(bounds(f.inds$res))>1 && strip)
    warning("NOTE: Truncated pdbs at non-terminus positions. \n\t - Fluctuations at neighboring positions may be affected.")
  
  if(fit) {
    #if(!is.null(xyz))
    #  warning("'pdbs$xyz' issued to re-fitting since 'fit=TRUE'")

    xyz <- fit.xyz(fixed = pdbs$xyz[1, ], mobile = pdbs,
                   fixed.inds = f.inds$pos, mobile.inds = f.inds$pos)
                   ##pdb.path = ".", pdbext = "", outpath = "core_fitlsq", full.pdbs = TRUE, het2atom = TRUE)
  }
  
  if(strip) {
    flucts <- matrix(NA, nrow=nrow(gaps.res$bin), ncol=length(f.inds$res))
  }
  else {
    flucts <- matrix(NA, nrow=nrow(gaps.res$bin), ncol=ncol(gaps.res$bin))
  }

  ## List object to store each modes object
  all.modes <- NULL

  ## 3D array- containing the modes vectors for each structure
  modes.array <- NULL
  
  if(full) {
    ## 3N x Num-modes x Num-structs
    modes.array <- array(NA, dim=c(length(f.inds$pos), keep, nrow(gaps.res$bin)))
    all.modes <- list()
  }
    
  if(is.null(outpath))
    fname <- tempfile(fileext = "pdb")

  ## Loop through each structure in 'pdbs'
  for ( i in 1:nrow(pdbs$xyz) ) {
    if(!strip) {
      f.inds <- NULL
      f.inds$res <- which(gaps.res$bin[i,]==0)
      f.inds$pos <- atom2xyz(f.inds$res)
    }
    
    resno <- pdbs$resno[i,f.inds$res]
    chain <- pdbs$chain[i,f.inds$res]
    resid <- aa123(pdbs$ali[i,f.inds$res])
    
    if(is.null(xyz))
      tmp.xyz <- pdbs$xyz[i, f.inds$pos]
    else
      tmp.xyz <- xyz[i, f.inds$pos]

    if(!is.null(outpath))
      fname <- file.path(outpath, basename(pdbs$id[i]))

    ## Make the PDB object
    write.pdb(pdb=NULL, xyz=tmp.xyz, resno=resno, chain=chain,
              resid=resid, file=fname)
    tmp.pdb <- read.pdb(fname)
    
    print(length(tmp.xyz))
    
    if(is.null(outpath))
       unlink(fname)

    ## Calculate normal modes
    modes <- nma(tmp.pdb, ...)

    if(full) {
      new.modes <- matrix(NA, ncol=keep, nrow=ncol(pdbs$xyz))
      
      j <- 1
      for(k in 7:(keep+6)) {
        new.modes[f.inds$pos, j] <- modes$U[,k]
        j <- j+1
      }
      
      ##print(dim(new.modes))
      
      all.modes[[i]] <- modes
      if(strip)
        modes.array[,,i] <- new.modes[f.inds$pos,]
      else
        modes.array[,,i] <- new.modes[gaps.pos$f.inds,]
      
    }
    
    if(!strip)
      flucts[i, f.inds$res] <- modes$fluctuations
    else
      flucts[i, ] <- modes$fluctuations
  }
  
  row.names(flucts) <- row.names(pdbs$xyz)
  out <- list(fluctuations=flucts, modes.array=modes.array, all.modes=all.modes)
      
  return(out)
}
