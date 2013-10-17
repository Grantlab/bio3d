## Two options - both calculating modes on the FULL structure:
## 1 -  use k <- kaa - ((kaq %*% kqq.inv) %*% kqa) to derive hessian for core atoms
## 2 - return the full objects

"nma.pdbs" <- function(pdbs, fit=TRUE, full=FALSE, 
                       rm.gaps=TRUE, outpath = "pdbs_nma", ...) {

  if(class(pdbs)!="3dalign")
    stop("input 'pdbs' should be a list object as obtained from 'read.fasta.pdb'")
  
  if(!is.null(outpath))
    dir.create(outpath, FALSE)

  ## Passing arguments to functions aa2mass and nma
  am.names <- names(formals( aa2mass ))
  nm.names <- names(formals( nma ))
  
  dots <- list(...)
  am.args <- dots[names(dots) %in% am.names]
  nm.args <- dots[names(dots) %in% nm.names]

  ## Limiting input 
  if("mass" %in% names(nm.args))
    mass <- nm.args$mass
  else
    mass <- TRUE
  if("ff" %in% names(nm.args))
    ff <- nm.args$ff
  else
    ff <- 'calpha'
  if("temp" %in% names(nm.args))
    temp <- nm.args$temp
  else
    temp <- 300
  if("keep" %in% names(nm.args))
    nm.keep <- nm.args$temp
  else
    nm.keep <- NULL
  
  if(!all((names(nm.args) %in% c("mass", "ff", "temp", "keep")))) {
    war <- paste(names(nm.args)[! names(nm.args) %in% c("mass", "ff", "temp", "keep") ], collapse=", ")
    warning(paste("ignoring arguments:", war))
  }
  
  ## Set indicies
  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  ## Use for later indexing
  pdbs$inds <- matrix(NA, ncol=ncol(pdbs$resno), nrow=nrow(pdbs$resno))
  
  ## Number of modes to keep  
  keep <- 20
  if (length(gaps.pos$f.inds) < (keep+6))
    keep <- length(gaps.pos$f.inds)*3 - 6

  ## Coordiantes - fit or not
  if(fit) {
    xyz <- fit.xyz(fixed = pdbs$xyz[1, ], mobile = pdbs,
                   fixed.inds = gaps.pos$f.inds, mobile.inds = gaps.pos$f.inds)
                   ##pdb.path = ".", pdbext = "", outpath = "core_fitlsq", full.pdbs = TRUE, het2atom = TRUE)
  }
  else
    xyz <- pdbs$xyz
  
  ## Fluctuations for each structure
  if(rm.gaps)
    flucts <- matrix(NA, nrow=nrow(gaps.res$bin), ncol=length(gaps.res$f.inds))
  else
    flucts <- matrix(NA, nrow=nrow(gaps.res$bin), ncol=ncol(gaps.res$bin))
  
  ## List object to store each modes object
  if(full)
    all.modes <- list()
  else
    all.modes <- NULL

  ## 3D array- containing the modes vectors for each structure
  if(rm.gaps)
    modes.array <- array(NA, dim=c(length(gaps.pos$f.inds), keep, nrow(gaps.res$bin)))
  else
    modes.array <- array(NA, dim=c(length(pdbs$xyz), keep, nrow(gaps.res$bin)))
  
  if(is.null(outpath))
    fname <- tempfile(fileext = "pdb")
  
  pb <- txtProgressBar(min=0, max=nrow(pdbs$xyz), style=3)

  ## Loop through each structure in 'pdbs'
  for ( i in 1:nrow(pdbs$xyz) ) {

    ## Set indices for this structure only
    f.inds <- NULL
    f.inds$res <- which(gaps.res$bin[i,]==0)
    f.inds$pos <- atom2xyz(f.inds$res)

    ## similar to $resno but sequential indices
    pdbs$inds[i, f.inds$res] <- seq(1, length(f.inds$res))
    
    ## Indices to extract from Hessian
    inds.inc <- pdbs$inds[i, gaps.res$f.inds]
    inds.exc <- pdbs$inds[i, gaps.res$t.inds][ !is.na(pdbs$inds[i, gaps.res$t.inds]) ]

    inds.inc.xyz <- atom2xyz(inds.inc)
    inds.exc.xyz <- atom2xyz(inds.exc)

    ## Generate content of PDB object
    tmp.xyz <- xyz[i, f.inds$pos]
    resno   <- pdbs$resno[i,f.inds$res]
    chain   <- pdbs$chain[i,f.inds$res]

    ## Fix for missing chain IDs
    chain[is.na(chain)] <- ""
    
    ## Check if protein is 'complete'
    if(nrow(bounds(as.numeric(pdbs$resno[i,])))>1)
      warning(paste(basename(pdbs$id[i]), "might have missing residue(s) in structure:\n",
                    "  Fluctuations at neighboring positions may be affected."))
    
    ## Workaround for unknown residue types
    if(any(pdbs$ali[i,]=="X") && mass==TRUE) {
      if(!file.exists(pdbs$id[i])) {
        cat("\n")
        stop(paste("Non-standard residue type found in", basename(pdbs$id[i]), "\n",
                   "  attempt to re-read PDB file failed."))
      }
      
      pdb   <- read.pdb(pdbs$id[i])
      resid <- pdb$atom[atom.select(pdb, 'calpha', verbose=FALSE)$atom, "resid"]
    }
    else {
      resid <- aa123(pdbs$ali[i,f.inds$res])
    }
    
    ## Make a PDB object by writing to disk and re-read file
    if(!is.null(outpath))
      fname <- file.path(outpath, basename(pdbs$id[i]))

    write.pdb(pdb=NULL, xyz=tmp.xyz, resno=resno, chain=chain,
              resid=resid, file=fname)
    tmp.pdb <- read.pdb(fname)
        
    if(is.null(outpath))
      unlink(fname)

    if(mass) {
      masses <- try(
                    do.call('aa2mass', c(list(pdb=tmp.pdb, inds=NULL), am.args)),
                    silent=TRUE
                    )
      
      if(inherits(masses, "try-error")) {
        hmm <- attr(masses,"condition")
        cat("\n\n")
        stop(paste(hmm$message, "in file", basename(pdbs$id[i])))
      }
    }
    else
      masses <- NULL
    
    pfc.fun <- load.enmff(ff)
    
    if(rm.gaps) {
      ## Build the hessian of the complete structure
      hess    <- build.hessian(tmp.pdb$xyz, pfc.fun, aa.mass=masses)

      ## Effective hessian for atoms in the aligned core
      if(length(inds.exc)>0) {
        kaa     <- hess[inds.inc.xyz, inds.inc.xyz]
        kqq.inv <- solve(hess[inds.exc.xyz, inds.exc.xyz])
        kaq     <- hess[inds.inc.xyz, inds.exc.xyz]
        kqa     <- t(kaq)
        k <- kaa - ((kaq %*% kqq.inv) %*% kqa)
      }
      else {
        k <- hess
      }
      
      ## Second PDB - containing only the aligned atoms
      tmp2.pdb <- NULL
      tmp2.pdb$atom <- tmp.pdb$atom[inds.inc,]
      tmp2.pdb$xyz <-  tmp.pdb$xyz[inds.inc.xyz]
      tmp2.pdb$calpha <- (tmp2.pdb$atom[,"elety"]=="CA") && (tmp2.pdb$atom[,"resid"]!="CA")
      class(tmp2.pdb) <- "pdb"

      if(mass)
        masses <- do.call('aa2mass', c(list(pdb=tmp2.pdb, inds=NULL), am.args))
      else
        masses <- NULL

      ## Calculate the modes
      invisible(capture.output( modes <- nma(tmp2.pdb, ff=ff, mass=mass, temp=temp, keep=nm.keep, hessian=k, aa.mass=masses)))
    }
    else {
      ## Calculate the modes
      invisible(capture.output( modes <- nma(pdb=tmp.pdb, ff=ff, mass=mass, temp=temp, keep=nm.keep, aa.mass=masses)))
    }

    if(rm.gaps)
      modes.mat <- matrix(NA, ncol=keep, nrow=nrow(modes$U))
    else
      modes.mat <- matrix(NA, ncol=keep, nrow=ncol(pdbs$xyz))
    
    j <- 1
    for(k in 7:(keep+6)) {
      if(rm.gaps)
        modes.mat[, j] <- modes$U[,k]
      else
        modes.mat[f.inds$pos, j] <- modes$U[,k]
      j <- j+1
    }

    if(full)
      all.modes[[i]] <- modes
    modes.array[,,i] <- modes.mat
   

    if(rm.gaps)
      flucts[i, ] <- modes$fluctuations
    else
      flucts[i, f.inds$res] <- modes$fluctuations

    setTxtProgressBar(pb, i)
  }
  close(pb)

  calc.rmsip <- function(x) {
    n <- dim(x)[3]
    mat <- matrix(NA, n, n)
    inds <- rbind(pairwise(n),
                  matrix(rep(1:n,each=2), ncol=2, byrow=T))
    
    for(i in 1:nrow(inds)) { 
      mat[inds[i,1], inds[i,2]] <- rmsip(x[,,inds[i,1]],
                                         x[,,inds[i,2]])$rmsip
    }
    mat[ inds[,c(2,1)] ] = mat[ inds ]
    ##diag(mat) <- 1 ## better to calculate it
    return(round(mat, 4))
  }

  rmsip.map <- NULL
  if(rm.gaps) {
    rmsip.map <- calc.rmsip(modes.array)
    rownames(rmsip.map) <- basename(rownames(pdbs$xyz))
    colnames(rmsip.map) <- basename(rownames(pdbs$xyz))

    if(!fit)
      warning("rmsip calculated on non-fitted structures:
               ensure that your input coordinates are pre-fitted.")
  }
  
  rownames(flucts) <- basename(rownames(pdbs$xyz))
  out <- list(fluctuations=flucts, rmsip=rmsip.map,
              U.subs=modes.array, full.nma=all.modes )
      
  return(out)
}
