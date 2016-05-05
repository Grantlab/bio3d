#' @param ... (in \code{gnm.pdbs}) additional arguments passed to \code{gnm.pdb}. 
#' @inheritParams aanma.pdbs
#' @rdname gnm
gnm.pdbs <- function(x, fit=TRUE, full=FALSE, subspace=NULL, rm.gaps=TRUE,
                     gc.first=TRUE, ncore=NULL, ...) {

  pdbs <- x
 
  if(!inherits(pdbs, "pdbs"))
    stop("input 'pdbs' should be a list object as obtained from 'read.fasta.pdb'")

  ## Log the call
  cl <- match.call()
  
  ## Parallelized by parallel package
  ncore <- setup.ncore(ncore)

  if(ncore>1) {
    prev.warn <- getOption("warn")
    options(warn=1)
    on.exit(options(warn=prev.warn))
  }

  dots <- list(...)
  if('outmodes' %in% names(dots)) {
     warning('Customized "outmodes" is not supported')
     dots$outmodes <- NULL
  }
  if("keep" %in% names(dots))
    nm.keep <- dots$keep
  else
    nm.keep <- NULL
  dots$check.connect <- FALSE

  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  ## Check connectivity
  con <- inspect.connectivity(pdbs, cut=4.05)
  if(!all(con)) {
    warning(paste(paste(basename(pdbs$id[which(!con)]), collapse=", "),
                  "might have missing residue(s) in structure:\n",
                  "  Fluctuations at neighboring positions may be affected."))
  }

  ## Number of modes to store in U.subspace
  if(is.null(subspace)) {
    keep <- length(gaps.res$f.inds)-1
  }
  else {
    keep <- subspace
    if (length(gaps.res$f.inds) < (keep+1))
      keep <- length(gaps.res$f.inds)-1
  }
  if(!is.null(nm.keep) && keep > nm.keep)
    keep <- nm.keep

  if(fit) {
     xyz <- fit.xyz(pdbs$xyz[1, ], pdbs, gaps.pos$f.inds, gaps.pos$f.inds, ncore=ncore)
     pdbs$xyz[,] <- xyz
  }

  #### Prepare for NMA calculation ####      
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
    modes.array <- array(NA, dim=c(length(gaps.res$f.inds), keep, nrow(gaps.res$bin)))
  else
    modes.array <- array(NA, dim=c(ncol(pdbs$ali), keep, nrow(gaps.res$bin)))

  ## store eigenvalues of the first modes
  L.mat <- matrix(NA, ncol=keep, nrow=nrow(gaps.res$bin))

  ### Memory usage ###
  dims <- dim(modes.array)
  mem.usage <- sum(c(as.numeric(object.size(modes.array)),
                     as.numeric(object.size(L.mat)),
                     as.numeric(object.size(flucts)),
                     as.numeric(object.size(matrix(NA, ncol=dims[3], nrow=dims[3]))) ))*2

  if(full) {
    if(is.null(nm.keep))
      tmpncol <- dims[2]
    else
      tmpncol <- nm.keep

    size.mat <- object.size(matrix(0.00000001, ncol=tmpncol, nrow=dims[1]))
    size.vec <- object.size(vector(length=dims[1], 'numeric'))

    tot.size <- ((size.mat * 2) + (size.vec * 4)) * length(pdbs$id)
    mem.usage <- mem.usage+tot.size
  }
  mem.usage=round(mem.usage/1048600,1)

  #### Print overview of scheduled calcualtion ####
  cat("\nDetails of Scheduled Calculation:\n")
  cat(paste("  ...", length(pdbs$id), "input structures", "\n"))
  if(keep>0)
    cat(paste("  ...", "storing", keep, "eigenvectors for each structure", "\n"))
  if(keep>0)
    cat(paste("  ...", "dimension of x$U.subspace: (",
              paste(dims[1], dims[2], dims[3], sep="x"), ")\n"))

  if(fit)
    cat(paste("  ...", "coordinate superposition prior to NM calculation", "\n"))

  if(full)
    cat(paste("  ... individual complete 'gnm' objects will be stored", "\n"))

  if(rm.gaps)
    cat(paste("  ... aligned eigenvectors (gap containing positions removed) ", "\n"))

  if(mem.usage>0)
    cat(paste("  ...", "estimated memory usage of final 'eGNM' object:", mem.usage, "Mb \n"))

  cat("\n")

  ##### Start modes calculation #####

  ## Initialize progress bar
  pb <- .init.pb(ncore, min=0, max=length(pdbs$id))

  all.modes <- mclapply(1:length(pdbs$id), function(i) {

     if(gc.first) gc()

     pdb <- pdbs2pdb(pdbs, i, rm.gaps = FALSE)[[1]]
     sele <- match(gaps.res$f.inds, which(!is.gap(pdbs$ali[i, ])))
     sele <- list(atom=sele, xyz=atom2xyz(sele))
     class(sele) <- 'select'

     if(rm.gaps)
        modes <- try(do.call(gnm, c(list(x=pdb, outmodes=sele), dots)))
     else
        modes <- try(do.call(gnm, c(list(x=pdb), dots)))

     if(inherits(modes, 'try-error')) {
        .close.pb(ncore, pb)
        stop(paste('Encounter errors in ', i, 'th structure', sep=''))
     }

     .update.pb(ncore, pb, i)

     modes$call <- NULL
     return( modes )

  }, mc.cores=ncore)

  ## Finish progress bar
  .close.pb(ncore, pb)

  ##### Finalize calculation #####
  for(i in 1:length(all.modes)) {
     if(rm.gaps) {
       flucts[i, ] <- all.modes[[i]]$fluctuations
       modes.array[,,i] <- all.modes[[i]]$U[, 2:(keep+1)]
     } else {
       flucts[i, !is.gap(pdbs$ali[i, ])] <- all.modes[[i]]$fluctuations
       modes.array[!is.gap(pdbs$ali[i, ]),, i] <- all.modes[[i]]$U[, 2:(keep+1)]
     }
     L.mat[i, ] <- all.modes[[i]]$L[2:(keep+1)]
  }
  if(!full) all.modes <- NULL

  ##### RMSIP ######
  rmsip.map <- NULL
  if(rm.gaps) {
    rmsip.map <- .calcRMSIP(modes.array, ncore=ncore)
    rownames(rmsip.map) <- basename(rownames(pdbs$xyz))
    colnames(rmsip.map) <- basename(rownames(pdbs$xyz))

    if(!fit)
      warning("rmsip calculated on non-fitted structures:
               ensure that your input coordinates are pre-fitted.")
  }

  rownames(flucts) <- basename(rownames(pdbs$xyz))
  out <- list(fluctuations=flucts, rmsip=rmsip.map,
              U.subspace=modes.array, L=L.mat, full.nma=all.modes, call=cl)

  class(out) <- c("egnm", "enma")

  return(out)
}

