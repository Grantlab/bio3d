#' Ensemble Normal Mode Analysis with All-Atom ENM
#'
#' Perform normal mode analysis (NMA) on an ensemble of aligned protein
#' structures using all-atom elastic network model (aaENM).
#'
#' @details This function builds elastic network model (ENM) using all heavy 
#'    atoms and performs subsequent normal mode analysis (NMA) on a set of 
#'    aligned protein structures obtained with function \code{\link{read.all}}.
#'    The main purpose is to automate ensemble normal mode analysis using 
#'    all-atom ENMs.
#'
#'    By default, the effective Hessian for all C-alpha atoms is calculated 
#'    based on the Hessian built from all heavy atoms (including ligand atoms if 
#'    \code{ligand=TRUE}). Returned values include aligned mode vectors and 
#'    (when \code{full=TRUE}) a list containing the full \sQuote{nma} objects
#'    one per each structure. When \sQuote{rm.gaps=TRUE} the unaligned atoms 
#'    are ommited from output. With default arguments \sQuote{rmsip} provides 
#'    RMSIP values for all pairwise structures.
#'
#'    When \code{outmodes} is provided and is not \sQuote{calpha} 
#'    (e.g. \sQuote{noh}. See \code{\link{aanma}} for more details), the 
#'    function simply returns a list of \sQuote{nma} objects, one per each 
#'    structure, and no aligned mode vector is returned. In this case, the 
#'    arguments \code{full}, \code{subspace}, and \code{rm.gaps} are ignored. 
#'    This is equivalent to a wrapper function repeatedly calling 
#'    \code{\link{aanma}}.
#'
#' @param pdbs an \sQuote{pdbs} object as obtained from \code{\link{read.all}}. 
#' @param fit logical, if TRUE C-alpha coordinate based superposition is 
#'    performed prior to normal mode calculations. 
#' @param full logical, if TRUE return the complete, full structure,
#'    \sQuote{nma} objects.
#' @param subspace number of eigenvectors to store for further analysis.
#' @param rm.gaps logical, if TRUE obtain the hessian matrices for only
#'    atoms in the aligned positions (non-gap positions in all aligned
#'    structures). Thus, gap positions are removed from output.
#' @param ligand logical, if TRUE ligand molecules are also included in the 
#'    calculation.
#' @param outpath character string specifing the output directory to
#'    which the PDB structures should be written.
#' @param gc.first logical, if TRUE will call gc() first before mode calculation
#'    for each structure. This is to avoid memory overload when 
#'    \code{ncore > 1}.
#' @param ncore number of CPU cores used to do the calculation.
#' @param ... additional arguments to \code{\link{aanma}}.
#'
#' @return Returns a list of \sQuote{nma} objects (\code{outmodes} is provided 
#'    and is not \sQuote{calpha}) or an \sQuote{enma} object with the following 
#'    components:
#'    \item{fluctuations }{ a numeric matrix containing aligned atomic
#'      fluctuations with one row per input structure. }
#'    \item{rmsip }{ a numeric matrix of pair wise RMSIP values (only the ten
#'      lowest frequency modes are included in the calculation). }
#'    \item{U.subspace }{ a three-dimensional array with aligned
#'      eigenvectors  (corresponding to the subspace defined by the first N
#'      non-trivial eigenvectors (\sQuote{U}) of the \sQuote{nma} object). }
#'    \item{L }{ numeric matrix containing the raw eigenvalues with one row
#'      per input structure. }
#'    \item{full.nma }{ a list with a \code{nma} object for each input
#'      structure (available only when \code{full=TRUE}). }
#'
#' @seealso 
#'      For normal mode analysis on single structure PDB:
#'      \code{\link{aanma}}
#'
#'      For conventional C-alpha based normal mode analysis:
#'      \code{\link{nma}}, \code{\link{nma.pdbs}}.
#'
#'      For the analysis of the resulting \sQuote{eNMA} object:
#'      \code{\link{mktrj.enma}}, \code{\link{dccm.enma}},
#'      \code{\link{plot.enma}}, \code{\link{cov.enma}}.
#'    
#'      Similarity measures:
#'      \code{\link{sip}}, \code{\link{covsoverlap}},
#'      \code{\link{bhattacharyya}}, \code{\link{rmsip}}.
#'    
#'      Related functionality:
#'      \code{\link{read.all}}.
#'
#' @author Xin-Qiu Yao & Lars Skjaerven
#' 
#' @examples
#' \donttest{
#'   # Needs MUSCLE installed - testing excluded
#'   if(check.utility("muscle")) {
#'
#'     ## Fetch PDB files and split to chain A only PDB files
#'     ids <- c("1a70_A", "1czp_A", "1frd_A", "1fxi_A", "1iue_A", "1pfd_A")
#'     files <- get.pdb(ids, split = TRUE, path = tempdir())
#'     
#'     ## Sequence Alignement
#'     aln <- pdbaln(files, outfile = tempfile())
#'     
#'     ## Read all pdb coordinates
#'     pdbs <- read.all(aln)
#'
#'     ## Normal mode analysis on aligned data
#'     modes <- aanma(pdbs, rm.gaps=TRUE)
#'     
#'     ## Plot fluctuation data
#'     plot(modes, pdbs=pdbs)
#'     
#'     ## Cluster on Fluctuation similariy
#'     sip <- sip(modes)
#'     hc <- hclust(dist(sip))
#'     col <- cutree(hc, k=3)
#'     
#'     ## Plot fluctuation data
#'     plot(modes, pdbs=pdbs, col=col)
#'     
#'     ## RMSIP is pre-calculated
#'     heatmap(1-modes$rmsip)
#'     
#'     ## Bhattacharyya coefficient
#'     bc <- bhattacharyya(modes)
#'     heatmap(1-bc)
#'
#'   }
#' }
aanma.pdbs <- function(pdbs, fit=TRUE, full=FALSE, subspace=NULL, rm.gaps=TRUE, 
  ligand=FALSE, outpath=NULL, gc.first=TRUE, ncore=NULL, ...) {

  if(!inherits(pdbs, "pdbs") || is.null(pdbs$all))
    stop("input 'pdbs' should be a list object as obtained from 'read.all()'")

  ## Log the call
  cl <- match.call()

  if(!is.null(outpath))
    dir.create(outpath, FALSE)

  ## Parallelized by parallel package
  ncore = setup.ncore(ncore)
  if(ncore>1) {
    prev.warn <- getOption("warn")
    options(warn=1)
    on.exit(options(warn=prev.warn))
  }

  dots <- list(...)

  aligned.modes <- TRUE
  if('outmodes' %in% names(dots)) {
     if(dots$outmodes == 'noh') {
       warning(paste('Non C-alpha atoms are selected for output modes.', 
        'A plain list of "nma" objects will be returned.'))
       aligned.modes <- FALSE
     } else if(dots$outmodes == 'calpha') {
       ## use default select, i.e. 'non-gap' C-alpha positions.
       dots$outmodes <- NULL
     } else {
       stop('Unsupported "outmodes"')
     }
  }

  if("keep" %in% names(dots))
    nm.keep <- dots$keep
  else
    nm.keep <- NULL

  if('rtb' %in% names(dots)) {
    rtb <- dots$rtb
    if(isTRUE(rtb)) {
      if('nmer' %in% names(dots) && dots$nmer!=1)
        stop( paste('Currently only nmer=1 is supported for', 
         'RTB based ensemble normal mode analysis') )
    }
  } else {
    rtb <- formals(aanma.pdb)$rtb
  }

  if('reduced' %in% names(dots)) 
    reduced <- dots$reduced
  else
    reduced <- formals(aanma.pdb)$reduced

  gaps.res <- gap.inspect(pdbs$ali)
  gaps.pos <- gap.inspect(pdbs$xyz)

  ## Number of modes to store in U.subspace
  if(is.null(subspace)) {
    keep <- length(gaps.pos$f.inds)-6
  }
  else {
    keep <- subspace
    if (length(gaps.pos$f.inds) < (keep+6))
      keep <- length(gaps.pos$f.inds)-6
  }
  if(!is.null(nm.keep) && keep > nm.keep) 
    keep <- nm.keep

  ## Convert from pdbs to pdb
  all.pdb <- pdbs2pdb(pdbs, all.atom=TRUE, ncore=ncore)

  ## Check if some structures are unavailable
  keep.inds <- which(sapply(all.pdb, is.pdb))
  if(length(keep.inds) == 0) 
     stop('No pdb file is found.')
  if(length(keep.inds) != length(all.pdb)) {
     warning(paste('Following pdb coordinates are not found: ', pdbs$id[-keep.inds], 
               sep=''))

     pdbs <- trim(pdbs, row.inds=keep.inds, col.inds=1:ncol(pdbs$ali))
     all.pdb <- all.pdb[keep.inds]
  }

  ## Remove ligand
  if(!ligand) 
    all.pdb <- mclapply(all.pdb, trim, 'protein', mc.cores=ncore)

  ## Non-gap position for each pdb
  nogap.inds <- mclapply(1:length(all.pdb), function(i)  {
     pdb <- all.pdb[[i]]
     pdb2aln.ind(pdbs, pdb, aln.id = pdbs$id[i], gaps.res$f.inds, file=NULL)$b
  }, mc.cores=ncore)

  if(fit) {
     cat('Fitting pdb structures')
     all.pdb <- mclapply(1:length(all.pdb), function(i) {
        cat('.')
        pdb <- all.pdb[[i]]
        xyz <- fit.xyz(pdbs$xyz[1, ], pdb$xyz, gaps.pos$f.inds, nogap.inds[[i]]$xyz)
        if(!is.null(outpath))
           ofile <- file.path(outpath, basename(pdbs$id[i]))
        else
           ofile <- tempfile()
        write.pdb(pdb, xyz=xyz, file=ofile)
        read.pdb(ofile)
     }, mc.cores=ncore, mc.allow.recursive=FALSE)
     cat('done\n')

  } else {

     if(!is.null(outpath)) 
        mclapply(1:length(all.pdb), function(i)
           write.pdb(all.pdb[[i]], file=file.path(outpath, basename(pdbs$id[i]))), 
           mc.cores = ncore)
  }

  if(aligned.modes) {

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
      modes.array <- array(NA, dim=c(length(gaps.pos$f.inds), keep, nrow(gaps.res$bin)))
    else
      modes.array <- array(NA, dim=c(ncol(pdbs$xyz), keep, nrow(gaps.res$bin)))

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
      cat(paste("  ... individual complete 'nma' objects will be stored", "\n"))

    if(rm.gaps)
      cat(paste("  ... aligned eigenvectors (gap containing positions removed) ", "\n"))

    if(reduced)
      cat(paste("  ... reduced all-atom ENM will be employed", "\n"))

    if(rtb)
      cat(paste("  ... rotation-translation block (RTB) approximation will be applied", "\n"))

    if(mem.usage>0)
      cat(paste("  ...", "estimated memory usage of final 'eNMA' object:", mem.usage, "Mb \n"))

    cat("\n")

  }

  ##### Start modes calculation #####
  ## Initialize progress bar
  pb <- .init.pb(ncore, min=0, max=length(pdbs$id))

  all.modes <- mclapply(1:length(all.pdb), function(i) {

     if(gc.first) gc()

     pdb <- all.pdb[[i]]
     nogap.inds <- nogap.inds[[i]]

     if(aligned.modes) {
 
       if(rm.gaps) 
          capture.output( modes <- try(do.call(aanma, c(list(pdb=pdb, 
              outmodes=nogap.inds), dots))) )
       else
          capture.output( modes <- try(do.call(aanma, c(list(pdb=pdb), dots))) )

     } else {

       capture.output( modes <- try(do.call(aanma, c(list(pdb=pdb), dots))) )

     } 

     if(inherits(modes, 'try-error')) {
        .close.pb(pb)
        stop(paste('Encounter errors in ', i, 'th structure', sep=''))
     }

     .update.pb(pb)

     modes$call <- NULL
     return( modes )

  }, mc.cores=ncore)

  ## Finish progress bar
  .close.pb(pb)

  if(!aligned.modes) return( all.modes )

  ##### Finalize calculation #####
  for(i in 1:length(all.modes)) {
     if(rm.gaps) {
       flucts[i, ] <- all.modes[[i]]$fluctuations
       modes.array[,,i] <- all.modes[[i]]$U[, 7:(keep+6)]
     } else {
       flucts[i, !is.gap(pdbs$ali[i, ])] <- all.modes[[i]]$fluctuations
       modes.array[!is.gap(pdbs$xyz[i, ]),, i] <- all.modes[[i]]$U[, 7:(keep+6)]
     }
     L.mat[i, ] <- all.modes[[i]]$L[7:(keep+6)]
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

  if(fit) {
    xyz <- fit.xyz(fixed = pdbs$xyz[1, ], mobile = pdbs,
                   fixed.inds = gaps.pos$f.inds, mobile.inds = gaps.pos$f.inds,
                   ncore = ncore)
  }
  else {
     xyz <- pdbs$xyz
  }

  rownames(flucts) <- basename(rownames(pdbs$xyz))
  out <- list(fluctuations=flucts, rmsip=rmsip.map,
              U.subspace=modes.array, L=L.mat, full.nma=all.modes, 
              xyz=xyz, call=cl)

  class(out) <- "enma"

  return(out)
}

.init.pb <- function(ncore, min=0, max=1) {
  if(ncore == 1) {

     return ( txtProgressBar(min=min, max=max, style=3) )

  } else if(ncore > 1) {

     mcparallel <- get("mcparallel", envir = getNamespace("parallel"))
     mccollect <- get("mccollect", envir = getNamespace("parallel"))

     fpb <- fifo(tempfile(), open = "w+b", blocking = T)

     # spawn a child process for message printing
     child <- mcparallel({
        pb <- txtProgressBar(min=min, max=max, style=3)
        progress <- 0.0
        while(progress < max && !isIncomplete(fpb)) {
           msg <- readBin(fpb, "double")
           progress <- progress + as.numeric(msg)
           setTxtProgressBar(pb, progress)
        }
        close(pb)
     } )

     names(fpb) <- child$pid
     return(fpb)
  }
}
.update.pb <- function(pb, step=1) {

  if(inherits(pb, "txtProgressBar")) {
     i <- getTxtProgressBar(pb)
     setTxtProgressBar(pb, i+step)
  }
  else {
     if(inherits(pb, "fifo"))
        writeBin(step, pb)
  }

}
.close.pb <- function(pb) {
   if(inherits(pb, "fifo")) {
      mccollect <- get("mccollect", envir = getNamespace("parallel"))
      mccollect(as.numeric(names(pb)))
#      mccollect(as.numeric(names(pb)))
   }
   close(pb)
}

