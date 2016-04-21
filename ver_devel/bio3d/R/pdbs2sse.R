"pdbs2sse" <- function(pdbs, ind=NULL, rm.gaps=TRUE, resno=TRUE, pdb=FALSE, ...) {
  ## Log the call
  cl <- match.call()
  by.resno <- resno

  if(is.null(ind))
    ind <- 1:length(pdbs$id)

  gaps.res <- gap.inspect(pdbs$ali)

  ## Use SSE information from pdbs object
  if(!is.null(pdbs$sse) & !pdb) {
      message("Extracting SSE from pdbs$sse")

      if(rm.gaps) {
          sse <- pdbs$sse[ind, gaps.res$f.inds, drop=FALSE]
          resno <- pdbs$resno[ind[1], gaps.res$f.inds]
          chain <- pdbs$chain[ind[1], gaps.res$f.inds]
      }
      else {
          sse <- pdbs$sse[ind,, drop=FALSE]
          resno <- pdbs$resno[ind[1], ]
          chain <- pdbs$chain[ind[1], ]
      }

      if(nrow(sse) > 1) {
          h.inds <- which(apply(sse, 2, function(x) sum(x=="H")) == length(ind))
          e.inds <- which(apply(sse, 2, function(x) sum(x=="E")) == length(ind))
      }
      else {
          h.inds <- which(sse == "H")
          e.inds <- which(sse == "E")
      }

      if(by.resno) {
          h <- bounds( resno[h.inds], pre.sort=FALSE )
          e <- bounds( resno[e.inds], pre.sort=FALSE )
      }
      else {
          h <- bounds( h.inds, pre.sort=FALSE )
          e <- bounds( e.inds, pre.sort=FALSE )
      }

      out <- list()
      out$sse <- sse
      out$helix$start <- h[, "start"]
      out$helix$end <- h[, "end"]
      out$helix$length <- h[, "length"]
      out$helix$chain <- chain[ bounds(h.inds)[, "start"] ]

      out$sheet$start <- e[, "start"]
      out$sheet$end <- e[, "end"]
      out$sheet$length <- e[, "length"]
      out$sheet$chain <- chain[ bounds(e.inds)[, "start"] ]

      out$call <- cl
      class(out) <- "sse"
      return(out)
  }

  ind <- ind[1]
  message(paste("Re-reading PDB (", basename.pdb(pdbs$id[ind]), ") to extract SSE", sep=""))
    
  if(file.exists(pdbs$id[ind]))
    id <- pdbs$id[ind]

  sse.aln <- NULL
  pdb.ref <- try(read.pdb(id), silent=TRUE)

  if(inherits(pdb.ref, "try-error"))
    pdb.ref <- try(read.pdb(substr(basename(id), 1, 4)), silent=TRUE)

  sse.ref <- NULL
  if(!inherits(pdb.ref, "try-error"))
    sse.ref <- try(dssp(pdb.ref, ...), silent=TRUE)

  if(!inherits(sse.ref, "try-error") & !inherits(pdb.ref, "try-error")) {
    if(rm.gaps) {
      resno <- pdbs$resno[ind, gaps.res$f.inds]
      chain <- pdbs$chain[ind, gaps.res$f.inds]
    }
    else {
      resno <- pdbs$resno[ind, ]
      chain <- pdbs$chain[ind, ]
    }

    resid <- paste0(resno, chain)

    ## Helices
    if(length(sse.ref$helix$start) > 0) {
      resid.helix <- unbound(sse.ref$helix$start, sse.ref$helix$end)
      resid.helix <- paste0(resid.helix, rep(sse.ref$helix$chain, sse.ref$helix$length))
      inds        <- which(resid %in% resid.helix)

      ## inds points to the position in the alignment where the helices are
      if(by.resno) {
          resno.sse <- resno[ inds ]
          new.sse <- bounds(resno.sse, pre.sort=FALSE)
      }
      else {
          sids <- 1:length(resid)
          resno.sse <- sids[inds]
          new.sse <- bounds(resno.sse, pre.sort=FALSE)
      }
      chain.sse <- chain[ bounds(inds, pre.sort=FALSE)[, "start"] ]

      if(length(new.sse) > 0) {
        sse.aln$helix$start  <- new.sse[,"start"]
        sse.aln$helix$end    <- new.sse[,"end"]
        sse.aln$helix$length <- new.sse[,"length"]
        sse.aln$helix$chain <- chain.sse
      }
    }

    ## Sheets
    if(length(sse.ref$sheet$start) > 0) {
      resid.sheet <- unbound(sse.ref$sheet$start, sse.ref$sheet$end)
      resid.sheet <- paste0(resid.sheet, rep(sse.ref$sheet$chain, sse.ref$sheet$length))
      inds        <- which(resid %in% resid.sheet)

      if(by.resno) {
          resno.sse <- resno[ inds ]
          new.sse <- bounds(resno.sse, pre.sort=FALSE)
      }
      else {
          sids <- 1:length(resid)
          resno.sse <- sids[inds]
          new.sse <- bounds(resno.sse, pre.sort=FALSE)
      }
      chain.sse <- chain[ bounds(inds, pre.sort=FALSE)[, "start"] ]

      if(length(new.sse) > 0) {
        sse.aln$sheet$start  <- new.sse[,"start"]
        sse.aln$sheet$end    <- new.sse[,"end"]
        sse.aln$sheet$length <- new.sse[,"length"]
        sse.aln$sheet$chain <- chain.sse
      }
    }

    ## SSE vector
    sse <- rep(" ", length(resid))
    if(length(sse.aln$helix$start) > 0) {
      for(i in 1:length(sse.aln$helix$start))
        sse[sse.aln$helix$start[i]:sse.aln$helix$end[i]] <- "H"
    }

    if(length(sse.aln$sheet$start) > 0) {
      for(i in 1:length(sse.aln$sheet$start))
        sse[sse.aln$sheet$start[i]:sse.aln$sheet$end[i]] <- "E"
    }

    sse.aln$sse <- sse
  }
  else {
    msg <- NULL
    if(inherits(pdb.ref, "try-error"))
      msg = c(msg, paste("File not found:", pdbs$id[1]))
    if(inherits(sse.ref, "try-error"))
      msg = c(msg, "Launching external program 'DSSP' failed")

    warning(paste("SSE failed, ", msg, sep="\n  "))
  }

  sse.aln$call <- cl
  class(sse.aln) <- "sse"
  return(sse.aln)
}
