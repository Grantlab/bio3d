trim <- function(s) {
  ##- Remove leading and trailing spaces from character strings
  s <- sub("^ +", "", s)
  s <- sub(" +$", "", s)
  s[(s=="")]<-""
  return(s)
}

randstr <- function() {
  return(basename(tempfile(pattern="")))
}


"pdbs2sse" <- function(pdbs, ind=1, rm.gaps=FALSE) {
  ind <- ind[1]
  if(file.exists(pdbs$id[ind]))
    id <- pdbs$id[ind]
  else if(file.exists(rownames(pdbs$ali)[ind]))
    id <- rownames(pdbs$ali)[ind]

  sse.aln <- NULL
  pdb.ref <- try(read.pdb(id), silent=TRUE)

  if(inherits(pdb.ref, "try-error"))
    pdb.ref <- try(read.pdb(substr(basename(id), 1, 4)), silent=TRUE)

  gaps.res <- gap.inspect(pdbs$ali)

  sse.ref <- NULL
  if(!inherits(pdb.ref, "try-error"))
    sse.ref <- try(dssp(pdb.ref), silent=TRUE)

  if(!inherits(sse.ref, "try-error") & !inherits(pdb.ref, "try-error")) {
    if(rm.gaps) {
      resid <- paste0(pdbs$resno[ind, gaps.res$f.inds], pdbs$chain[ind, gaps.res$f.inds])
    }
    else {
      resid <- paste0(pdbs$resno[ind, ], pdbs$chain[ind, ])
    }

    ## Helices
    resid.helix <- unbound(sse.ref$helix$start, sse.ref$helix$end)
    resid.helix <- paste0(resid.helix, rep(sse.ref$helix$chain, sse.ref$helix$length))
    inds        <- which(resid %in% resid.helix)

    ## inds points now to the position in the alignment where the helices are
    new.sse <- bounds( seq(1, length(resid))[inds] )
    if(length(new.sse) > 0) {
      sse.aln$helix$start  <- new.sse[,"start"]
      sse.aln$helix$end    <- new.sse[,"end"]
      sse.aln$helix$length <- new.sse[,"length"]
    }

    ## Sheets
    resid.sheet <- unbound(sse.ref$sheet$start, sse.ref$sheet$end)
    resid.sheet <- paste0(resid.sheet, rep(sse.ref$sheet$chain, sse.ref$sheet$length))
    inds        <- which(resid %in% resid.sheet)

    new.sse <- bounds( seq(1, length(resid))[inds] )
    if(length(new.sse) > 0) {
      sse.aln$sheet$start  <- new.sse[,"start"]
      sse.aln$sheet$end    <- new.sse[,"end"]
      sse.aln$sheet$length <- new.sse[,"length"]
    }

    ## SSE vector
    sse <- rep(" ", length(resid))
    for(i in 1:length(sse.aln$helix$start))
      sse[sse.aln$helix$start[i]:sse.aln$helix$end[i]] <- "H"

    for(i in 1:length(sse.aln$sheet$start))
      sse[sse.aln$sheet$start[i]:sse.aln$sheet$end[i]] <- "E"

    sse.aln$sse <- sse
  }
  else {
    msg <- NULL
    if(inherits(pdb.ref, "try-error"))
      msg = c(msg, paste("File not found:", pdbs$id[1]))
    if(inherits(sse.ref, "try-error"))
      msg = c(msg, "Launching external program 'DSSP' failed")

    warning(paste("SSE cannot be drawn", msg, sep="\n  "))
  }

  return(sse.aln)
}

'col2hex' <- function(cname) {
    colMat <- col2rgb(cname)
    rgb(red = colMat[1, ]/255,
        green = colMat[2, ]/255,
        blue = colMat[3, ]/255)
}

