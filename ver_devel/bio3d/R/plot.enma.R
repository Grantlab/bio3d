"plot.enma" <-
  function(x,
           pdbs=NULL,
           xlab=NULL, ylab="Fluctuations",
           ...) {

    if(!(inherits(x, "enma")))
        stop("provide a enma object as obtained from 'nma.pdbs'")
    yval <- x$fluctuations

    ## check for gaps
    gaps <- gap.inspect(yval)
    if(any(gaps$col>0))
      rm.gaps <- FALSE
    else
      rm.gaps <- TRUE


    ## check if pdbs match enma object
    gaps.pdbs <- NULL
    if(!is.null(pdbs)) {
      if(!inherits(pdbs, "pdbs")) {
        warning("argument 'pdbs' is not a 'pdbs' object (as obtained from pdbaln())")
        pdbs <- NULL
      }
      else {
        gaps.pdbs <- gap.inspect(pdbs$ali)

        if(rm.gaps)
          dims.pdbs <- dim(pdbs$ali[, gaps.pdbs$f.inds, drop=FALSE])
        else
          dims.pdbs <- dim(pdbs$ali)

        if(!identical(dim(yval), dims.pdbs)) {
          warning("dimenension mismatch between modes and pdbs object")
          pdbs <- NULL
        }
      }
    }

    ## trim pdbs object if rm.gaps=TRUE
    if(!is.null(pdbs)) {
        if(rm.gaps)
            pdbs <- trim(pdbs, col.inds=gaps.pdbs$f.inds)
    }

    if(!is.null(pdbs)) {
        resno <- pdbs$resno[1, ]

        if(is.null(xlab)) {
            xlab <- paste0('Residue number (reference PDB: ',
                           basename.pdb(pdbs$id[1]), ')')
        }
    }
    else {
        resno <- NULL

        if(is.null(xlab))
            xlab <- "Alignment Position"
    }

    ## SSE information
    sse.aln <- NULL
    dots <- list(...)
    if( "sse" %in% names(dots) ) {
      sse.aln <- dots$sse
      dots$sse <- NULL
    }

    if(!is.null(pdbs) & is.null(sse.aln)) {
        rm.gaps <- TRUE ## required by plotb3, see lines 70-76
        sse.aln <- pdbs2sse(pdbs, ind=1, rm.gaps=rm.gaps, resno=TRUE)
        print(sse.aln)
    }

    if( "rm.gaps" %in% names(dots) ) {
        warning("'rm.gaps=TRUE' might result in incorrect SSE annotation")
    }

    resno2 = resno[!is.na(resno)] ## req by plotb3 line 25

    ## Plot fluctuations
    do.call('plot.fluct', c(list(x=yval,
                                 resno=resno2, sse=sse.aln,
                                 xlab=xlab, ylab=ylab),
                            dots))

}

