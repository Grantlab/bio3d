"plot.enma" <-
  function(x,
           pdbs=NULL,
           xlab=NULL, ylab="Fluctuations",
           ...) {

    if(!(inherits(x, "enma")))
        stop("provide a enma object as obtained from 'nma.pdbs'")

      yval <- x$fluctuations


      ## SSE information
      sse.aln <- NULL
      dots <- list(...)

      ## reference structure for SSE and resno
      if("ind" %in% names(dots)) {
          ref.ind <- dots[["ind"]]
          dots[["ind"]] <- NULL

          ## Note - plot.fluct use plotb3(yval[1,]) for base plotting
          ## Lines 54-61 checks for NAs in yval[1,] to trim SSE
          ## SSE must therefore correspond to yval[1,]
          ## ref.ind != 1 might give wrong sse annotation in plot
          
          warning("reference structure can not be set")
          ref.ind <- 1
      }
      else {
          ref.ind <- 1
      }

      ## use first non-NA in col as ref.ind
      if("col" %in% names(dots)) {
          col <- dots[["col"]]

          if(length(col) != nrow(yval))
              stop("length of col doesn't match dimension of x")
          
          if(any(is.na(col))) {
              ref.ind <- which(!is.na(col))[1]
          }
      }

  
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
        resno <- pdbs$resno[ref.ind, ]

        if(is.null(xlab)) {
            xlab <- paste0('Residue number (reference PDB: ',
                           basename.pdb(pdbs$id[ref.ind]), ')')
        }
    }
    else {
        resno <- NULL

        if(is.null(xlab))
            xlab <- "Alignment Position"
    }

      
    if( "sse" %in% names(dots) ) {
        sse.aln <- dots$sse
        dots$sse <- NULL
    }
      
    if(!is.null(pdbs) & is.null(sse.aln)) {
        sse.aln <- pdbs2sse(pdbs, ind=ref.ind, rm.gaps=rm.gaps, resno=TRUE)
        sse.aln$sse[ is.na(yval[ref.ind,]) ] <- NA
        ## see lines 54-61 in plotb3.R
    }

    if( "rm.gaps" %in% names(dots) ) {
        warning("'rm.gaps=TRUE' might result in incorrect SSE annotation")
    }

    ## Plot fluctuations
    do.call('plot.fluct', c(list(x=yval,
                                 resno=resno, sse=sse.aln,
                                 xlab=xlab, ylab=ylab),
                            dots))

}

