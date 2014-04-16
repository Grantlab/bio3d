"plot.enma" <-
  function(x, y="fluctuations",
           pdbs=NULL, entropy=FALSE, variance=FALSE,
           col=NULL, signif=FALSE,
           pcut=0.005, qcut=0.04,
           xlab="Residue Position", ylab="Fluctuations",
           xlim=NULL, ylim=NULL,
           mar = c(4, 5, 2, 2),
           ...) {

    if(!inherits(x, "enma"))
      stop("provide a enma object as obtained from 'nma.pdbs'")

    if(is.null(pdbs) && entropy) {
      entropy=FALSE
      warning("forcing 'entropy=FALSE': entropy plot requires the corresponding 'pdbs' object")
    }

    if(is.null(x$call$rm.gaps))
      rm.gaps <- TRUE
    else if(x$call$rm.gaps=="T" || x$call$rm.gaps=="TRUE")
      rm.gaps <- TRUE
    else
      rm.gaps <- FALSE

    if(y=="deformations")
      yval <- x$deform
    else
      yval <- x$fluctuations

    if(!is.null(col) && any(is.na(col)))
      inds.plot <- which(!is.na(col))
    else
      inds.plot <- 1:nrow(yval)

    if(is.null(col))
      col <- seq(1, nrow(yval))

    if(is.null(ylim))
      ylim <- c(0,max(yval, na.rm=TRUE))
    
    if(is.null(xlim))
      xlim <- c(0,ncol(yval))

    dots <- list(...)
    sse.aln <- NULL
    if(!is.null(pdbs)) {
      if( "sse" %in% names(dots) )
        warning(paste("Incompatible arguments: SSE information from 'pdbs'\n",
                      "  will not be generated when 'sse' is provided"))

      pdb.ref <- try(read.pdb(pdbs$id[1]), silent=TRUE)
      if(inherits(pdb.ref, "try-error"))
        pdb.ref <- try(read.pdb(substr(basename(pdbs$id[1]), 1, 4)), silent=TRUE)

      sse.ref <- NULL
      if(!inherits(pdb.ref, "try-error"))
        sse.ref <- try(dssp(pdb.ref), silent=TRUE)

      if(!inherits(sse.ref, "try-error") && !inherits(pdb.ref, "try-error")) {
        if(rm.gaps) {
          gaps.res <- gap.inspect(pdbs$ali)
          resnos <- pdbs$resno[1, gaps.res$f.inds]
        }
        else {
          resnos <- pdbs$resno[1, ]
        }

        ## Helices
        resno.helix <- unbound(sse.ref$helix$start, sse.ref$helix$end)
        inds <- which(resnos %in% as.character(resno.helix))

        ## inds points now to the position in the alignment where the helices are
        new.sse <- bounds( seq(1, length(resnos))[inds] )
        if(length(new.sse) > 0) {
           sse.aln$helix$start <- new.sse[,"start"]
           sse.aln$helix$end <- new.sse[,"end"]
        }

        ## Sheets
        resno.sheet <- unbound(sse.ref$sheet$start, sse.ref$sheet$end)
        inds <- which(resnos %in% as.character(resno.sheet))

        new.sse <- bounds( seq(1, length(resnos))[inds] )
        if(length(new.sse) > 0) {
           sse.aln$sheet$start <- new.sse[,"start"]
           sse.aln$sheet$end <- new.sse[,"end"]
        }
      }
      else {
          msg <- NULL
          if(inherits(pdb.ref, "try-error"))
              msg = c(msg, paste("File not found:", pdbs$id[1]))
          if(inherits(sse.ref, "try-error"))
              msg = c(msg, "Launching external program 'DSSP' failed")

          warning(paste("SSE cannot be drawn", msg, sep="\n  "))
      }
    }

    if( !"sse" %in% names(dots) ) {
      dots$sse <- sse.aln
    }

    ## Perform test of significance
    ns <- levels(as.factor(col))
    if((length(ns) !=2) && signif) {
      warning("Number of states is not equal to 2. Ignoring significance test")
      signif <- FALSE
    }

    sig <- NULL
    if(signif) {
      inds1 <- which(col==ns[1])
      inds2 <- which(col==ns[2])

      if(length(inds1)>1 & length(inds2)>1) {
        p <- NULL; q <- NULL
        for(i in 1:ncol(yval)) {
          p <- c(p, t.test(yval[inds1,i],
                           yval[inds2,i],
                           alternative="two.sided")$p.value)
          m <- mean(yval[inds1,i])
          n <- mean(yval[inds2,i])
          q <- c(q, abs(m-n))
        }
        sig <- which(p<pcut & q>qcut)
      }

      ## Plot significance as shaded blocks
      if(is.null(sig))
        warning("Too few data points. Ignoring significance test")

      if(length(sig)==0) {
        ##warning("No significant differences found")
        sig <- NULL
      }
    }


    ## Configure plot
    nrows <- 1
    if(entropy)
      nrows=nrows+1
    if(variance)
      nrows=nrows+1

    op <- par(no.readonly=TRUE)
    on.exit(par(op))

    if(nrows>1)
      par(mfrow=c(nrows,1), mar=mar)
    else
      par(mar=mar)

    plot.new()
    plot.window(xlim=xlim, ylim=ylim, ...)
  
    ## If significance test was performed successfully
    if(!is.null(sig)) {
      ##maxy <- max(yval, na.rm=TRUE)
      bds <- bounds(sig)
      ii <- 1:nrow(bds)
      rect(bds[ii,1], rep(0, length(ii)), bds[ii,2],
           rep(ylim[2], length(ii)),
           col=rep("lightblue", length(ii)), border=NA)
    }

    ## Plot fluctuations / deformations
    par(new=TRUE)
    do.call('plot.bio3d', c(list(x=yval[inds.plot[1],], xlab=xlab, ylab=ylab,
                                 ylim=ylim, xlim=xlim, col=col[1]), type='h',
                            dots))
    
    ## Plot all lines (col==NA will not be plotted)
    for(i in 1:nrow(yval)) {
      lines( yval[i,], col=col[i], lwd=2, ... )
    }

    ## Fluctuation / deformations variance
    if (variance) {
      fluct.sd <- apply(yval, 2, var, na.rm=T)
      do.call('plot.bio3d', c(list(x=fluct.sd,
                                   xlab="Residue position",
                                   ylab="Fluct. variance",
                                   ##ylim=ylim,
                                   xlim=xlim,
                                   col=1), dots))
    }

    ## Sequence Entropy
    if (entropy) {
      if(rm.gaps) {
        h   <- entropy(pdbs$ali[,gaps.res$f.inds])
      }
      else {
        h   <- entropy(pdbs)
      }
      ##H <- h$H.10.norm
      H <- h$H.norm

      do.call('plot.bio3d', c(list(x=H,
                                   ylab="Seq. entropy",
                                   xlab="Residue position",
                                   col=1), dots))
    }


    out <- list(signif=sig, sse=sse.aln)
    invisible(out)                
  }
