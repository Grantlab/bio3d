"plot.enma" <-
  function(x, pdbs=NULL, entropy=FALSE, col=NULL, signif=FALSE,
           pcut=0.005, qcut=0.04,
           xlab="Residue Position", ylab="Fluctuations",
           mar = c(4, 5, 2, 2), 
           ...) {
    
    if(!inherits(x, "enma"))
      stop("provide a enma object as obtained from 'nma.pdbs'")

    if(is.null(x$call$rm.gaps))
      rm.gaps <- TRUE
    else if(x$call$rm.gaps=="T" || x$call$rm.gaps=="TRUE")
      rm.gaps <- TRUE
    else
      rm.gaps <- FALSE
    
    dots <- list(...)
    sse.aln <- NULL
    if(!is.null(pdbs)) {
      if( "sse" %in% names(dots) )
        warning(paste("incompatible arguments: SSE information from 'pdbs'\n", 
                      "  will not be generated when 'sse' is provided"))

      pdb.ref <- try(read.pdb(pdbs$id[1]), silent=TRUE)
      if(inherits(pdb.ref, "try-error")) {
        ## Try more
        pdb.ref <- try(read.pdb(substr(basename(pdbs$id[1]), 1, 4)), silent=TRUE)
      }
      if(!inherits(pdb.ref, "try-error")) {
        ## Reference PDB and SSE
        
        sse.ref <- dssp(pdb.ref)
        
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
        sse.aln$helix$start <- new.sse[,"start"]
        sse.aln$helix$end <- new.sse[,"end"]
        
        ## Sheets
        resno.sheet <- unbound(sse.ref$sheet$start, sse.ref$sheet$end)
        inds <- which(resnos %in% as.character(resno.sheet))
        
        new.sse <- bounds( seq(1, length(resnos))[inds] )
        sse.aln$sheet$start <- new.sse[,"start"]
        sse.aln$sheet$end <- new.sse[,"end"]
      }
      else {
        warning(paste("can not draw SSE info. \n\tFile not found:", pdbs$id[1]))
      }
    }
    
    if( !"sse" %in% names(dots) ) {
      dots$sse <- sse.aln
    }
    
    if(is.null(col))
      col <- seq(1, nrow(x$fluctuations))

    op <- par(no.readonly=TRUE)
    on.exit(par(op))
    if(entropy)
      par(mfrow=c(3,1), mar=mar)
    else
      par(mar=mar)
    
    if(signif) {
      # Do student's t-test
      ns <- levels(as.factor(col))
      if(length(ns) !=2) {
        warning("Number of states is not equal to 2. Ignore plot of significance")
      } else {
         inds1 <- which(col==ns[1])
         inds2 <- which(col==ns[2])
         p <- NULL; q <- NULL
         for(i in 1:ncol(x$fluctuations)) {
           p <- c(p, t.test(x$fluctuations[inds1,i],
                         x$fluctuations[inds2,i],
                         alternative="two.sided")$p.value)
           m <- mean(x$fluctuations[inds1,i])
           n <- mean(x$fluctuations[inds2,i])
           q <- c(q, abs(m-n))
         }
         sig <- which(p<pcut & q>qcut)
      
         # Plot significance as shaded blocks
         maxy <- max(x$fluctuations, na.rm=TRUE)
         do.call('plot.bio3d', c(list(x=x$fluctuations[1,], xlab=xlab, ylab=ylab,
                                      ylim=c(0,maxy),
                                      col=col[1]), type="n", dots))
         bds <- bounds(sig)
         ii <- 1:nrow(bds) 
         rect(bds[ii,1], rep(0, length(ii)), bds[ii,2], 
              rep(maxy, length(ii)), 
              col=rep("lightblue", length(ii)), border=NA)
#         ii <- which(bds[, 3]<=1)
#         lines(bds[ii, 1], rep(maxy, length(ii)), type="h", col="lightblue")
         lines( x$fluctuations[1,], col=col[1], type="h" )
      }
    } else {
    
       ## Plot fluctuations plot
       do.call('plot.bio3d', c(list(x=x$fluctuations[1,], xlab=xlab, ylab=ylab,
                                    ylim=c(0,max(x$fluctuations, na.rm=TRUE)),
                                    col=col[1]), dots))
    }
    for(i in 2:nrow(x$fluctuations)) {
      lines( x$fluctuations[i,], col=col[i] )
    }


    ## Entropy and fluct.variance
    if (entropy) {
      ## Fluctuation variance
      fluct.sd <- apply(x$fluctuations, 2, var, na.rm=T)
      
      ## Sequence Entropy
      if(rm.gaps) {
        h   <- entropy(pdbs$ali[,gaps.res$f.inds])
      }
      else {
        h   <- entropy(pdbs)
      }
      H <- h$H.norm
      
      mp <- barplot(fluct.sd, ylab="Fluct. variance")
      axis(side=1, at=mp[ seq(1, nrow(mp), by=50) ],
           labels=seq(0,length(H),by=50))
      box()
      
      
      mp <- barplot(H, border=NA, ylab = "Seq. entropy",
                    ylim=c(0,1))
      axis(side=2, at=c(0.2,0.4, 0.6, 0.8))
      axis(side=1, at=mp[ seq(1, nrow(mp), by=50) ],
           labels=seq(0,length(H),by=50))
      box()
      
    }
    
  }
