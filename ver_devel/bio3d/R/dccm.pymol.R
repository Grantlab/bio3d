"dccm.pymol" <-
  function(dccm, xyz, step=0.2, omit=0.2) {
    
    if(missing(dccm))
      stop("correlation matrix must be provided")
    if(missing(xyz))
      stop("cartesian coordinates missing")

    dims <- dim(dccm)
    if((length(xyz)/3)!=dims[1L])
      stop("unequal vector lengths")

    ## make temp-files
    pdbfile <- tempfile(fileext = ".pdb")
    outfile <- tempfile(fileext = ".py")

    ## mask lower tri of correlation matrix
    dccm[lower.tri(dccm, diag=TRUE)] <- NA

    ## start script
    scr <- c("from pymol import cmd")
    scr <- c(scr, "from pymol.cgo import *")
    scr <- c(scr, paste("cmd.load('", pdbfile, "', 'prot')", sep=""))
    scr <- c(scr, "cmd.show('cartoon')")
    scr <- c(scr, "cmd.set('cartoon_trace_atoms', 1)")

    ## define color range
    blues <- colorRamp(c("white", "blue"))
    reds <- colorRamp(c("white", "red"))
    w <- 0.15
    
    lims <- c(-1, 1)
    intervals <- seq(lims[1], lims[2], by=step)

    ## get rid of interval around 0
    if(!is.null(omit)) {
      i <- which(intervals>(omit-0.001))
      j <- which(intervals<(-omit+0.001))
      inds <- sort(c(i,j))
      intervals <- sort(intervals[inds])
    }
  
    for ( i in 1:(length(intervals)-1) ) {
      lower <- intervals[i]
      upper <- intervals[i+1]
      
      if(lower<0 && upper>0)
        next
      
      ##cat(lower, "\t", upper, "\n")
    
      sele <- intersect( which(dccm>lower), which(dccm<=upper) )
      if(length(sele)==0)
        next
      
      f <- matrix(FALSE, ncol(dccm), nrow(dccm))
      f[sele] <- TRUE
      inds <- which(f, arr.ind=TRUE)
      
      scr <- c(scr, "obj=[]")
      
      for ( j in 1:nrow(inds) ) {
        x <- inds[j,1]; y <- inds[j,2]; 
        
        if(x==y)
          next
        
        k <- atom2xyz(inds[j,1]) ## resi 1
        l <- atom2xyz(inds[j,2]) ## resi 2
      
        a <- paste(xyz[k], collapse=",")
        b <- paste(xyz[l], collapse=",")
        
        val <- dccm[x,y]
        if(val<=0)
          col <- blues( abs(val) )
        else
          col <- reds( abs(val) )
        
        col <- round(col/256,4)
        col <- paste(col, collapse=", ")
        
        str <- paste("obj.extend([CYLINDER", a, b, w, col, col, "])", sep=", ")
        scr <- c(scr, str)
      }
      ##name <- paste("obj", i, collapse="_")
      tmpa <- gsub("\\.", "", as.character(lower))
      tmpb <- gsub("\\.", "", as.character(upper))
      name <- paste("cor_", tmpa, "_", tmpb, sep="")
      
      str <- paste("cmd.load_cgo(obj, '", name, "')", sep="")
      scr <- c(scr, str)
    }
    
    ## write pymol script file and PDB file
    write.pdb(xyz=xyz, file=pdbfile)
    write(scr, file=outfile, sep="\n")

    ## Open pymol
    cmd <- paste('pymol', outfile)
    system(cmd)
  }
