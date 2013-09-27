"view.dccm" <-
  function(dccm, pdb, step=0.2, omit=0.2, type="pymol",
           outprefix="corr") {

    if(class(pdb)=="pdb") {
      ca.inds <- atom.select(pdb, 'calpha', verbose=FALSE)
      xyz <- pdb$xyz[ca.inds$xyz]
    }
    else {
      xyz <- pdb
    }
    
    if(missing(dccm))
      stop("correlation matrix must be provided")
    if(missing(xyz))
      stop("cartesian coordinates missing")

    if(type!="pdb" && type!="pymol")
      stop("provide type 'pdb' or 'pymol'")

    dims <- dim(dccm)
    if((length(xyz)/3)!=dims[1L])
      stop("unequal vector lengths")
    
    ## make temp-files
    if(is.null(outprefix)) {
      pdbfile <- tempfile(fileext = ".inpcrd.pdb")
      if(type=="pymol")
        outfile <- tempfile(fileext = ".py")
      else
        outfile <- tempfile(fileext = ".pdb")
    }
    else {
      pdbfile <- paste(outprefix, ".inpcrd.pdb", sep="")
      outfile <- paste(outprefix, ".pdb", sep="")
    }

    ## Build the new PDB or pymol script in a vector
    scr <- c()
    
    if(type=="pymol") {
      ## start pymol script
      scr <- c("from pymol import cmd")
      scr <- c(scr, "from pymol.cgo import *")
      scr <- c(scr, paste("cmd.load('", pdbfile, "', 'prot')", sep=""))
      scr <- c(scr, "cmd.show('cartoon')")

      if(class(pdb)!="pdb")
        scr <- c(scr, "cmd.set('cartoon_trace_atoms', 1)")
      
      ## define color range 
      blues <- colorRamp(c("white", "blue"))
      reds  <- colorRamp(c("white", "red"))
      w <- 0.15
    }
    else {
      chains <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      m <- 0
    }

    ## mask lower tri of correlation matrix
    dccm[lower.tri(dccm, diag=TRUE)] <- NA
    
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
      
      sele <- intersect( which(dccm>lower), which(dccm<=upper) )
      if(length(sele)==0)
        next
      
      f       <- matrix(FALSE, ncol(dccm), nrow(dccm))
      f[sele] <- TRUE
      inds    <- which(f, arr.ind=TRUE)

      if(type=="pymol") {
        scr <- c(scr, "obj=[]")
      }
      else {
        m <- m+1 
        chain <- substr(chains, m,m)
      }
      
      for ( j in 1:nrow(inds) ) {
        x <- inds[j,1]; y <- inds[j,2]; 
        
        if(x==y)
          next

        val <- dccm[x,y]         ## corr coeff
        k   <- atom2xyz(inds[j,1]) ## resi 1
        l   <- atom2xyz(inds[j,2]) ## resi 2
        
        if(type=="pymol") {
          a <- paste(xyz[k], collapse=",")
          b <- paste(xyz[l], collapse=",")
          
          if(val<=0)
            col <- blues( abs(val) )
          else
            col <- reds( abs(val) )
        
          col <- round(col/256,4)
          col <- paste(col, collapse=", ")

          str <- paste("obj.extend([CYLINDER", a, b, w, col, col, "])", sep=", ")
          scr <- c(scr, str)
        }
        else {
          a <- paste(format(xyz[k], justify="right", width=8), collapse="")
          b <- paste(format(xyz[l], justify="right", width=8), collapse="")
      
          val <- format(round(val,2), justify="right", width=6)
          res.str <- format(x, justify="right", width=4)
          str <- paste("ATOM   ", res.str, "  CA  ALA ", chain, res.str, "    ", a, "  0.00", val, sep="")
          scr <- c(scr, str)
          
          res.str <- format(y, justify="right", width=4)
          str <- paste("ATOM   ", res.str, "  CA  ALA ", chain, res.str, "    ", b, "  0.00", val,  sep="")
          scr <- c(scr, str)
          
          str <- paste("CONECT",
                       format(x, justify="right", width=5),
                       format(y, justify="right", width=5), sep="")
          scr <- c(scr, str)
        }
      }
      
      if(type=="pymol") {
        tmpa <- gsub("\\.", "", as.character(lower))
        tmpb <- gsub("\\.", "", as.character(upper))
        name <- paste("cor_", tmpa, "_", tmpb, sep="")
        str <- paste("cmd.load_cgo(obj, '", name, "')", sep="")
        scr <- c(scr, str)
      }
      else {
        str <- "TER"
        scr <- c(scr, str)
      }
    }
    
    ## write pymol script file and PDB file
    if(class(pdb)=="pdb")
      write.pdb(pdb, file=pdbfile)
    else
      write.pdb(xyz=xyz, file=pdbfile)
    write(scr, file=outfile, sep="\n")

    ## Open pymol
    cmd <- paste('pymol', outfile)

    os1 <- .Platform$OS.type
    if (os1 == "windows") {
      shell(shQuote(cmd))
    }
    else {
      system(cmd)
    }
  }
