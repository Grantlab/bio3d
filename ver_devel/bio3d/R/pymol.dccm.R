pymol.dccm <- function(dccm, pdb, file=NULL,
                       step=0.2, omit=0.2, radius = 0.15,
                       type="script", exefile = "pymol", ...) {
  
  allowed <- c("session", "script", "launch", "pdb")
  if(!type %in% allowed) {
    stop(paste("input argument 'type' must be either of:",
               paste(allowed, collapse=", ")))
  }

  if(type %in% "pdb")
    pymol <- FALSE
  else
    pymol <- TRUE

  ## Check if the program is executable
    if(type %in% c("session", "launch")) {
        
        ## determine path to exefile
        exefile1 <- .get.exepath(exefile)
        
        ## Check if the program is executable
        success <- .test.exefile(exefile1)
        
        if(!success) {
            stop(paste("Launching external program failed\n",
                       "  make sure '", exefile, "' is in your search path", sep=""))
        }
        exefile <- exefile1
    }
  
  if(is.pdb(pdb)) {
    ca.inds <- atom.select(pdb, 'calpha', verbose=FALSE)
    bb.inds <- atom.select(pdb, 'backbone', verbose=FALSE)
    xyz <- pdb$xyz[ca.inds$xyz]
    
    ## If more than CA atoms are provided, assume its enough to draw cartoon in pymol
    if(length(pdb$xyz[bb.inds$xyz])==length(xyz))
      ca.pdb <- TRUE
    else
      ca.pdb <- FALSE
  }
  else {
    xyz <- pdb
  }
  
  if(missing(dccm))
    stop("correlation matrix must be provided")
  if(missing(xyz))
    stop("cartesian coordinates missing")

  if(is.matrix(radius)) {
    if(!all(dim(dccm) == dim(radius)))
      stop("dimension mismatch. provide a 'radius' matrix with the same dimensions as 'dccm'")
  }    
  
  dims <- dim(dccm)
  if((length(xyz)/3)!=dims[1L])
    stop("unequal vector lengths")

  if(is.null(file)) {
    if(type=="session")
      file <- "R.pse"
    if(type=="script")
      file <- "R.py"
    if(type=="pdb")
      file <- "R.pdb"
  }

  ## use temp-dir unless we output a PML script
  if(type %in% c("session", "launch"))
    tdir <- tempdir()
  else
    tdir <- "."

  psefile <- tempfile(tmpdir=tdir, fileext=".pse")
  pyfile <- tempfile(tmpdir=tdir, fileext=".py")
  pdbfile <- tempfile(tmpdir=tdir, fileext=".pdb")
  
  ## Build the new PDB or pymol script in a vector
  scr <- c()
    
  if(pymol) {
    ## start pymol script
    scr <- c("from pymol import cmd")
    scr <- c(scr, "from pymol.cgo import *")
    scr <- c(scr, paste("cmd.load('", 
      normalizePath(pdbfile, winslash='/', mustWork=FALSE), 
       "', 'prot')", sep=""))
    scr <- c(scr, "cmd.show('cartoon')")
    
    if(!is.pdb(pdb) || ca.pdb)
      scr <- c(scr, "cmd.set('cartoon_trace_atoms', 1)")
    
    ## define color range 
    blues <- colorRamp(c("white", "blue"))
    reds  <- colorRamp(c("white", "red"))
    w <- radius 
  }
  else {
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
   
    ## make the positive and negative distributions symmetric
    if(lower<0)
      sele <- intersect( which(dccm>=lower), which(dccm<upper) )
    else
      sele <- intersect( which(dccm>lower), which(dccm<=upper) )

    if(length(sele)==0)
      next
    
    f       <- matrix(FALSE, ncol(dccm), nrow(dccm))
    f[sele] <- TRUE
    inds    <- which(f, arr.ind=TRUE)
    
    if(pymol) {
      scr <- c(scr, "obj=[]")
    }
    else {
      m <- m+1
      chain <- LETTERS[m]
    }
    
    for ( j in 1:nrow(inds) ) {
      x <- inds[j,1]; y <- inds[j,2]; 
      
      if(x==y)
        next
      
      val <- dccm[x,y]           ## corr coeff
      k   <- atom2xyz(inds[j,1]) ## resi 1
      l   <- atom2xyz(inds[j,2]) ## resi 2

      if(is.matrix(w))
        w2 <- abs(w[x, y])
      else
        w2 <- w
      
      if(pymol) {
        a <- paste(xyz[k], collapse=",")
        b <- paste(xyz[l], collapse=",")
        
        if(val<=0)
          col <- blues( abs(val) )
        else
          col <- reds( abs(val) )
        
        col <- round(col/256,4)
        col <- paste(col, collapse=", ")
        
        str <- paste("obj.extend([CYLINDER", a, b, w2, col, col, "])", sep=", ")
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
    
    if(pymol) {
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

  if(type == "session") {
    scr <- c(scr, paste0("cmd.save('", 
      normalizePath(psefile, winslash='/', mustWork=FALSE),
      "')"))
    }
  
  ## Write python script or PDB with conect records
  if(pymol) {
    write(scr, file=pyfile, sep="\n")
    
    ## Write PDB structure file
    if(is.pdb(pdb))
      write.pdb(pdb, file=pdbfile)
    else
      write.pdb(xyz=xyz, file=pdbfile)
  }
  else {
    write(scr, file=pdbfile, sep="\n")
  }
  
  if(type %in% c("session", "launch")) {
    if(type == "session")
      args <- "-cq"
    else
      args <- ""
    
    ## Open pymol
    cmd <- paste(exefile, args, pyfile)

    os1 <- Sys.info()["sysname"]
    if (os1 == "Windows") {
        status <- shell(paste(shQuote(exefile), args, pyfile))
    }
    else {
        status <- system(cmd)
    }
    
    if(!(status %in% c(0,1))) {
        stop(paste("An error occurred while running command\n '",
                   exefile, "'", sep=""))
    }

  }

  if(type == "session") {
    file.copy(psefile, file, overwrite=TRUE)
    message(paste("PyMOL session written to file", file))
    invisible(file)
  }
  
  if(type == "script") {
    file.copy(pyfile, file, overwrite=TRUE)
    unlink(pyfile)
    message(paste("PyMOL script written to file", file))
    invisible(file)
  }

  if(type == "pdb") {
    file.copy(pdbfile, file, overwrite=TRUE)
    unlink(pdbfile)
    message(paste("PDB written to file", file))
    invisible(file)
  }
  
}
