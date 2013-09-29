"view.modes" <-
  function(modes, mode=NULL, outprefix="mode_vecs",
           scale=5, dual=FALSE, launch=FALSE) {

    if(missing(modes))
      stop("modes object missing")
    
    if("nma" %in% class(modes)) {
      if(is.null(mode))
        mode <- 7
      xyz <- modes$xyz
      mode.vecs <- matrix(modes$modes[,mode], ncol=3, byrow=T)
    }

    else if (class(modes)=="pca") {
      if(is.null(mode))
        mode <- 1
      xyz <- modes$mean
      mode.vecs <- matrix(modes$U[,mode], ncol=3, byrow=T)
    }

    ## calc all vec lengths (for coloring later)
    all.lens <- apply(mode.vecs, 1, function(x) sqrt(sum(x**2)))
    
    ## make temp-files
    if(is.null(outprefix)) {
      pdbfile <- tempfile(fileext = ".inpcrd.pdb")
      outfile <- tempfile(fileext = ".py")
    }
    else {
      pdbfile <- paste(outprefix, ".inpcrd.pdb", sep="")
      outfile <- paste(outprefix, ".py", sep="")
    }

    ## start building pymol script
    scr <- c("from pymol import cmd")
    scr <- c(scr, "from pymol.cgo import *")
    scr <- c(scr, paste("cmd.load('", pdbfile, "', 'prot')", sep=""))
    scr <- c(scr, "cmd.show('cartoon')")
    scr <- c(scr, "cmd.set('cartoon_trace_atoms', 1)")
      
    ## define color range 
    blues <- colorRamp(c("blue", "white", "red"))

    ## Arrow widths
    w.body <- 0.15; w.head <- 0.2

    scr <- c(scr, "obj=[]")
    for ( i in 1:nrow(mode.vecs)) {
      inds    <- atom2xyz(i) 
      coords  <- xyz[inds]

      ## For coloring (longest vec has length=1)
      tmp.len <- sqrt(sum((mode.vecs[i,]/max(all.lens))**2))
      if(tmp.len>1)
        tmp.len <- 1

      col <- blues(tmp.len)
      col <- round(col/256,4)
      col <- paste(col, collapse=", ")

      ## Main vector
      tmp.vec <- mode.vecs[i,] * scale
      
      ## For arrow head
      if(sqrt(sum(tmp.vec**2))<1)
        norm.vec <-  tmp.vec
      else
        norm.vec <- normalize.vector(mode.vecs[i,])

      ## Set vectors
      arrow.vec.a <- (coords + tmp.vec)
      head.vec.a  <- (arrow.vec.a + (norm.vec))

      arrow.vec.b <- (coords - tmp.vec)
      head.vec.b <- (arrow.vec.b - (norm.vec))

      a  <- paste(coords,      collapse=",")
      b1 <- paste(arrow.vec.a, collapse=",")
      c1 <- paste(head.vec.a,  collapse=",")
      b2 <- paste(arrow.vec.b, collapse=",")
      c2 <- paste(head.vec.b,  collapse=",")


      ## Arrow body
      scr <- c(scr, paste("obj.extend([CYLINDER", a, b1, w.body, col, col, "])", sep=", "))
      if(dual)
        scr <- c(scr, paste("obj.extend([CYLINDER", a, b2, w.body, col, col, "])", sep=", "))
      
      ## Arrow heads
      scr <- c(scr, paste("obj.extend([CONE", b1, c1, w.head, 0.0, col, col, 1.0, 1.0,"])", sep=", "))
      if(dual)
        scr <- c(scr, paste("obj.extend([CONE", b2, c2, w.head, 0.0, col, col, 1.0, 1.0,"])", sep=", "))
    }
      
    name <- "vecs"
    scr <- c(scr, paste("cmd.load_cgo(obj, '", name, "')", sep=""))
    
    ## Write PDB structure file
    write.pdb(xyz=xyz, file=pdbfile)

    ## Write python script or PDB with conect records
    write(scr, file=outfile, sep="\n")
    
    if(launch) {
      ## Open pymol
      cmd <- paste('pymol', outfile)
      
      os1 <- .Platform$OS.type
      if (os1 == "windows") {
        shell(shQuote(cmd))
      }
      else {
        if(Sys.info()["sysname"]=="Darwin") {
          system(paste("open -a MacPyMOL", outfile))
        } else { system(cmd) }
      }
    }
  }
