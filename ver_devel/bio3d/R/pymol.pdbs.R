.resno2str <- function(res, sep=c("+", "-")) {
  res <- res[!is.na(res)]
  if(!length(res)>0){
    return(NULL)
  }
  else {
    res1 <- bounds(res)
    res2 <- paste(res1[,"start"], res1[,"end"], sep=sep[2])
    inds <- res1[,"start"] == res1[,"end"]
    res2[inds] <- res1[inds, "start"]
    res3 <- paste(res2, collapse=sep[1])
    return(res3)
  }
}

pymol <- function(...)
  UseMethod("pymol")

pymol.pdbs <- function(pdbs, col=NULL, file=NULL,
                       type="script", exefile = "pymol") {
  
  allowed <- c("session", "script", "launch")
  if(!type %in% allowed) {
    stop(paste("input argument 'type' must be either of:",
               paste(allowed, collapse=", ")))
  }

  ## output file name
  if(is.null(file)) {
    if(type=="session")
      file <- "R.pse"
    if(type=="script")
      file <- "R.pml"
  }
  
  ## Check if the program is executable
  if(type %in% c("session", "launch")) {
    ver <- "-cq"
    os1 <- .Platform$OS.type
    status <- system(paste(exefile, ver),
                     ignore.stderr = TRUE, ignore.stdout = TRUE)
    
    if(!(status %in% c(0,1)))
      stop(paste("Launching external program failed\n",
                 "  make sure '", exefile, "' is in your search path", sep=""))
  }

  ## use temp-dir unless we output a PML script
  if(type %in% c("session", "launch"))
    tdir <- tempdir()
  else
    tdir <- "."

  pdbdir <- paste(tdir, "pdbs", sep="/")
  if(!file.exists(pdbdir))
    dir.create(pdbdir)
  
  pmlfile <- tempfile(tmpdir=tdir, fileext=".pml")
  psefile <- tempfile(tmpdir=tdir, fileext=".pse")
  ids <- basename.pdb(pdbs$id)

  ## include stuff in the b-factor column
  bf <- NULL
  if(!is.null(col)) {
    if(col[1] == "rmsf" | col[1] == "putty") {
      bf <- rmsf(pdbs$xyz)
    }
    if(col[1] == "index2") {
      bf <- 1:ncol(pdbs$ali)/ncol(pdbs$ali)
    }
  }
  
  ## use all all-atom PDBs if they exist
  if(all(file.exists(pdbs$id))) {
    allatom <- TRUE
    files <- pdbs$id

    ## align all-atom PDBs to pdbs$xyz
    for(i in 1:length(pdbs$id)) {
      pdb <- read.pdb(files[i])
      sele <- atom.select(pdb, "calpha")
      gaps <- is.gap(pdbs$xyz[i,])
      pdb$xyz <- fit.xyz(pdbs$xyz[i, !gaps], pdb$xyz,
                         fixed.inds = 1:length(pdbs$xyz[i, !gaps]),
                         mobile.inds = sele$xyz)
      fn <- paste0(pdbdir, "/", ids[i], ".pdb")

      ## store new b-factor column to PDB
      tmpbf <- NULL
      if(!is.null(bf)) {
        gaps <- is.gap(pdbs$ali[i,])
        tmpbf <- pdb$atom$b*0
        tmpbf[sele$atom] <- bf[!gaps]
      }

      write.pdb(pdb, b=tmpbf, file=fn)
      files[i] <- fn
    }
  }
  else {
    ## use pdbs$xyz to build CA-atom PDBs
    allatom <- FALSE
    files <- rep(NA, length(pdbs$id))
    for(i in 1:length(pdbs$id)) {
      pdb <- pdbs2pdb(pdbs, inds=i)[[1]]
      fn <- paste0(pdbdir, "/", ids[i], ".pdb")

      ## store new b-factor column to PDB
      tmpbf <- NULL
      if(!is.null(bf)) {
        gaps <- is.gap(pdbs$ali[i,])
        tmpbf <- bf[!gaps]
      }
        
      write.pdb(pdb=pdb, b=tmpbf, file=fn)
      files[i] <- fn
    }
  }

  ## load PDBs
  lines <- rep(NA, 5*length(pdbs$id))
  for(i in 1:length(files)) {
    lines[i] <- paste("load", files[i])
  }

  ## line pointer
  l <- i
  
  ## Coloring
  if(!is.null(col)) {
    if(inherits(col, "core")) {
      core <- col
      l <- l+1
      lines[l] <- "color grey50"
      
      for(j in 1:length(files)) {
        res <- .resno2str(pdbs$resno[j, core$atom])
        if(!is.null(res)) {
          selname <- paste0(ids[j], "-core")
          lines[l+1] <- paste0("select ", selname, ", ", ids[j], " and resi ", res)
          lines[l+2] <- paste0("color red, ", selname)
          l <- l+2
        }
      }
    }
    
    if(col[1] == "gaps") {
      l <- l+1
      lines[l] <- "color grey50"
      
      gaps <- gap.inspect(pdbs$ali)
      for(j in 1:length(files)) {
        res <- .resno2str(pdbs$resno[j, gaps$t.inds])
        if(!is.null(res)) {
          selname <- paste0(ids[j], "-gap")
          lines[l+1] <- paste0("select ", selname, ", ", ids[j], " and resi ", res)
          lines[l+2] <- paste0("color red, ", selname)
          l <- l+2
        }
      }
    }
    
    if(length(col) > 1 & is.vector(col)) {
      if(length(col) != length(files))
        stop("col must be a vector with length equal to the number of structures in input pdbs")

      ## add more colors here
      cols <- c("black", "red", "green", "blue", "cyan", "purple", "yellow", "grey50")
      for(j in 1:length(files)) {
        lines[l+1] <- paste0("color ", cols[col[j]], ", ", ids[j])
        l <- l+1
      }
    }

    if(col[1] == "putty") {
      lines[l+1] <- "cartoon putty"
      lines[l+2] <- "as cartoon"
      lines[l+3] <- "unset cartoon_smooth_loops"
      lines[l+4] <- "unset cartoon_flat_sheets"
      lines[l+5] <- "spectrum b, rainbow"
      lines[l+6] <- "set cartoon_putty_radius, 0.2"
      l <- l+6
    }
    
    if(col[1] == "rmsf") {
      l <- l+1
      lines[l] <- "spectrum b, rainbow"
    }

    if(col[1] == "index") {
      for(i in 1:length(pdbs$id)) {
        l <- l+1
        lines[l] <- paste("spectrum count, rainbow,", ids[i], "and name C*")
      }
    }

    if(col[1] == "index2") {
      for(i in 1:length(pdbs$id)) {
        l <- l+1
        lines[l] <- paste("spectrum b, rainbow,", ids[i])
      }
    }

    
  } ## coloring ends
  
  if(!allatom) {
    lines[l+1] <- "as cartoon"
    lines[l+2] <- "set cartoon_trace_atoms, 1"
    l <- l+2
  }

  if(type == "session")
    lines[l+1] <- paste("save", psefile)
  lines <- lines[!is.na(lines)]
  print(lines)
  write.table(lines, file=pmlfile, append=FALSE, quote=FALSE, sep="\n",
              row.names=FALSE, col.names=FALSE)

  if(type %in% c("session", "launch")) {
    if(type == "session")
      args <- "-cq"
    else
      args <- ""
    
    ## Open pymol
    cmd <- paste('pymol', args, pmlfile)
    
    os1 <- .Platform$OS.type
    if (os1 == "windows") {
      success <- shell(shQuote(cmd))
    }
    else {
      if(Sys.info()["sysname"]=="Darwin") {
        success <- system(paste("open -a MacPyMOL", pmlfile))
      }
      else {
        success <- system(cmd)
      }
    }
    
    if(success!=0)
      stop(paste("An error occurred while running command\n '",
                 exefile, "'", sep=""))
  }

  if(type == "session") {
    file.copy(psefile, file, overwrite=TRUE)
    unlink(pmlfile)
    unlink(psefile)
    message(paste("PyMOL session written to file", file))
    invisible(file)
  }

  if(type == "script") {
    file.copy(pmlfile, file, overwrite=TRUE)
    unlink(pmlfile)
    message(paste("PyMOL script written to file", file))
    invisible(file)
  }
}


