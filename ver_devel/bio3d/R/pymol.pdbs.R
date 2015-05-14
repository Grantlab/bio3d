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

pymol.pdbs <- function(pdbs, col=NULL, file="R.pse") {
  
  tdir <- tempdir()
  pmlfile <- tempfile(tmpdir=tdir, fileext=".pml")
  if(is.null(file))
    psefile <- tempfile(tmpdir=tdir, fileext=".pse")
  else
    psefile <- file
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
      fn <- paste0(tdir, "/", ids[i], ".pdb")

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
      fn <- paste0(tdir, "/", ids[i], ".pdb")

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
  
  lines[l+1] <- paste("save", psefile)
  lines <- lines[!is.na(lines)]
  write.table(lines, file=pmlfile, append=FALSE, quote=FALSE, sep="\n",
              row.names=FALSE, col.names=FALSE)
  system(paste("pymol -cq", pmlfile))
  
  message(paste("PyMOL session written to file", psefile))
  invisible(psefile)
}


